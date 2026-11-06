#!/usr/bin/env python3
"""Finite NEW matched joins/local-field normalization oracle; no material run.

SymPy differentiates an independently assembled real-space scalar. SciPy
propagates that finite dense system. This does not replay ClampedArray's suite
or establish material/density/size/nuclear convergence.
"""
from pathlib import Path
import argparse
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
for directory in (ROOT/'src',HERE,HERE.parent/'0025',HERE.parent/'0024'):
    sys.path.insert(0,str(directory))

import numpy as np
import sympy as sp
from scipy.linalg import expm, sqrtm
from matched_clamped_case import build_clamped
from matched_clamped_isotopes import isotope_masses, nist_hydrogen_masses
from matched_clamped_readouts import connected_readouts
from matched_clamped_sources import source_bindings
from run_full_field_moments_material import pack
from substrate_framework.material_modes import AMU_RY

TOL = 1e-10


def fixture():
    onsite = np.array([[4.,.13,-.08],[.13,6.,.11],[-.08,.11,8.]])
    offsite = np.array([[.2,-.07,.03],[-.07,.3,.02],[.03,.02,.4]])
    return dict(common_force_lattice_ry_bohr2=np.stack((onsite,offsite)).reshape(2,1,1,1,3,1,3),
        atomic_masses_amu=np.array([1.008]),atom_names=np.array(['H']),atom_indices=np.array([9]),
        valence_charges=np.array([1.3]),cell_bohr=np.diag([2.,3.,4.]),
        contact=np.array([[.7,.04,-.02],[.04,.9,.03],[-.02,.03,1.1]]),
        mixed_ionic_field=np.array([[.12,-.09,.02],[.04,.17,-.06],[-.03,.08,.14]]),
        ionic_force_ry_bohr=np.array([[.3,-.2,.1]]),
        positive_gradient_ry_bohr=np.array([[.07,.11,-.04]]),
        field_gradient=np.array([.06,-.03,.02]))


def real_scalar(arrays, atomic_mass, independent):
    """Direct coordinates: (u0,p0,u1,p1), then each (a,Pi); no q formulas."""
    nf = 2 if independent else 1
    size = 12+6*nf
    x = sp.symbols('x:'+str(size))
    mass = atomic_mass*AMU_RY
    charge = float(arrays['valence_charges'][0])
    volume = abs(np.linalg.det(arrays['cell_bohr']))
    field_mass = volume*(1 if independent else 2)/(8*np.pi)
    force = (arrays['ionic_force_ry_bohr']-arrays['positive_gradient_ry_bohr']).ravel()
    c = arrays['common_force_lattice_ry_bohr2'].reshape(2,3,3)
    energy = 0
    for r in range(2):
        a = x[12+6*(r if independent else 0):15+6*(r if independent else 0)]
        u,p = x[6*r:6*r+3],x[6*r+3:6*r+6]
        energy += sum((p[i]-charge*a[i])**2/(2*mass)-force[i]*u[i] for i in range(3))
        energy += sum(u[i]*arrays['mixed_ionic_field'][i,j]*a[j] for i in range(3) for j in range(3))
        energy += sum(arrays['field_gradient'][i]*a[i] for i in range(3))
        energy += sum(a[i]*arrays['contact'][i,j]*a[j]/2 for i in range(3) for j in range(3))
        for s in range(2):
            if not independent or r == s:
                energy += sum(u[i]*c[(r-s)%2,i,j]*x[6*s+j]/2 for i in range(3) for j in range(3))
    for f in range(nf):
        energy += sum(x[15+6*f+i]**2/(2*field_mass) for i in range(3))
    zero = dict.fromkeys(x,0)
    h = np.asarray(sp.hessian(energy,x),float)
    linear = np.asarray([sp.diff(energy,v).subs(zero) for v in x],float)
    j,w = np.zeros((size,size)),np.zeros((size,size),complex)
    widths = np.sqrt(mass*np.diag(c[0]))
    omega = np.sqrt(np.diag(c[0])/mass)
    for r in range(2):
        for i in range(3):
            u,p = 6*r+i,6*r+3+i
            j[u,p],j[p,u] = 1,-1
            w[u,u],w[p,p],w[u,p],w[p,u] = 1/(2*widths[i]),widths[i]/2,.5j,-.5j
    stiffness = (1 if independent else 2)*(arrays['contact']+np.eye(3)*charge**2/mass)
    root = np.asarray(sqrtm(stiffness),float)
    for f in range(nf):
        a,pi = slice(12+6*f,15+6*f),slice(15+6*f,18+6*f)
        j[a,pi],j[pi,a] = np.eye(3),-np.eye(3)
        w[a,a] = np.linalg.inv(root)/(2*np.sqrt(field_mass))
        w[pi,pi] = root*np.sqrt(field_mass)/2
        w[a,pi],w[pi,a] = .5j*np.eye(3),-.5j*np.eye(3)
    return h,linear,j,w,omega,field_mass


def evolve(h,linear,j,time):
    size = len(h)
    augmented = np.zeros((size+1,size+1))
    augmented[:size,:size] = j@h
    augmented[:size,size] = j@linear
    value = expm(time*augmented)
    return value[:size,:size],value[:size,size]


def checks_for_case(arrays, variant, control, close, killed):
    prefix = variant+'/'+control+'/'
    masses,mass_record = isotope_masses(arrays,variant)
    expected_mass = {'original':1.008,'H1':1.00782503223,'H2':2.01410177812}[variant]
    close(prefix+'original-or-isotope-mass',masses,[expected_mass])
    if mass_record['original_H_exactly_NIST_H1']:
        raise AssertionError('standard atomic weight mislabeled H1')
    model,preparation = build_clamped(arrays,masses,authority='finite independent oracle',control=control)
    independent = control == 'independent-regions'
    h,linear,j,w,omega,field_mass = real_scalar(arrays,expected_mass,independent)
    close(prefix+'field-volume-mass',preparation['field_mass'],field_mass)
    close(prefix+'trap-frequency',model.omega,omega)
    close(prefix+'ionic-product-preparation',preparation['ionic_ordered'],w[:6,:6])
    close(prefix+'field-product-preparation',preparation['field_ordered'],w[12:18,12:18])
    left,right = model.evolve(.23),model.evolve(.11)
    st,bt = evolve(h,linear,j,.23)
    ss,_ = evolve(h,linear,j,.11)
    equal,unequal,comm = st@w@st.T,st@w@ss.T,st@(1j*j)@ss.T
    anchor,coordinate = (1,0,0),1
    local,common = model.kicks(anchor,coordinate)
    for name,kick in (('baseline',None),('local',local),('common',common)):
        initial = np.zeros(len(h))
        amplitude = np.sqrt(2*expected_mass*AMU_RY*omega[coordinate])
        if name == 'local':
            initial[6+3+coordinate] = amplitude
        elif name == 'common':
            initial[3+coordinate] = initial[6+3+coordinate] = amplitude/np.sqrt(2)
        actual = model.mean(left,kick)
        actual_real = np.concatenate((actual['ions'].ravel(),actual['field'].ravel()))
        expected = st@initial+bt
        close(prefix+name+'-affine-mean',actual_real,expected)
        diagnostic = model.diagnostics(left,actual['q'])
        direct_energy = np.sum(h*equal.T)/2+expected@h@expected/2+linear@expected
        close(prefix+name+'-total-energy',diagnostic['energy_ry'],direct_energy)
        zero = model.evolve(0.)
        e0 = model.diagnostics(zero,model.mean(zero,kick)['q'])['mean_energy_ry']
        close(prefix+name+'-same-added-energy',e0,0 if name == 'baseline' else omega[coordinate])
    values = connected_readouts(model,left,right,anchor,[coordinate,3+coordinate])
    target = [6+coordinate,6+3+coordinate]
    for name,matrix in (('spatial_equal_time_ordered',equal),('spatial_two_time_ordered',unequal),
                        ('spatial_two_time_commutator',comm)):
        expected = np.asarray([matrix[np.ix_([6*r+coordinate,6*r+3+coordinate],target)] for r in range(2)])
        close(prefix+name,values[name],expected)
    field_target = np.arange(18,24) if independent else np.arange(12,18)
    for name,matrix in (('field_equal_time_ordered',equal),('field_two_time_ordered',unequal)):
        expected = np.asarray([matrix[np.ix_(np.arange(12+6*f,18+6*f),field_target)]
                               for f in range(2 if independent else 1)])
        close(prefix+name,values[name],expected)
    expected_mixed = np.asarray([unequal[np.ix_([6*r+coordinate,6*r+3+coordinate],field_target)] for r in range(2)])
    close(prefix+'mixed-two-time',values['ion_field_two_time_ordered'],expected_mixed)
    # Mutations expose the specific joins, not general ClampedArray equations.
    wrong = linear.copy()
    for r in range(2):
        wrong[6*r:6*r+3] = -arrays['ionic_force_ry_bohr'].ravel()
    killed(prefix+'bare-ionic-force-substitution',evolve(h,wrong,j,.23)[1],bt)
    wrong = linear.copy(); wrong[12:15] *= -1
    killed(prefix+'field-gradient-sign',evolve(h,wrong,j,.23)[1],bt)
    wrong_h = h.copy(); wrong_h[:3,:3] += np.eye(3)*.4
    killed(prefix+'double-counted-onsite-curvature',evolve(wrong_h,linear,j,.23)[0]@w,equal@np.linalg.inv(st.T))
    killed(prefix+'unequal-common-energy',2*omega[coordinate],omega[coordinate])
    if independent:
        killed(prefix+'retained-global-field-mass',2*field_mass,field_mass)
        killed(prefix+'retained-global-field-vacuum',preparation['field_ordered'][:3,:3]/2,w[12:15,12:15])
        shared = real_scalar(arrays,expected_mass,False)
        gs,_ = evolve(shared[0],shared[1],shared[2],.23)
        shared_cross = (gs@shared[3]@gs.T)[0:6,6:12]
        killed(prefix+'onsite-or-global-shared-not-uncoupled',shared_cross,equal[0:6,6:12])
    if variant != 'original':
        original_mass = arrays['atomic_masses_amu'][0]*AMU_RY
        local_h = model.local.h[0] if independent else model.h[0]
        coupling = local_h[3:6,6:9]
        expected_coupling = -np.eye(3)*arrays['valence_charges'][0]/(expected_mass*AMU_RY)*(1 if independent else np.sqrt(2))
        close(prefix+'isotope-canonical-coupling',coupling,expected_coupling)
        killed(prefix+'stale-original-canonical-mass',-np.eye(3)*arrays['valence_charges'][0]/original_mass*(1 if independent else np.sqrt(2)),expected_coupling)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    sources = source_bindings()
    errors,mutations = {},{}
    def close(name,actual,expected):
        error = float(np.linalg.norm(np.asarray(actual)-expected)/max(1.,np.linalg.norm(expected)))
        errors[name] = error
        if not np.isfinite(error) or error > TOL:
            raise AssertionError((name,error))
    def killed(name,wrong,expected):
        error = float(np.linalg.norm(np.asarray(wrong)-expected)/max(1.,np.linalg.norm(expected)))
        mutations[name] = error
        if not np.isfinite(error) or error <= 100*TOL:
            raise AssertionError(('insensitive mutation',name,error))
    masses = nist_hydrogen_masses()
    close('saved-NIST-H1',masses['H1'],1.00782503223)
    close('saved-NIST-H2',masses['H2'],2.01410177812)
    for variant in ('original','H1','H2'):
        for control in ('periodic','independent-regions'):
            checks_for_case(fixture(),variant,control,close,killed)
    result = dict(schema='P254-matched-clamped-oracle-v1',status='passed',source_bindings=sources,
        tolerance=TOL,errors=errors,mutation_separations=mutations,
        scope='Finite independent real-space scalar oracle for new matched joins and local-field normalization only; no material evidence.')
    with args.output.open('x') as stream:
        json.dump(pack(result),stream,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(status='passed',checks=len(errors),mutations=len(mutations),output=str(args.output))))


if __name__ == '__main__':
    main()
