#!/usr/bin/env python3
"""Independent SymPy/Fock/DFT checks for the conditional ionic preparation.

Root runs this finite mathematical oracle; it does not read material arrays.
The finite Fock vacuum moments are exact, not a claim of finite-dimensional CCR
as an operator identity. A report is emitted only after every assertion passes.
"""
from pathlib import Path
from types import SimpleNamespace
import argparse
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
for directory in (ROOT/'src', HERE.parent/'0025', HERE.parent/'0024'):
    sys.path.insert(0, str(directory))

import numpy as np
import sympy as sp
from material_run_contract import binding
from substrate_framework.ionic_gaussian_preparation import prepare_local_trap_vacuum, MissingTrapPreparation
from substrate_framework.physical_bilinear_preparation import IonicPreparation
from run_ionic_gaussian_material import source_paths


def close(actual, expected):
    if not np.allclose(actual, expected, rtol=2e-12, atol=2e-12):
        raise AssertionError(('independent oscillator oracle mismatch', actual, expected))


def rejected(call, exception=ValueError):
    try:
        call()
    except exception:
        return True
    raise AssertionError('invalid preparation accepted')


def symbolic_block():
    # SI L/ Ry with tau=t Ry/hbar and x=a0*u; derive the Ry inertia,
    # independently of the existing QE helper's implementation.
    msi, a0, hb, me, ry, velocity = sp.symbols('m_si a0 hbar m_e Ry v', positive=True)
    kinetic = msi*(a0*ry/hb*velocity)**2/(2*ry)
    converted = sp.simplify(kinetic.subs(ry, hb**2/(2*me*a0**2)))
    assert sp.diff(converted, velocity, 2) == msi/(2*me)
    m, k, x = sp.symbols('M kappa x', positive=True)
    beta = sp.sqrt(m*k)
    psi = (beta/sp.pi)**sp.Rational(1,4)*sp.exp(-beta*x*x/2)
    hpsi = -sp.diff(psi,x,2)/(2*m)+k*x*x*psi/2
    assert sp.simplify(hpsi/psi-sp.sqrt(k/m)/2) == 0
    # x is unrestricted for integration; the positive mass/stiffness determine
    # the convergent Gaussian, not a numerical production formula.
    y = sp.Symbol('y', real=True)
    wave = psi.subs(x,y)
    uu = sp.simplify(sp.integrate(y*y*wave**2,(y,-sp.oo,sp.oo)))
    pp = sp.simplify(sp.integrate(-wave*sp.diff(wave,y,2),(y,-sp.oo,sp.oo)))
    up = sp.simplify(sp.integrate(-sp.I*wave*y*sp.diff(wave,y),(y,-sp.oo,sp.oo)))
    assert sp.simplify(uu-1/(2*sp.sqrt(m*k))) == 0
    assert sp.simplify(pp-sp.sqrt(m*k)/2) == 0
    assert up == sp.I/2 and sp.simplify(uu*pp) == sp.Rational(1,4)
    return dict(SI_to_Ry_inertia=str(sp.diff(converted,velocity,2)),
                uu=str(uu), pp=str(pp), up=str(up), vacuum_energy=str(sp.sqrt(k/m)/2))


def fixture():
    lattice = np.zeros((4,4,4,2,3,2,3), complex)
    stiffness = np.array([2.,3.,5.,7.,11.,13.])
    lattice[0,0,0] = np.diag(stiffness).reshape(2,3,2,3)
    # Large negative intercell terms make Gamma unstable. The preparation must
    # use the clamped R=0 diagonal and leave the dynamics completely unchanged.
    lattice[1,0,0] = (-20*np.eye(6)).reshape(2,3,2,3)
    lattice[3,0,0] = (-20*np.eye(6)).reshape(2,3,2,3)
    mass_amu = np.array([1.25, 7.5])
    return lattice, mass_amu, stiffness


def fock_block():
    lattice, atomic, stiffness = fixture()
    original = lattice.copy()
    state = prepare_local_trap_vacuum(lattice, atomic, state_authority='synthetic trap-vacuum oracle')
    assert np.array_equal(lattice,original)
    # Independent source constants as decimal SymPy rationals, not AMU_RY import.
    conversion = sp.Rational('1.66053906660e-27')/(2*sp.Rational('9.1093837015e-31'))
    mass = np.repeat(atomic*float(conversion),3)
    close(state.masses_ry,mass)
    close(state.omega_ry,np.sqrt(stiffness/mass))
    a = sp.zeros(4)
    for j in range(1,4):
        a[j-1,j] = sp.sqrt(j)
    vacuum = sp.Matrix([1,0,0,0])
    n = state.ndof
    expected = np.zeros((2*n,2*n),complex)
    for j in range(n):
        m = sp.Rational(str(atomic[j//3]))*conversion
        k = sp.Rational(str(stiffness[j]))
        width = sp.sqrt(m*k)
        x = (a+a.T)/sp.sqrt(2*width)
        p = -sp.I*sp.sqrt(width/2)*(a-a.T)
        for i, left in enumerate((x,p)):
            for h, right in enumerate((x,p)):
                expected[i*n+j,h*n+j] = complex((vacuum.T*left*right*vacuum)[0].evalf())
        assert sp.simplify((vacuum.T*(x*p-p*x)*vacuum)[0]) == sp.I
    close(state.local.ordered,expected)
    v = state.symmetric
    jmat = np.block([[np.zeros((n,n)),np.eye(n)],[-np.eye(n),np.zeros((n,n))]])
    close(state.local.ordered-state.local.ordered.T,1j*jmat)
    close(v@jmat@v,jmat/4)
    assert np.linalg.eigvalsh(state.local.ordered).min() > -1e-11
    # Explicit non-diagonal electronic density tensor ionic vacuum: means of
    # electrons may be nonzero, but centered mixed moments vanish in both orders.
    rhoe = sp.Matrix([[sp.Rational(2,3),sp.I/6],[-sp.I/6,sp.Rational(1,3)]])
    electronic = sp.Matrix([[1,sp.I],[-sp.I,2]])
    rho = sp.kronecker_product(rhoe,vacuum*vacuum.T)
    e = sp.kronecker_product(electronic,sp.eye(4))
    for ion in (a+a.T,-sp.I*(a-a.T)):
        z = sp.kronecker_product(sp.eye(2),ion)
        assert sp.trace(rho*e*z)-sp.trace(rho*e)*sp.trace(rho*z) == 0
        assert sp.trace(rho*z*e)-sp.trace(rho*e)*sp.trace(rho*z) == 0
    cross = state.product_cross([SimpleNamespace(h0=np.eye(2))])
    assert cross.matrices[0].shape == (2*n,2,2) and not np.any(cross.matrices[0])
    mutations = dict(
        reversed_up_sign=rejected(lambda:IonicPreparation(expected.T,'wrong sign')),
        erased_up=rejected(lambda:IonicPreparation(v,'classical covariance is not quantum vacuum')),
        wrong_Hartree_mass_exposure=float(np.max(np.abs(mass-2*mass))),
        wrong_mass_vacuum_exposure=float(np.max(np.abs(
            np.diag(v)[:n]-1/(2*np.sqrt(2*mass*stiffness))))),
        wrong_frequency_exposure=float(np.max(np.abs(state.omega_ry-np.sqrt(stiffness/(2*mass))))))
    assert all(value is True or value > 1e-6 for value in mutations.values())
    for name, value in [('negative',-1.),('zero',0.),('complex',2.+.01j)]:
        changed = lattice.copy()
        changed[0,0,0,0,0,0,0] = value
        before = changed.copy()
        mutations[name+'_stiffness'] = rejected(lambda:prepare_local_trap_vacuum(
            changed,atomic,state_authority='synthetic refusal'),MissingTrapPreparation)
        assert np.array_equal(changed,before)
    small = lattice.copy()
    small[0,0,0,0,0,0,0] += 1e-12j
    retained = prepare_local_trap_vacuum(small,atomic,state_authority='explicit tolerance')
    assert retained.stiffness_raw[0].imag == 1e-12
    close(retained.local.ordered,state.local.ordered)
    assert np.linalg.eigvalsh(np.fft.fftn(lattice,axes=(0,1,2))[0,0,0].reshape(6,6)).max() < 0
    return state, dict(Fock_vacuum_matches=True,uncertainty_saturated=True,
        product_cross_exact=True,dynamically_unstable_fixture_unchanged=True,mutations=mutations)


def fourier_block(state):
    shape = state.grid_shape
    indices = np.array(list(np.ndindex(shape)))
    count = len(indices)
    # Direct character sum independently checks numpy FFT normalization/sign.
    transform = np.exp(-2j*np.pi*(indices/np.asarray(shape))@indices.T)/np.sqrt(count)
    basis = np.eye(count).reshape(shape+(count,))
    fft = np.fft.fftn(basis,axes=(0,1,2),norm='ortho').reshape(count,count)
    close(fft,transform)
    ordered_delta = transform@transform.T
    dagger_delta = transform@transform.conj().T
    w = state.local.ordered
    for i,q in enumerate(indices):
        opposite = (-q)%np.asarray(shape)
        other_index = np.flatnonzero(np.all(indices==opposite,axis=1))[0]
        close(transform[i].conj(),transform[other_index])
        for h,r in enumerate(indices):
            close(state.q_pair(q,r),ordered_delta[i,h]*w)
            close(state.q_pair(q,r,dagger=True),dagger_delta[i,h]*w)
    close(state.real_pair((0,0,0),(1,0,0)),0*w)
    close(state.real_pair((0,0,0),(4,0,0)),w)
    mean_weights = np.ones(count)/count
    close(state.mean_ordered,(mean_weights@mean_weights)*w)
    close(state.mean_ordered-state.mean_ordered.T,(w-w.T)/64)
    assert not np.any(state.q_pair((1,0,0),(1,0,0)))
    close(state.q_pair((1,0,0),(3,0,0)),w)
    mutations = dict(q_not_minus_q_exposure=float(np.max(np.abs(
        state.q_pair((1,0,0),(1,0,0),dagger=True)-state.q_pair((1,0,0),(1,0,0))))),
        omitted_cell_factor_exposure=float(np.max(np.abs(w-state.mean_ordered))),
        duplicated_cell_factor_exposure=float(np.max(np.abs(state.mean_ordered-state.mean_ordered/64))),
        noncanonical_mean_rejected=rejected(lambda:IonicPreparation(state.mean_ordered,'wrong mean CCR')))
    assert all(value is True or value > 1e-5 for value in mutations.values())
    payload = state.archive_arrays()
    close(payload['ordered_up'],.5j*np.eye(state.ndof))
    close(payload['symmetrized_up'],np.zeros((state.ndof,state.ndof)))
    assert np.array_equal(payload['q_minus_index'],(-payload['q_index'])%np.asarray(shape))
    return dict(periodic_cells=count,independent_DFT_FFT=True,q_conjugacy=True,
                mean_covariance_factor='1/64',mutations=mutations)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    sources = [binding(path) for path in source_paths()]
    symbolic = symbolic_block()
    state, fock = fock_block()
    fourier = fourier_block(state)
    result = dict(schema='P254-ionic-gaussian-oracle-v1',status='passed',
        scope='finite independent mathematical oracle, no scientific material execution',
        symbolic=symbolic,fock=fock,fourier=fourier,source_bindings=sources)
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
