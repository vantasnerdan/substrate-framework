#!/usr/bin/env python3
"""Independent dense/SymPy full coupled oracle against the new disk runtime.

Synthetic complete four-sector/two-spin model, not material production. Dense
PW-square matrices are confined to this independent three-orbital oracle.
"""
import argparse
from contextlib import contextmanager
import json
from pathlib import Path
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT/'src'))
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from substrate_framework.prepared_reference import PreparedReference
from substrate_framework.material_joint_field_scalar import JointFieldCoefficients
from substrate_framework.physical_field_environment import UniformMaxwellField
from substrate_framework.streamed_material_state import DiskMean, file_binding
from substrate_framework.streamed_joint_state import DiskJointMean, JointBudget, JointBudgetExceeded
from substrate_framework.streamed_joint_reduction import JointStage
from substrate_framework.streamed_joint_mean import StreamedJointMean, JointConvergenceError


def fixture():
    rng = np.random.default_rng(743927)
    def hermitian(scale):
        a = rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
        return scale*(a+a.conj().T)
    frames = [(rng.normal(size=(3,2))+1j*rng.normal(size=(3,2)))*.3 for _ in range(4)]
    ds = []
    for i in range(4):
        x = rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
        if i == 0:
            x[:,1] = 0  # Singular original D is retained without inversion.
        ds.append(x @ x.conj().T*.13)
    hs = np.asarray([hermitian(.16) for _ in range(4)])
    gs = np.asarray([[hermitian(.12) for _ in range(3)] for _ in range(4)])
    ls = np.asarray([[hermitian(.17) for _ in range(3)] for _ in range(4)])
    measure = .7
    nbar = np.full((3,2), .08)
    kernel = np.array([[.7,-.13],[-.13,.6]])
    mass = np.array([2.,3.,4.])
    c = np.array([[.8,.1,-.04],[.1,.6,.03],[-.04,.03,.9]])
    dc = np.diag([.23,-.13,.19])
    contact = np.array([[.9,.1,-.04],[.1,.7,.08],[-.04,.08,1.1]])
    mixed = np.array([[.2,-.1,.03],[.05,.17,.04],[-.03,.08,.21]])
    force0 = np.array([.13,-.06,.08])
    field = UniformMaxwellField(6*measure,np.eye(3),'synthetic two-cell full Cartesian Maxwell')
    coefficients = JointFieldCoefficients(dc,contact,mixed,'synthetic fixed scalar coefficients')
    class Coupling:
        npw, nr, ndof = 3,3,3
        grid_measure_bohr3, cell_volume_bohr3 = measure, 3*measure
        def __init__(self, i):
            self.i = i
        def moments(self, left, right, d):
            rho = left @ d @ right.conj().T
            return np.diag(rho)/measure, np.einsum('iab,ba->i',gs[self.i],rho)
        def potential_action(self, potential, v):
            return potential[:,None]*v
        def ionic_action(self, u, v):
            return np.einsum('i,iab->ab',u,gs[self.i]) @ v
    class Vertex:
        def __init__(self, i):
            self.i = i
        def moments(self, psi, d):
            rho = psi @ d @ psi.conj().T
            return dict(field_gradient=np.einsum('iab,ba->i',ls[self.i],rho))
        def prepare(self, a):
            matrix = np.einsum('i,iab->ab',a,ls[self.i])
            class Action:
                def apply(self, v):
                    return matrix @ v
            return Action()
    active, completed, barriers = [], [], []
    @contextmanager
    def factory(slot):
        i = int(slot.identity)
        active.append(i)
        try:
            yield Coupling(i), lambda v: hs[i] @ v, Vertex(i)
        finally:
            active.remove(i)
            completed.append(i)
    def checked_kernel(n):
        assert not active and completed[-4:] == [0,1,2,3]
        barriers.append(len(completed))
        return n @ kernel
    reference = PreparedReference(nbar,force0,np.zeros_like(nbar),np.zeros(3),c,mass,
        np.zeros(3),checked_kernel,2,'synthetic scalar','synthetic ions','all-four-slots',
        'zero support','atom-major xyz','fixed original complex D')
    return dict(frames=frames, ds=ds, hs=hs, gs=gs, ls=ls, measure=measure, nbar=nbar,
        kernel=kernel, mass=mass, c=c, dc=dc, contact=contact, mixed=mixed, force0=force0,
        field=field, coefficients=coefficients, factory=factory, reference=reference,
        barriers=barriers, charge=1.3)


def independent_scalar(f):
    """Symbolic scalar differentiated independently, not coordinate_correction."""
    u = sp.Matrix(sp.symbols('u0:3', real=True))
    p = sp.Matrix(sp.symbols('p0:3', real=True))
    a = sp.Matrix(sp.symbols('a0:3', real=True))
    pi = sp.Matrix(sp.symbols('pi0:3', real=True))
    g = sp.Matrix(sp.symbols('g0:3', real=True))
    ell = sp.Matrix(sp.symbols('l0:3', real=True))
    coords = list(u)+list(p)+list(a)+list(pi)
    scalar = 2*((u.dot(g-sp.Matrix(f['force0'])))
        +(u.T*sp.Matrix(f['c']+f['dc'])*u)[0]/2+a.dot(ell)
        +(a.T*sp.Matrix(f['contact'])*a)[0]/2+(u.T*sp.Matrix(f['mixed'])*a)[0]
        +sum((p[i]-f['charge']*a[i])**2/(2*f['mass'][i]) for i in range(3)))
    scalar += pi.dot(pi)/(2*f['field'].mass)
    # Exact identity: retain each declared binary coefficient as a rational,
    # rather than simplify an inexact floating polynomial to literal zero.
    scalar = scalar.xreplace({x: sp.Rational(x) for x in scalar.atoms(sp.Float)})
    rhs = [sp.diff(scalar,p[i])/2 for i in range(3)]+[-sp.diff(scalar,u[i])/2 for i in range(3)]
    rhs += [sp.diff(scalar,pi[i]) for i in range(3)]+[-sp.diff(scalar,a[i]) for i in range(3)]
    energy = sp.lambdify(coords+list(g)+list(ell), scalar, 'numpy')
    derivative = sp.lambdify(coords+list(g)+list(ell), rhs, 'numpy')
    # Symbolic Hamiltonian conservation of the complete canonical coordinate block.
    assert sp.simplify(sum(sp.diff(scalar,x)*v for x,v in zip(coords,rhs))) == 0
    return energy, derivative


def pack(frames, coords):
    return np.concatenate([np.asarray(frames).ravel(), np.asarray(coords).ravel()]).astype(complex)


def unpack(y):
    return y[:24].reshape(4,3,2), y[24:].real.reshape(4,3)


def dense(f, functions, y, mutation=None):
    frames, (u,p,a,pi) = unpack(y)
    rhos = [psi @ d @ psi.conj().T for psi,d in zip(frames,f['ds'])]
    weight = 1. if mutation == 'weights' else .5
    density = np.zeros((3,2), complex)
    g, ell, bare = np.zeros(3,complex), np.zeros(3,complex), 0j
    for i,rho in enumerate(rhos):
        density[:,i%2] += weight*np.diag(rho)/f['measure']
        g += weight*np.einsum('iab,ba->i',f['gs'][i],rho)
        ell += weight*np.einsum('iab,ba->i',f['ls'][i],rho)
        bare += weight*np.trace(f['hs'][i] @ rho)
    delta = density.real-f['nbar']
    potential = delta @ f['kernel']
    coordinate_rhs = np.array(functions[1](*u,*p,*a,*pi,*g.real,*ell.real))
    if mutation == 'field-backreaction':
        coordinate_rhs[9:] = 0
    if mutation == 'minimal-coupling':
        coordinate_rhs[:3] = p/f['mass']
        coordinate_rhs[9:] -= 2*f['charge']*(p-f['charge']*a)/f['mass']
    if mutation == 'DeltaC':
        coordinate_rhs[3:6] += f['dc'] @ u
    hamiltonians = [f['hs'][i]+np.diag(potential[:,i%2])
        +np.einsum('i,iab->ab',u,f['gs'][i])+np.einsum('i,iab->ab',a,f['ls'][i]) for i in range(4)]
    rhs = pack([-1j*h @ psi for h,psi in zip(hamiltonians,frames)], coordinate_rhs)
    energy = float((2*bare+f['measure']*np.sum(delta*potential)).real
        +functions[0](*u,*p,*a,*pi,*g.real,*ell.real))
    return rhs, energy, hamiltonians


def budget(**kwargs):
    return JointBudget(**dict(dict(max_actions=1000000,max_seconds=600.,max_bytes=100000000),**kwargs))


def state_vector(state):
    return pack([state.arrays(s)[0] for s in state.slots],
                [state.u,state.p,state.coordinate,state.momentum])


def run(directory):
    f = fixture()
    functions = independent_scalar(f)
    coords = np.array([[.07,-.09,.11],[.12,-.08,.15],[.04,-.06,.03],[-.02,.01,.025]])
    records = [(dict(identity=str(i),component=i%2,weight=.5,basis_identity='basis-'+str(i),
        gauge_identity='raw-'+str(i),source_identity='source-'+str(i)), psi,d)
        for i,(psi,d) in enumerate(zip(f['frames'],f['ds']))]
    original = DiskMean.create(directory/'original',records,identity='OLD-scalar',time=0.,u=coords[0],p=coords[1])
    original_binding = file_binding(original.path)
    state = DiskJointMean.wrap_original(original_binding,directory/'joint',identity='NEW-joint-scalar',
        coordinate=coords[2],momentum=coords[3],expected_slots=tuple(str(i) for i in range(4)),
        evidence=dict(new_scalar='independent coupled synthetic scalar',field_preparation='explicit means'),budget=budget())
    kwargs = dict(model_identity=state.identity,expected_slots=tuple(s.identity for s in state.slots),grid_measure=f['measure'])
    def model(**options):
        return StreamedJointMean(f['reference'],f['factory'],f['coefficients'],f['field'],[f['charge']],
                                 **kwargs,budget=options.pop('allocation',budget()),**options)
    m = model()
    stage = m.reduce(state)
    y0 = state_vector(state)
    rhs, energy, hs = dense(f,functions,y0)
    observed = pack([-1j*h @ psi for h,psi in zip(hs,f['frames'])],
                    [stage.velocity,stage.force,stage.field_velocity,stage.field_force])
    assert np.max(np.abs(observed-rhs)) < 2e-12
    assert abs(stage.energy-energy) < 2e-12
    for i,slot in enumerate(state.slots):
        with m.operators(slot) as (c,h,v):
            assert np.max(np.abs(m.frozen_action(c,h,v,slot,stage)(np.eye(3))-hs[i])) < 2e-12
    interval = .12
    exact = solve_ivp(lambda t,y:dense(f,functions,y)[0],(0,interval),y0,method='DOP853',rtol=2e-12,atol=2e-14)
    assert exact.success
    errors, drifts, finals = [], [], []
    for steps in (4,8,16):
        end, history = m.evolve(state,directory/('steps-'+str(steps)),physical_interval=interval,
            steps=steps,max_steps=16,max_iterations=50,tolerance=2e-13,refinement_label=str(steps))
        restored = DiskJointMean(end.path,expected_identity=state.identity)
        assert np.array_equal(state_vector(restored),state_vector(end))
        assert all(a.coefficient == b.coefficient for a,b in zip(state.slots,end.slots))
        errors.append(float(np.linalg.norm(state_vector(end)-exact.y[:,-1])))
        drifts.append(abs(m.reduce(end).energy-energy))
        finals.append(end)
    assert 3.5 < errors[0]/errors[1] < 4.5 and 3.5 < errors[1]/errors[2] < 4.5
    assert errors[-1] < 2e-5
    final_stage = JointStage.restore(finals[-1].history[-1])
    assert final_stage.coordinate.shape == (3,)
    # Replaying the actual saved frozen action must reproduce the last full frame.
    prior = DiskJointMean(finals[-1].history[-1]['prior_boundary'])
    dt = finals[-1].history[-1]['dt_ry_inverse']
    replay_error = 0.
    for source,target in zip(prior.slots,finals[-1].slots):
        with m.operators(source) as (c,h,v):
            replay = m.propagate(c,h,v,source,final_stage,prior.arrays(source)[0],dt)
        replay_error = max(replay_error,float(np.linalg.norm(replay-finals[-1].arrays(target)[0])))
    assert replay_error < 2e-12
    # Injected dense exponential receives precisely the same frozen action.
    captured = []
    def dense_exponential(action, psi, dt):
        matrix = action(np.eye(len(psi),dtype=complex))
        assert np.linalg.norm(matrix-matrix.conj().T) < 2e-12
        captured.append(matrix)
        return expm(-1j*dt*matrix) @ psi
    injected = model(exponential=dense_exponential).step(state,directory/'injected',dt=.03,max_iterations=50,tolerance=2e-13)
    default = m.step(state,directory/'default',dt=.03,max_iterations=50,tolerance=2e-13)
    assert captured and np.max(np.abs(state_vector(injected)-state_vector(default))) < 2e-12
    mutations = {}
    for name in ('field-backreaction','minimal-coupling','DeltaC','weights'):
        altered = solve_ivp(lambda t,y:dense(f,functions,y,name)[0],(0,interval),y0,
                            method='DOP853',rtol=2e-12,atol=2e-14)
        assert altered.success
        mutations[name] = float(np.linalg.norm(altered.y[:,-1]-exact.y[:,-1]))
        assert mutations[name] > 20*errors[-1]
    stale = pack([expm(-1j*interval*h) @ psi for h,psi in zip(hs,f['frames'])],
                 coords+interval*rhs[24:].real.reshape(4,3))
    mutations['stale-stage'] = float(np.linalg.norm(stale-exact.y[:,-1]))
    assert mutations['stale-stage'] > 20*errors[-1]
    number_error = abs(m.reduce(finals[-1]).number-stage.number)
    gram_error = max(float(np.linalg.norm(finals[-1].arrays(s)[0].conj().T @ finals[-1].arrays(s)[0]
                          -psi.conj().T @ psi)) for s,psi in zip(finals[-1].slots,f['frames']))
    assert number_error < 2e-11 and gram_error < 2e-11 and f['barriers']
    refusals = []
    for name, allocation in [('action',budget(max_actions=1)),('storage',budget(max_bytes=1)),('time',budget())]:
        if name == 'time':
            allocation.started -= 1000
        try:
            model(allocation=allocation).step(state,directory/('cap-'+name),dt=.03,max_iterations=40,tolerance=2e-13)
        except JointBudgetExceeded:
            refusals.append(name)
        else:
            raise AssertionError('cap not enforced: '+name)
        assert not (directory/('cap-'+name)/'completed'/'boundary.json').exists()
    try:
        m.step(state,directory/'nonconverged',dt=.12,max_iterations=1,tolerance=1e-15)
    except JointConvergenceError:
        refusals.append('nonconvergence')
    else:
        raise AssertionError('incomplete midpoint was adopted')
    assert not (directory/'nonconverged'/'completed'/'boundary.json').exists()
    assert file_binding(original.path) == original_binding
    assert all(a.frame == b.frame and a.coefficient == b.coefficient for a,b in zip(original.slots,state.slots))
    return dict(energy_initial=energy,steps=[4,8,16],errors=errors,energy_drifts=drifts,
        convergence_ratios=[errors[0]/errors[1],errors[1]/errors[2]],mutations=mutations,
        saved_stage_replay_error=replay_error,number_error=number_error,gram_error=gram_error,
        cap_and_incomplete_stage_refusals=refusals,complete_barriers=len(f['barriers']),
        original_full_frames_D_unchanged=True,scalar='independent SymPy canonical derivatives and dense original-D traces',
        scope='synthetic full coupled homogeneous mean oracle; not executed material/local-q dynamics')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    args = p.parse_args()
    with tempfile.TemporaryDirectory(prefix='p254-joint-streamed-') as tmp:
        result = run(Path(tmp))
    paths = [HERE/'run_joint_streamed_material.py',Path(__file__),HERE/'joint-streamed-registration.json']
    paths += [ROOT/'src/substrate_framework'/name for name in (
        'streamed_joint_state.py','streamed_joint_reduction.py','streamed_joint_mean.py',
        'material_joint_field_scalar.py','material_physical_operators.py',
        'streamed_material_state.py','streamed_material_evolution.py')]
    result.update(schema='P254-joint-streamed-oracle-v1',status='passed',source_bindings=[file_binding(p) for p in paths])
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write('\n')


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()
