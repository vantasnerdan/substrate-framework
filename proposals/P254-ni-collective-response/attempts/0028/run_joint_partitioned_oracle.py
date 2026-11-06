#!/usr/bin/env python3
"""Source-bound independent dense/SymPy oracle of actual partitioned disk output.

Only the pre-existing synthetic fixture is reused, not its oracle results or
production scalar. This program must be launched by the allocating parent.
"""
import argparse
from contextlib import contextmanager
from copy import copy
from dataclasses import replace
import json
import os
from pathlib import Path
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE)]
for variable in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[variable] = '1'
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from run_joint_streamed_oracle import fixture
from run_joint_partitioned_material import implementation_bindings, source_identity
from substrate_framework.partitioned_joint_budget import PartitionBudget
from substrate_framework.partitioned_joint_mean import PartitionedJointMean
from substrate_framework.streamed_joint_mean import StreamedJointMean
from substrate_framework.streamed_joint_state import DiskJointMean, JointBudgetExceeded
from substrate_framework.streamed_joint_reduction import JointStage
from substrate_framework.streamed_material_state import DiskMean, file_binding


def pack(frames, coords):
    return np.concatenate([np.asarray(frames).ravel(), np.asarray(coords).ravel()]).astype(complex)


def vector(state):
    return pack([state.arrays(s)[0] for s in state.slots],
                [state.u,state.p,state.coordinate,state.momentum])


def independent(f):
    """Differentiate a separately assembled scalar, not coordinate_correction."""
    u,p,a,pi,g,ell = [sp.Matrix(sp.symbols(name+'0:3',real=True)) for name in ('u','p','a','pi','g','l')]
    scalar = 2*(u.dot(g-sp.Matrix(f['force0']))+(u.T*sp.Matrix(f['c']+f['dc'])*u)[0]/2
        +a.dot(ell)+(a.T*sp.Matrix(f['contact'])*a)[0]/2+(u.T*sp.Matrix(f['mixed'])*a)[0]
        +sum((p[i]-f['charge']*a[i])**2/(2*f['mass'][i]) for i in range(3)))
    scalar += pi.dot(pi)/(2*f['field'].mass)
    scalar = scalar.xreplace({x:sp.Rational(x) for x in scalar.atoms(sp.Float)})
    coordinates = list(u)+list(p)+list(a)+list(pi)
    rhs = [sp.diff(scalar,p[i])/2 for i in range(3)]+[-sp.diff(scalar,u[i])/2 for i in range(3)]
    rhs += [sp.diff(scalar,pi[i]) for i in range(3)]+[-sp.diff(scalar,a[i]) for i in range(3)]
    assert sp.simplify(sum(sp.diff(scalar,x)*v for x,v in zip(coordinates,rhs))) == 0
    scalar_fn = sp.lambdify(coordinates+list(g)+list(ell),scalar,'numpy')
    rhs_fn = sp.lambdify(coordinates+list(g)+list(ell),rhs,'numpy')

    def evaluate(y, mutation=None):
        frames = y[:24].reshape(4,3,2)
        u,p,a,pi = y[24:].real.reshape(4,3)
        density, gradient, optical, bare = np.zeros((3,2),complex), np.zeros(3,complex), np.zeros(3,complex), 0j
        for i,(psi,d) in enumerate(zip(frames,f['ds'])):
            rho = psi @ d @ psi.conj().T
            w = 1. if mutation == 'weights' else .5
            if mutation == 'missing-slot' and i == 2:
                continue
            density[:,i%2] += w*np.diag(rho)/f['measure']
            gradient += w*np.einsum('iab,ba->i',f['gs'][i],rho)
            optical += w*np.einsum('iab,ba->i',f['ls'][i],rho)
            bare += w*np.trace(f['hs'][i] @ rho)
        potential = (density.real-f['nbar']) @ f['kernel']
        coord_rhs = np.asarray(rhs_fn(*u,*p,*a,*pi,*gradient.real,*optical.real))
        matrices = [f['hs'][i]+np.diag(potential[:,i%2])+np.einsum('i,iab->ab',u,f['gs'][i])
                    +np.einsum('i,iab->ab',a,f['ls'][i]) for i in range(4)]
        if mutation == 'field':
            coord_rhs[9:] = 0
            matrices = [h-np.einsum('i,iab->ab',a,f['ls'][i]) for i,h in enumerate(matrices)]
        if mutation == 'H0':
            matrices = [h-f['hs'][i] for i,h in enumerate(matrices)]
        energy = float((2*bare+f['measure']*np.sum((density.real-f['nbar'])*potential)).real
                       +scalar_fn(*u,*p,*a,*pi,*gradient.real,*optical.real))
        return dict(rhs=pack([-1j*h @ psi for h,psi in zip(matrices,frames)],coord_rhs),
                    energy=energy, density=density, gradient=gradient, ell=optical, bare=bare,
                    number=f['measure']*density.sum(), hamiltonians=matrices)
    return evaluate


def allocation(**overrides):
    values = dict(max_actions=2000000,max_seconds=1800.,max_bytes=1000000000,
                  max_memory_bytes=12000000000,worker_memory_bytes=4000000000)
    values.update(overrides)
    return PartitionBudget(**values)


def rejected(action, exceptions=(ValueError,)):
    try:
        action()
    except exceptions:
        return True
    raise AssertionError('mutation/incomplete request was accepted')


def run(directory, cpus):
    registration = json.loads((HERE/'joint-partitioned-registration.json').read_text())
    assert registration['status'] == 'prospective-source-registration-not-numerical-evidence'
    bindings = implementation_bindings()
    source = source_identity(bindings)
    f = fixture()
    kernel_calls = []
    def kernel(n):
        kernel_calls.append(1)
        return n @ f['kernel']
    reference = replace(f['reference'], kernel=kernel)
    evaluate = independent(f)
    coords = np.array([[.07,-.09,.11],[.12,-.08,.15],[.04,-.06,.03],[-.02,.01,.025]])
    records = [(dict(identity=str(i),component=i%2,weight=.5,basis_identity='basis-'+str(i),
        gauge_identity='raw-'+str(i),source_identity='source-'+str(i)),psi,d)
        for i,(psi,d) in enumerate(zip(f['frames'],f['ds']))]
    old = DiskMean.create(directory/'original',records,identity='old',time=0.,u=coords[0],p=coords[1])
    state = DiskJointMean.wrap_original(file_binding(old.path),directory/'joint',identity='joint',
        coordinate=coords[2],momentum=coords[3],expected_slots=tuple(str(i) for i in range(4)),
        evidence=dict(new_scalar='full coupled synthetic',field_preparation='explicit'),budget=allocation())
    def serial(budget=None, factory=None):
        return StreamedJointMean(reference,factory or f['factory'],f['coefficients'],f['field'],[f['charge']],
            model_identity=state.identity,expected_slots=tuple(str(i) for i in range(4)),
            grid_measure=f['measure'],budget=budget or allocation())
    def partition(name, budget=None, factory=None, selected=None):
        return PartitionedJointMean.from_serial(serial(budget,factory),directory=directory/name,
            cpus=cpus if selected is None else selected,source_identity=source)
    m = partition('work')
    stage = m.reduce(state)
    direct = serial().reduce(state)
    initial = evaluate(vector(state))
    for name in ('density','gradient','field_gradient_electronic','potential','force','field_force','velocity'):
        assert np.max(np.abs(getattr(stage,name)-getattr(direct,name))) < 2e-12
    assert abs(stage.energy-direct.energy) < 2e-12 and abs(stage.energy-initial['energy']) < 2e-12
    assert abs(stage.number-initial['number']) < 2e-12
    assert np.max(np.abs(stage.density-initial['density'])) < 2e-12
    assert np.max(np.abs(stage.gradient-initial['gradient'])) < 2e-12
    assert np.max(np.abs(stage.field_gradient_electronic-initial['ell'])) < 2e-12
    observed = pack([-1j*h @ psi for h,psi in zip(initial['hamiltonians'],f['frames'])],
                    [stage.velocity,stage.force,stage.field_velocity,stage.field_force])
    assert np.max(np.abs(observed-initial['rhs'])) < 2e-12
    key, work, _ = m.executor.request(state)
    saved = m.executor.collect(work,key,state)
    bare = sum(s.weight*np.load(r['arrays']['one_body']['path']).item() for s,r in zip(state.slots,saved))
    assert abs(bare-initial['bare']) < 2e-12 and abs(bare) > 1e-4
    # Reusing a completed reduction performs exactly zero new operator actions.
    actions = m.budget.actions
    assert m.reduce(state).energy == stage.energy and m.budget.actions == actions
    worker_cpus = sorted({r['cpu'] for r in saved})
    assert set(worker_cpus).issubset(cpus) and all(r['blas_threads'] == 1 for r in saved)

    # Force the real disk/process path to stop after two completed slots.
    marker = directory/'fail-slot-2'
    marker.write_text('injected independent oracle failure')
    @contextmanager
    def failing_factory(slot):
        if slot.identity == '2' and marker.exists():
            raise RuntimeError('oracle intentional slot interruption')
        with f['factory'](slot) as operators:
            yield operators
    partial = partition('partial-work',factory=failing_factory,selected=(cpus[0],))
    before = len(kernel_calls)
    rejected(lambda:partial.reduce(state), (RuntimeError,))
    assert len(kernel_calls) == before
    pk, pd, _ = partial.executor.request(state)
    completed = [pd/('slot-%04d' % i)/'complete.json' for i in (0,1)]
    preserved = [(file_binding(p),p.stat().st_mtime_ns) for p in completed]
    rejected(lambda:partial.executor.collect(pd,pk,state))
    assert not list(pd.glob('**/boundary.json'))
    marker.unlink()
    resumed = partial.reduce(state)
    assert abs(resumed.energy-stage.energy) < 2e-12
    assert [(file_binding(p),p.stat().st_mtime_ns) for p in completed] == preserved
    assert partial.budget.actions == 4  # two before interruption, two after.

    # Missing/stale completed receipts are not silently accepted by a barrier.
    receipt = work/'slot-0002'/'complete.json'
    held = receipt.with_name('held.json')
    receipt.rename(held)
    rejected(lambda:m.executor.collect(work,key,state))
    held.rename(receipt)
    changed = replace(stage, coordinate=stage.coordinate+np.array([.03,0,0]),
                      shift=stage.shift+np.array([.03,0,0]))
    fk, fd, _ = m.executor.request(state,frozen_stage=stage,interval=.01)
    ck, cd, _ = m.executor.request(state,frozen_stage=changed,interval=.01)
    assert fk != ck
    target = cd/'slot-0000'
    target.mkdir()
    (target/'complete.json').write_bytes((work/'slot-0000'/'complete.json').read_bytes())
    rejected(lambda:m.executor.collect(cd,ck,state))
    (target/'complete.json').rename(target/'rejected-stale.json')
    altered_stage = m.reduce(state,frozen_stage=changed,interval=.01)
    normal_stage = m.reduce(state,frozen_stage=stage,interval=.01)
    field_request_difference = float(np.linalg.norm(altered_stage.density-normal_stage.density))
    assert field_request_difference > 1e-7
    weighted = copy(state)
    weighted.slots = (replace(state.slots[0],weight=.6),)+state.slots[1:]
    rejected(lambda:m.reduce(weighted))
    missing = copy(state)
    missing.slots = state.slots[:-1]
    rejected(lambda:m.reduce(missing))

    interval = .12
    exact = solve_ivp(lambda t,y:evaluate(y)['rhs'],(0,interval),vector(state),method='DOP853',rtol=2e-12,atol=2e-14)
    assert exact.success
    errors, finals, energies = [], [], []
    for steps in (2,4,8):
        final, history = m.evolve(state,directory/('steps-'+str(steps)),physical_interval=interval,
            steps=steps,max_steps=8,max_iterations=50,tolerance=2e-13,refinement_label=str(steps))
        restored = DiskJointMean(final.path,expected_identity=state.identity)
        assert np.array_equal(vector(final),vector(restored))
        errors.append(float(np.linalg.norm(vector(final)-exact.y[:,-1])))
        energies.append(abs(m.reduce(final).energy-stage.energy))
        finals.append(final)
    ratios = [errors[i]/errors[i+1] for i in range(2)]
    assert all(3.3 < ratio < 4.7 for ratio in ratios) and errors[-1] < 1e-4
    # A complete evolve restart returns the same immutable boundary and no actions.
    before = m.budget.actions
    again, _ = m.evolve(state,directory/'steps-8',physical_interval=interval,steps=8,max_steps=8,
        max_iterations=50,tolerance=2e-13,refinement_label='8')
    assert again.path == finals[-1].path and m.budget.actions == before
    last = finals[-1]
    stage_record = last.history[-1]
    frozen = JointStage.restore(stage_record)
    prior = DiskJointMean(stage_record['prior_boundary'])
    replay = 0.
    for i,(s,t) in enumerate(zip(prior.slots,last.slots)):
        h = (f['hs'][i]+np.diag(frozen.potential[:,i%2])+np.einsum('i,iab->ab',frozen.u,f['gs'][i])
             +np.einsum('i,iab->ab',frozen.shift,f['ls'][i]))
        predicted = expm(-1j*stage_record['dt_ry_inverse']*h) @ prior.arrays(s)[0]
        replay = max(replay,float(np.linalg.norm(predicted-last.arrays(t)[0])))
    assert replay < 2e-12
    gram = max(float(np.linalg.norm(last.arrays(s)[0].conj().T @ last.arrays(s)[0]-psi.conj().T @ psi))
               for s,psi in zip(last.slots,f['frames']))
    number = abs(m.reduce(last).number-stage.number)
    assert gram < 2e-11 and number < 2e-11
    mutations = {}
    for name in ('missing-slot','weights','field','H0'):
        altered = solve_ivp(lambda t,y:evaluate(y,name)['rhs'],(0,interval),vector(state),method='DOP853',rtol=2e-12,atol=2e-14)
        assert altered.success
        mutations[name] = float(np.linalg.norm(altered.y[:,-1]-exact.y[:,-1]))
        assert mutations[name] > 20*errors[-1]
    stale = pack([expm(-1j*interval*h) @ psi for h,psi in zip(initial['hamiltonians'],f['frames'])],
                 coords+interval*initial['rhs'][24:].real.reshape(4,3))
    mutations['stale-stage'] = float(np.linalg.norm(stale-exact.y[:,-1]))
    assert mutations['stale-stage'] > 20*errors[-1]
    refusals = []
    for name, cap in [('actions',allocation(max_actions=1)),('storage',allocation(max_bytes=1)),
                      ('time',allocation()),('memory',allocation(max_memory_bytes=1,worker_memory_bytes=1))]:
        if name == 'time':
            cap.started -= 2000
        limited = partition('cap-'+name,budget=cap,selected=(cpus[0],))
        rejected(lambda:limited.reduce(state),(RuntimeError,JointBudgetExceeded))
        refusals.append(name)
    assert all(a.frame == b.frame and a.coefficient == b.coefficient for a,b in zip(old.slots,state.slots))
    assert all(a.coefficient == b.coefficient for a,b in zip(state.slots,last.slots))
    return dict(schema='P254-joint-partitioned-oracle-v1',status='passed',source_bindings=bindings,
        energy_initial=stage.energy, physical_H0_trace=bare.real, step_counts=[2,4,8], errors=errors,
        convergence_ratios=ratios, energy_drifts=energies, mutations=mutations,
        field_shift_saved_output_difference=field_request_difference, saved_history_dense_replay_error=replay,
        number_error=number, original_gram_error=gram, partial_receipts_reused_without_recomputation=True,
        restart_performed_zero_new_actions=True, strict_missing_slot_and_stale_receipt_refusals=True,
        cap_refusals=refusals, assigned_cpus=list(cpus), observed_worker_cpus=worker_cpus,
        completed_kernel_barriers=len(kernel_calls), budget=m.budget.snapshot(),
        scope='actual synthetic partitioned full coupled disk path; not material numerical evidence')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cpus',type=int,nargs='+',choices=range(6,11),required=True)
    p.add_argument('--work-directory',type=Path)
    args = p.parse_args()
    os.sched_setaffinity(0,{args.cpus[0]})
    if args.work_directory:
        args.work_directory.mkdir(parents=True,exist_ok=False)
        result = run(args.work_directory,tuple(args.cpus))
    else:
        with tempfile.TemporaryDirectory(prefix='p254-joint-partitioned-') as tmp:
            result = run(Path(tmp),tuple(args.cpus))
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write('\n')


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()
