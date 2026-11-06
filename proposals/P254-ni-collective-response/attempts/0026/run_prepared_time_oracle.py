#!/usr/bin/env python3
"""New finite-time oracle; symbolic checks FIRST, no old oracle or material run."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'src'))
CHECKS = []
RTOL, ATOL = 1e-8, 1e-10


def bind(path):
    path = Path(path).resolve()
    raw = path.read_bytes()
    return dict(path=str(path), bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())


def source_access():
    local = [Path(__file__), HERE/'prepared_time_fixture.py', HERE/'continuation-registration.json',
             HERE/'prepared-model-implementation-handoff.json', HERE/'reference-preparation-grounding.json',
             HERE/'action-dynamics-grounding.json', HERE/'new-implementation-source-binding-001.json']
    local += [ROOT/'src/substrate_framework'/name for name in (
        'prepared_time_evolution.py','prepared_time_observables.py','prepared_factor_dynamics.py',
        'prepared_covariance_source.py','prepared_reference.py','physical_factor_preparation.py',
        'physical_bilinear_preparation.py')]
    # Distribution metadata and original on-disk solver bytes BEFORE importing
    # or invoking SciPy. DOP853 coefficient table is independently bound as well.
    distribution = importlib.metadata.distribution('scipy')
    paths = ['scipy/__init__.py','scipy/version.py','scipy/integrate/_ivp/ivp.py',
             'scipy/integrate/_ivp/base.py','scipy/integrate/_ivp/common.py',
             'scipy/integrate/_ivp/rk.py','scipy/integrate/_ivp/dop853_coefficients.py']
    return dict(local=[bind(p) for p in local], scipy=dict(version=distribution.version,
                sources=[bind(distribution.locate_file(p)) for p in paths]),
                access_order='source bytes bound before SciPy import/reference solve')


def exact(name, value):
    import sympy as s
    residual = s.expand(value)
    passed = residual == 0
    CHECKS.append(dict(name=name, passed=passed, symbolic_residual=str(residual)))
    if not passed:
        raise AssertionError(name)


def symbolic_first():
    import sympy as s
    a = s.Matrix([[0,0,0,0],[s.Rational(1,2),0,0,0],[0,s.Rational(1,2),0,0],[0,0,1,0]])
    b = s.Matrix([[s.Rational(1,6),s.Rational(1,3),s.Rational(1,3),s.Rational(1,6)]])
    one = s.ones(4,1)
    c = a*one
    conditions = [('order1',(b*one)[0]-1), ('order2',(b*c)[0]-s.Rational(1,2)),
                  ('order3_bush',(b*c.applyfunc(lambda x:x*x))[0]-s.Rational(1,3)),
                  ('order3_chain',(b*a*c)[0]-s.Rational(1,6)),
                  ('order4_bush',(b*c.applyfunc(lambda x:x**3))[0]-s.Rational(1,4)),
                  ('order4_branch',(b*s.diag(*c)*a*c)[0]-s.Rational(1,8)),
                  ('order4_fork',(b*a*c.applyfunc(lambda x:x*x))[0]-s.Rational(1,12)),
                  ('order4_chain',(b*a*a*c)[0]-s.Rational(1,24))]
    for name, residual in conditions:
        exact('RK4_'+name, residual)
    h,v,p,q,d,l,r = s.symbols('H V Psi PsiAdj D L R', commutative=False)
    rho, delta = p*d*q, l*q+p*r
    actual = -s.I*(h*l+v*p*d)*q+l*(s.I*q*h)+(-s.I*h*p)*r+p*(s.I*(r*h+d*q*v))
    exact('moving_connection_original_commutator', actual+s.I*(h*delta-delta*h+v*rho-rho*v))
    z = s.symbols('z')
    exact('RK4_stability_polynomial', (1+z*(b*(s.eye(4)+z*a+z*z*a*a+z**3*a**3)*one)[0])
          -(1+z+z*z/2+z**3/6+z**4/24))


def compare(name, actual, expected, *, rtol=RTOL, atol=ATOL):
    import numpy as np
    actual, expected = np.asarray(actual), np.asarray(expected)
    residual = float(np.max(np.abs(actual-expected), initial=0))
    passed = bool(np.all(np.isfinite(actual)) and np.allclose(actual,expected,rtol=rtol,atol=atol))
    CHECKS.append(dict(name=name, passed=passed, max_absolute_error=residual, rtol=rtol, atol=atol))
    if not passed:
        raise AssertionError(name+': '+str(residual))
    return residual


def sensitive(name, wrong, expected):
    import numpy as np
    residual = float(np.max(np.abs(np.asarray(wrong)-expected), initial=0))
    passed = not np.allclose(wrong, expected, rtol=RTOL, atol=ATOL)
    CHECKS.append(dict(name=name,passed=bool(passed),max_absolute_difference=residual))
    if not passed:
        raise AssertionError('insensitive mutation: '+name)


def packed_boundary(f, boundary):
    import numpy as np
    from prepared_time_fixture import factor_vector
    return np.array([factor_vector(f.model,boundary.mean)] +
                    [factor_vector(f.model,boundary.mean,getattr(s,p)) for s in boundary.sources
                     for p in ('electronic','ionic','cross')])


def reference(f, initial, source_binding):
    import numpy as np
    import scipy
    from scipy.integrate import solve_ivp
    from prepared_time_fixture import original_rhs
    from scipy.integrate._ivp.rk import DOP853
    bound = {p['path']: p for p in source_binding['sources']}
    if scipy.__version__ != source_binding['version']:
        raise RuntimeError('loaded SciPy differs from pre-bound version')
    for path in (scipy.__file__, solve_ivp.__code__.co_filename,
                 DOP853.step.__code__.co_filename, DOP853._step_impl.__code__.co_filename):
        actual = bind(path)
        if bound.get(actual['path']) != actual:
            raise RuntimeError('loaded solver source differs from pre-bound bytes')
    # The complete original 3x3 matrices, not factor RHS, use a different
    # high-order adaptive method. Only this tiny oracle retains a trajectory.
    times = np.linspace(0., .1, 33)
    solution = solve_ivp(lambda t,y:original_rhs(f,y), (0.,.1), initial.ravel(),
                         method='DOP853', t_eval=times, rtol=1e-12, atol=1e-13)
    if not solution.success or solution.y.shape[1] != len(times):
        raise RuntimeError('original-basis reference failed: '+solution.message)
    return solution.y.T.reshape(33,*initial.shape), dict(method='DOP853', version=scipy.__version__,
                nfev=solution.nfev, rtol=1e-12, atol=1e-13, times=times.tolist())


def contracts(f, bindings):
    from substrate_framework.prepared_time_evolution import TimeLimits, RealityTolerance
    identity = dict(model='finite original 3-orbital mixed-spin scalar v1',
                    preparation='independent Bernoulli + Gaussian product and conditional displaced joint mixture',
                    coordinates=f.reference.coordinate_identity,
                    sources=hashlib.sha256(json.dumps(bindings,sort_keys=True).encode()).hexdigest(),
                    budget='new finite-time benchmark 16/32, each branch cumulative32-step cap')
    return identity, TimeLimits(32,128,1536,90.), RealityTolerance(atol=1e-12,rtol=1e-10)


def probes(f):
    import numpy as np
    from substrate_framework.prepared_covariance_source import total_force_probe
    from substrate_framework.prepared_time_observables import regional_number_probe
    return dict(force=total_force_probe(f.model,np.array([1.,-.3,.2]),identity='all-ion directional total force'),
                number=regional_number_probe(f.model,np.array([1.,.3]),identity='explicit weighted physical region'))


def compare_readouts(f, boundary, expected, physical):
    from prepared_time_fixture import original_probe
    from substrate_framework.prepared_time_observables import physical_readout
    record = physical_readout(f.model,boundary,physical)
    for name,item in record['probes'].items():
        prefix = name+'_t'+str(boundary.completed_steps)
        compare(prefix+'_mean',item['mean'],original_probe(f,expected[0],name))
        for source_index,connected in enumerate(item['connected']):
            for part_index,part in enumerate(('electronic','ionic','cross')):
                compare(prefix+'_'+str(source_index)+'_'+part,connected[part],
                        original_probe(f,expected[1+source_index*3+part_index],name,connected=True))
        for pair in (0,2):
            actual = item['connected'][pair]['total']-item['connected'][pair+1]['total']
            wanted = sum(original_probe(f,expected[1+pair*3+k]-expected[1+(pair+1)*3+k],name,connected=True) for k in range(3))
            compare(prefix+'_moving_ordered_commutator_'+str(pair),actual,wanted)
    return record


def diagnostics(f, boundary, expected, initial_diagnostics):
    import numpy as np
    from prepared_time_fixture import factor_vector, unpack, original_energy
    from substrate_framework.prepared_time_observables import boundary_diagnostics
    d = boundary_diagnostics(f.model,boundary.mean)
    compare('energy_original_t'+str(boundary.completed_steps),d['energy'],original_energy(f,expected[0]))
    rho,_,_ = unpack(factor_vector(f.model,boundary.mean))
    spectrum = np.array([np.linalg.eigvalsh(x) for x in rho])
    return dict(energy=d['energy'], energy_drift=abs(d['energy']-initial_diagnostics['energy']),
                trace_drift=float(np.max(np.abs(d['traces']-initial_diagnostics['traces']))),
                gram_drift=max(float(np.max(np.abs(a-b))) for a,b in zip(d['grams'],initial_diagnostics['grams'])),
                pauli_eigenvalues=spectrum, pauli_excursion=float(max(0.,-spectrum.min(),spectrum.max()-1.)))


def trajectory(f, sources, expected, steps, contract):
    import numpy as np
    from substrate_framework.prepared_time_evolution import PreparedTimeEvolution
    from substrate_framework.prepared_time_observables import boundary_diagnostics
    identity,limits,reality = contract
    evolution = PreparedTimeEvolution(f.model,f.mean,sources,identity=identity,limits=limits,reality=reality)
    initial_diagnostics = boundary_diagnostics(f.model,f.mean)
    physical = probes(f)
    records, deviations, rows = [], [], []
    checkpoint = None
    for boundary in evolution.stream(dt=.1/steps, steps=steps):
        wanted = expected[boundary.completed_steps*(32//steps)]
        packed = packed_boundary(f,boundary)
        deviations.append(compare('full_original_'+str(steps)+'_'+str(boundary.completed_steps),packed,wanted))
        records.append(dict(readout=compare_readouts(f,boundary,wanted,physical),
                            diagnostics=diagnostics(f,boundary,wanted,initial_diagnostics)))
        rows.append(packed)
        if steps == 16 and boundary.completed_steps == 8:
            checkpoint = json.loads(json.dumps(evolution.checkpoint(),allow_nan=False))
    active = len(sources)*3
    compare('reserved_executed_mean_'+str(steps),list(evolution.counts.values()),[4*steps,4*steps,4*steps*active,4*steps*active],rtol=0,atol=0)
    return evolution, rows, records, deviations, checkpoint


def recovery(f, checkpoint, remaining, contract):
    from substrate_framework.prepared_time_evolution import PreparedTimeEvolution
    identity,limits,reality = contract
    restored = PreparedTimeEvolution.restore(f.model,checkpoint,identity=identity,limits=limits,reality=reality)
    restored_records = []
    physical = probes(f)
    for index,boundary in enumerate(restored.stream(dt=.1/16,steps=8)):
        compare('restored_remaining_'+str(index),packed_boundary(f,boundary),remaining[index],rtol=0,atol=0)
        restored_records.append(compare_readouts(f,boundary,remaining[index],physical))
    return dict(checkpoint=checkpoint,counts=restored.counts,elapsed_seconds=restored.elapsed_seconds,
                remaining_readouts=restored_records,additional_executed_mean=32,additional_executed_tangent=384,
                scope='completed-boundary consumption; not in-flight crash recovery')


def null_parts(f, right):
    from substrate_framework.prepared_covariance_source import PreparedCovariance
    from substrate_framework.prepared_time_observables import physical_readout
    from substrate_framework.prepared_time_evolution import TimeBoundary
    partial = PreparedCovariance(f.model,electronic_authority='known quasifree marginals; ions and cross unknown').source(right)
    item = physical_readout(f.model,TimeBoundary(f.mean,(partial,),0.,0),probes(f))
    for probe in item['probes'].values():
        value = probe['connected'][0]
        if value['ionic'] is not None or value['cross'] is not None or value['total'] is not None or value['complete']:
            raise AssertionError('unknown initial part fabricated')
    CHECKS.append(dict(name='unknown_initial_parts_remain_null',passed=True))


def numerical(bindings):
    import numpy as np
    from prepared_time_fixture import fixture, initial_sources, factor_vector, original_rhs, original_energy
    from substrate_framework.prepared_time_evolution import TimeBoundary
    f = fixture()
    sources,direct,right = initial_sources(f)
    initial = np.vstack((factor_vector(f.model,f.mean),direct))
    compare('literal_Fock_canonical_joint_initial_sources',packed_boundary(f,TimeBoundary(f.mean,sources,0.,0)),initial)
    sensitive('collective_normalization',direct[0]/4,direct[0])
    sensitive('constructive_nonzero_cross',direct[8],np.zeros(24))
    sensitive('connection_transport',original_rhs(f,initial.ravel(),frozen_connection=True),original_rhs(f,initial.ravel()))
    null_parts(f,right)
    expected,ref_receipt = reference(f,initial,bindings['scipy'])
    compare('reference_trace_conservation',np.trace(expected[:,0,:18].reshape(33,2,3,3),axis1=-2,axis2=-1),
            np.broadcast_to(np.sum(f.occupations,axis=1),(33,2)))
    compare('reference_energy_conservation',[original_energy(f,x[0]) for x in expected],
            np.full(33,original_energy(f,initial[0])))
    contract = contracts(f,bindings)
    coarse,rows16,records16,error16,checkpoint = trajectory(f,sources,expected,16,contract)
    fine,rows32,records32,error32,_ = trajectory(f,sources,expected,32,contract)
    sensitive('Euler_instead_of_RK4',initial+.1/16*original_rhs(f,initial.ravel()).reshape(initial.shape),rows16[0])
    coarse_error, fine_error = max(error16),max(error32)
    passed = fine_error < .4*coarse_error+1e-13
    CHECKS.append(dict(name='half_step_RK4_behavior',passed=bool(passed),coarse_error=coarse_error,fine_error=fine_error,
                       error_ratio=coarse_error/fine_error if fine_error else None))
    if not passed:
        raise AssertionError('half-step refinement did not improve the finite trajectory')
    restored = recovery(f,checkpoint,rows16[8:],contract)
    return dict(reference=ref_receipt, interval=[0.,.1], trajectories={'16':records16,'32':records32},
                counts={'16':coarse.counts,'32':fine.counts},
                actual_executed_total=dict(mean=coarse.counts['executed_mean']+fine.counts['executed_mean']+32,
                                           tangent=coarse.counts['executed_tangent']+fine.counts['executed_tangent']+384),
                elapsed_seconds={'16':coarse.elapsed_seconds,'32':fine.elapsed_seconds},
                reality={'16':coarse.reality_diagnostics,'32':fine.reality_diagnostics},
                restoration=restored, contract=dict(identity=contract[0],limits=asdict(contract[1]),reality=asdict(contract[2])))


def json_default(value):
    if isinstance(value,complex):
        return {'real':value.real,'imag':value.imag}
    if hasattr(value,'tolist'):
        return value.tolist()
    raise TypeError(type(value).__name__)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--result',type=Path,required=True)
    args = parser.parse_args()
    if args.result.exists():
        raise FileExistsError('immutable new result path already exists')
    started = time.monotonic()
    result = dict(passed=False,scope='synthetic finite pilot only; no material or long-time accuracy claim')
    try:
        symbolic_first()
        result['source_access'] = source_access()
        result['finite'] = numerical(result['source_access'])
        result['passed'] = True
    except Exception as error:
        result['failure'] = dict(type=type(error).__name__,message=str(error),traceback=traceback.format_exc())
    result.update(checks=CHECKS,elapsed_seconds=time.monotonic()-started)
    with args.result.open('x') as handle:
        json.dump(result,handle,indent=2,default=json_default,allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(result=str(args.result),passed=result['passed'],elapsed_seconds=result['elapsed_seconds'])))
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
