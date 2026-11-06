#!/usr/bin/env python3
"""One original Native source0 readout regrouping; never prepares or evolves.

The existing restore method was inspected but cannot be called against a new
QWork: it attempts retirement in the original namespace. The read-only prepared
constructor is used with all restore identity checks and completed original
retirement authenticated, as explicitly directed by the allocation owner.
"""
import copy
import hashlib
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time
import traceback

START_WALL = time.time()
START_MONO = time.monotonic()
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DIAG = Path('/var/tmp/p254-a28-energy-wide-block-native-pilot-001-n6d4')
REPORT = HERE/'energy-wide-block-Native-pilot-001.json'
BLAS_KEYS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
CAP_MEMORY = 4*1024**3
CAP_OUTPUT = 512*1024**2
INPUTS = {
    'original_admission': dict(path='/var/tmp/p254-a28-saved-energy-positive-001/admission.json', bytes=1546663, sha256='82d2047bde8efea2a548388940cf6e76d7a8e84d8c920f96a534c04e711a649c'),
    'original_state': dict(path='/var/tmp/p254-a28-loaded-raw-mixed-shared-003/steps/000000/stage-3/final.json', bytes=221505, sha256='5f2b72b531b579be3568c1012c2d1f82c94762103b0a88f2ef5c540dbc171c39'),
    'original_prepared': dict(path='/var/tmp/p254-a28-saved-energy-positive-001/work/0d16cac3c09c3bf206c9518f5337e0c2a14bd4b50d49c7049ba55ad6e8874f3c/generation-422f68d35b6740ed84f965e48b5106f1/prepared.json', bytes=658588, sha256='53d6c6aa446950896cfb577d211ae2f04e19b06b1ff7f5d175cdb74261625d37'),
    'original_energy_request': dict(path='/var/tmp/p254-a28-saved-energy-positive-001/energy/request.json', bytes=127916, sha256='57ee646641418b046f210bd234392a6cc30ea008365a5ddb7e166b73b98413d9'),
    'original_source0_complete': dict(path='/var/tmp/p254-a28-saved-energy-positive-001/work/ea0c0042e7a7fc5152efd4d160c7611691f537525735a90e093d00f4c92a54f9/generation-42765db04bde403e8ded6e3b80a56f33/item-0000/complete.json', bytes=740, sha256='3760a7731c344ba001d82bcc484e90404c0ea27c4ed098db96d49bf97d3972c0'),
}


def binding(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    size = 0
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            size += len(block)
            digest.update(block)
    return dict(path=str(path), bytes=size, sha256=digest.hexdigest())


def authenticate_raw(expected):
    actual = binding(expected['path'])
    if actual != dict(path=expected['path'], bytes=int(expected['bytes']), sha256=expected['sha256']):
        raise ValueError('raw byte binding differs: '+expected['path'])
    return actual


def save(path, value):
    data = (json.dumps(value, indent=2, allow_nan=False)+'\n').encode()
    with Path(path).open('xb') as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())
    return binding(path)


def rss(pid):
    try:
        for line in Path('/proc/%s/status' % pid).read_text().splitlines():
            if line.startswith('VmRSS:'):
                return int(line.split()[1])*1024
    except FileNotFoundError:
        pass
    return 0


def worker(started_at):
    os.sched_setaffinity(0, {6})
    sys.dont_write_bytecode = True
    for path in (HERE, HERE.parent/'0025', HERE.parent/'0024', ROOT/'src',
                 Path('/home/administrator/.cache/dft/p254-a28-runtime')):
        sys.path.insert(0, str(path))
    result = dict(status='failed', inputs=INPUTS, invocation_started_unix=started_at,
                  diagnostic_program=binding(__file__), steps=[])
    budget = prepared = state = None
    try:
        if sorted(os.sched_getaffinity(0)) != [6] or any(os.environ.get(k) != '1' for k in BLAS_KEYS):
            raise ValueError('explicit CPU6 and BLASall5=1 allocation required')
        result['authenticated_inputs_before'] = {k: authenticate_raw(v) for k, v in INPUTS.items()}
        result['steps'].append('raw original input bytes authenticated')
        admission = json.loads(Path(INPUTS['original_admission']['path']).read_bytes())
        original = json.loads(Path(INPUTS['original_source0_complete']['path']).read_bytes())
        record = json.loads(Path(INPUTS['original_prepared']['path']).read_bytes())
        original_array = authenticate_raw(original['arrays']['mixed_bare_terms'])
        retirement_path = Path(record['moment_work']['directory'])/'retire.json'
        retirement_binding = binding(retirement_path)
        retirement = json.loads(retirement_path.read_bytes())
        if (retirement.get('schema') != 'q-retirement-v1'
                or retirement.get('barrier') != INPUTS['original_prepared']
                or retirement.get('owner') != record['moment_work']['directory']
                or any(Path(p['path']).exists() for p in retirement['payloads'])):
            raise ValueError('original moment retirement must already be complete, wholly read-only')
        result['original_moment_retirement'] = retirement_binding
        result['original_moment_retired_payload_count'] = len(retirement['payloads'])
        import numpy as np
        import sympy as sp
        from threadpoolctl import threadpool_info, threadpool_limits
        from mixed_energy_saved_runtime import build_saved_runtime
        from mixed_material_driver_support import ContinuationBudget
        from mixed_material_dynamics import NativeMixedPrepared, SCHEMA
        from mixed_response_storage import MappedMixedState
        import mixed_material_energy
        result['sources_before'] = {str(Path(m.__file__).resolve()): binding(m.__file__)
            for m in tuple(sys.modules.values()) if getattr(m, '__file__', None)
            and Path(m.__file__).suffix == '.py' and Path(m.__file__).resolve().is_relative_to(ROOT)}
        budget = ContinuationBudget(DIAG, started_at=started_at, max_actions=64,
            max_seconds=480., max_bytes=CAP_OUTPUT, max_memory_bytes=CAP_MEMORY,
            worker_memory_bytes=CAP_MEMORY)
        with threadpool_limits(limits=1):
            dynamics, work, evidence = build_saved_runtime(admission, budget=budget,
                cpus=[6], work_directory=DIAG/'work')
            state = MappedMixedState(INPUTS['original_state'], budget=budget)
            dynamics._validate_state(state)
            if (record.get('schema') != SCHEMA or record.get('complete') is not True
                    or record.get('identity') != dynamics._identity(state)
                    or dynamics.baseline._saved_frames(state.time) is None
                    or record.get('baseline_frames') is not None):
                raise ValueError('original prepared restore source/state/saved-baseline identity differs')
            prepared = NativeMixedPrepared(dynamics, INPUTS['original_prepared'], record, state)
            result['steps'].append('original saved Native runtime and prepared context restored read-only')
            result['runtime_evidence'] = evidence
            result['prepared_restore'] = dict(method='NativeMixedPrepared with exact restore identity checks',
                reason='restore() retires original moment namespace; no retirement permitted here',
                original_128_moments_repeated=False, prepared_recomputed=False,
                saved_baseline_endpoint=True, baseline_frames=None)
            sample_bindings = {k: state.leg_binding(k, 0) for k in ('plus', 'minus', 'mixed')}
            sample_bindings.update({'prepared-'+k: v for k, v in record['arrays'].items()})
            result['sample_bytes_before'] = {k: authenticate_raw(v) for k, v in sample_bindings.items()}
            result['source_geometry'] = dict(index=0, sources_in_model=dynamics.count,
                slot=dynamics.baseline.initial.slots[0].metadata(),
                raw_coefficient=dynamics.baseline.initial.slots[0].coefficient,
                state_shapes={c: state.record['shapes'][c][0] for c in ('plus','minus','mixed')})
            if (dynamics.block_size != 8 or dynamics.source_identity['block_size'] != 8
                    or any(h.block_size != 8 for h in dynamics.histories.values())):
                raise ValueError('original block8 model required')
            readout = copy.copy(dynamics)
            readout.block_size = 64
            if any(getattr(readout, k) is not v for k, v in vars(dynamics).items() if k != 'block_size'):
                raise ValueError('shallow copy changed more than readout grouping')
            before = budget.snapshot()
            readout_started, readout_cpu = time.monotonic(), time.process_time()
            print('P254_NATIVE_WIDE_N6D4_READOUT_READY', flush=True)
            values = mixed_material_energy._mixed_bare_item(0, state, readout, prepared.stage)['mixed_bare_terms']
            elapsed = time.monotonic()-readout_started
            cpu_seconds = time.process_time()-readout_cpu
            after = budget.snapshot()
            if values.shape != (4,) or values.dtype != np.complex128 or not np.isfinite(values).all():
                raise ValueError('four finite original complex128 dyad terms required')
            output_array = budget.array(DIAG/'source0-wide64-terms.npy', values)
            old = np.load(original_array['path'], allow_pickle=False)
            if old.shape != (4,) or old.dtype != np.complex128:
                raise ValueError('original source0 oracle must have four complex128 dyad terms')
            scale = float(max(np.max(np.abs(old)), np.max(np.abs(values)), np.finfo(float).tiny))
            atol, rtol = 1e-12*scale, 1e-10
            errors = np.abs(values-old)
            thresholds = atol+rtol*np.maximum(np.abs(old), np.abs(values))
            passed = bool(np.all(errors <= thresholds))
            # Independent exact-rational SymPy check of complex squared norms,
            # threshold arithmetic and a discriminating over-tolerance mutation.
            rational = lambda x: sp.Rational(float(x))
            sqnorm = lambda x: rational(x.real)**2+rational(x.imag)**2
            oracle = []
            mutation = values.copy()
            mutation[0] += 100*float(thresholds[0])
            mutation_checks = []
            for a,b,t,m in zip(values,old,thresholds,mutation):
                sq = (rational(a.real)-rational(b.real))**2+(rational(a.imag)-rational(b.imag))**2
                oracle.append(bool(sq <= rational(t)**2))
                msq = (rational(m.real)-rational(b.real))**2+(rational(m.imag)-rational(b.imag))**2
                mutation_checks.append(bool(msq <= rational(t)**2))
            if oracle != (errors <= thresholds).tolist() or mutation_checks[0]:
                raise ValueError('independent comparison oracle or discriminating mutation failed')
            terms = lambda a: [[float(x.real),float(x.imag)] for x in a]
            result['comparison'] = dict(passed=passed, dtype=str(values.dtype),
                term_order=['x D psi.H', 'psi D x.H', 'incoming_plus D incoming_plus.H', 'incoming_minus D incoming_minus.H'],
                original_terms=terms(old), wide64_terms=terms(values),
                absolute_differences=errors.tolist(),
                max_absolute_difference=float(errors.max()),
                relative_differences=(errors/np.maximum(np.abs(old),np.finfo(float).tiny)).tolist(),
                max_relative_difference=float(np.max(errors/np.maximum(np.abs(old),np.finfo(float).tiny))),
                tolerance=dict(rule='abs(new-old) <= 1e-12*global_scale + 1e-10*max(abs(new),abs(old))',
                    global_scale=scale, atol=atol, rtol=rtol, thresholds=thresholds.tolist(),
                    interpretation='declared floating comparison tolerance, not bitwise equality or physical-accuracy certification'),
                independent_sympy=dict(version=sp.__version__, exact_squared_complex_norm_pass=oracle,
                    mutation='add 100 times original term0 tolerance to term0 real part',
                    mutation_term0_rejected=not mutation_checks[0]),
                original_array=original_array, new_array=output_array)
            result['new_readout_scope'] = dict(block_size=64,
                hamiltonian_actions=after['hamiltonian_actions']-before['hamiltonian_actions'],
                action_columns=after['action_columns']-before['action_columns'],
                elapsed_seconds=elapsed, cpu_seconds=cpu_seconds,
                parent_dynamics_block_size=dynamics.block_size,
                original_model_block_size=dynamics.source_identity['block_size'],
                history_block_sizes={k:h.block_size for k,h in dynamics.histories.items()})
            result['historical_source0_scope'] = dict(hamiltonian_actions=original['hamiltonian_actions'],
                elapsed_seconds=original['elapsed_seconds'], cpu_seconds=original['cpu_seconds'],
                process_peak_rss_bytes=original['process_peak_rss_bytes'],
                accounting='old original176 actions are historical; not charged as new or silently credited')
            result['blas'] = dict(environment={k:os.environ[k] for k in BLAS_KEYS}, libraries=threadpool_info())
            result['affinity'] = sorted(os.sched_getaffinity(0))
            result['parent_model_unchanged'] = (dynamics.block_size == 8
                and dynamics.source_identity['block_size'] == 8
                and all(h.block_size == 8 for h in dynamics.histories.values()))
            result['sample_bytes_after'] = {k: authenticate_raw(v) for k,v in sample_bindings.items()}
            result['authenticated_inputs_after'] = {k:authenticate_raw(v) for k,v in INPUTS.items()}
            authenticate_raw(retirement_binding)
            authenticate_raw(original_array)
            for b in result['sources_before'].values():
                authenticate_raw(b)
            result['bytes_unchanged'] = True
            result['steps'].append('genuine original full-PW source0 comparison completed')
            result['status'] = 'complete' if passed and result['parent_model_unchanged'] else 'comparison_failed'
    except BaseException as exc:
        result['error'] = dict(type=type(exc).__name__, message=str(exc), traceback=traceback.format_exc())
        traceback.print_exc()
    finally:
        if prepared is not None:
            prepared.close()
        if state is not None:
            state.close()
        result['budget'] = budget.snapshot() if budget is not None else None
        result['worker_peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
        result['worker_elapsed_seconds'] = time.time()-started_at
        result['limitations'] = ['One original source0 readout regrouping only; not all128 energy sources or full energy completion.',
            'No state or mean evolution, no new DFT, no moment recomputation, no whole-run ETA.',
            'Comparison tolerance is floating readout agreement only, not physical convergence.']
        save(DIAG/'worker-result.json', result)
    return 0 if result['status'] == 'complete' else 1


def supervise():
    os.sched_setaffinity(0, {6})
    DIAG.mkdir(exist_ok=False)
    if REPORT.exists():
        raise FileExistsError(str(REPORT))
    argv = [sys.executable, str(Path(__file__).resolve()), '--worker', repr(START_WALL)]
    launch = dict(argv=argv, supervisor_argv=sys.argv, cwd=os.getcwd(),
        started_unix=START_WALL, cpu_affinity=sorted(os.sched_getaffinity(0)),
        environment={k:os.environ.get(k) for k in BLAS_KEYS}, program=binding(__file__),
        caps=dict(actions=64, seconds_from_invocation_start=480, external_deadline_seconds=540,
            aggregate_rss_bytes=CAP_MEMORY, worker_rss_bytes=CAP_MEMORY, logical_output_bytes=CAP_OUTPUT))
    save(DIAG/'launch.json', launch)
    with (DIAG/'stdout.log').open('xb') as out, (DIAG/'stderr.log').open('xb') as err:
        child = subprocess.Popen(argv, stdout=out, stderr=err, start_new_session=True)
        peak_aggregate = peak_parent = peak_worker = 0
        failure = None
        while child.poll() is None:
            parent_size, worker_size = rss(os.getpid()), rss(child.pid)
            peak_parent, peak_worker = max(peak_parent,parent_size), max(peak_worker,worker_size)
            peak_aggregate = max(peak_aggregate,parent_size+worker_size)
            logical = sum(p.stat().st_size for p in DIAG.rglob('*') if p.is_file())
            elapsed = time.monotonic()-START_MONO
            if elapsed >= 480 or parent_size+worker_size > CAP_MEMORY or worker_size > CAP_MEMORY or logical > CAP_OUTPUT:
                failure = dict(elapsed_seconds=elapsed, parent_rss=parent_size,
                    worker_rss=worker_size, logical_output_bytes=logical,
                    reason='invocation resource limit reached; no broader limits/retry')
                os.killpg(child.pid, signal.SIGKILL)
                break
            time.sleep(.05)
        exit_code = child.wait()
    result = dict(schema='P254-energy-wide-block-Native-pilot-001', status='failed', launch=launch,
        diagnostic_directory=str(DIAG), worker_exit=exit_code, limit_failure=failure,
        supervisor_elapsed_seconds=time.monotonic()-START_MONO,
        supervisor_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
        observed_peak_worker_rss_bytes=peak_worker, observed_peak_aggregate_rss_bytes=peak_aggregate,
        memory_scope='50ms samples: supervisor plus single Native worker; no other simulation workers launched',
        stdout=binding(DIAG/'stdout.log'), stderr=binding(DIAG/'stderr.log'),
        launch_record=binding(DIAG/'launch.json'))
    worker_path = DIAG/'worker-result.json'
    if worker_path.exists():
        result['worker_report'] = binding(worker_path)
        result['result'] = json.loads(worker_path.read_bytes())
        result['status'] = result['result']['status'] if exit_code == 0 and failure is None else 'failed'
    else:
        result['error'] = 'Native worker exited without result; retained exact stdout/stderr and durable action ledger'
    result['logical_diagnostic_output_bytes_before_report'] = sum(p.stat().st_size for p in DIAG.rglob('*') if p.is_file())
    report_binding = save(REPORT, result)
    print(json.dumps(dict(status=result['status'], report=report_binding, worker_exit=exit_code)), flush=True)
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--worker':
        raise SystemExit(worker(float(sys.argv[2])))
    if len(sys.argv) != 1:
        raise ValueError('no user-supplied cap/input overrides accepted')
    raise SystemExit(supervise())
