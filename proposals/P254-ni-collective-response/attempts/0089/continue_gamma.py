"""Advance only the independently ready all31/full114 Gamma Picard block.

No Sigma, electron bank, synthetic provider or full-G2 result is constructed.
Explicit resume skips only durable committed slots; retained uncommitted Native
staging is refused rather than automatically retried.
"""
from pathlib import Path
import argparse
import fcntl
import json
import resource
import sys
from time import perf_counter, time
import uuid
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), *(str(HERE.parent / n) for n in ('0086', '0063', '0062', '0061'))]
from pi_history import PiKernels
from gamma_result import (GammaResult, SCHEMA, validate_inputs, admitted_binary,
                          pi_identity, source_bindings)
from dressed_mean import DressedMean
from history_provider import PreviousHistory
from history_storage import atomic_json, authenticate, binding
from history_readout import gamma_readout, json_value, raw_difference
from production_transport import FinalBank, load_solution
from continue_history import solve_gamma_slot


def _open_bank(path, mean, sources, binary, order, beta):
    resume = path.exists()
    if resume and not (path / 'bank.json').is_file():
        raise ValueError('incomplete Gamma allocation retained; explicit recovery required')
    bank = FinalBank(path, mean.times, 114, 1, statistics=1, order=order,
                     beta=beta, resume=resume)
    try:
        if (bank.meta['shape'] != [1, 31, 31, 114, 114] or
                not np.array_equal(bank.times, mean.times) or
                bank.meta['statistics'] != 1 or bank.meta['order'] != order or
                bank.meta['beta'] != beta or set(bank.completed) - {0}):
            raise ValueError('recovered Gamma bank layout differs')
        for receipt in bank.completed.values():
            if (receipt['provenance']['sources'] != sources or
                    receipt['provenance']['binary'] != binary):
                raise ValueError('committed Gamma slot belongs to another invocation')
            authenticate(receipt['provenance']['native'])
        for stage in (path / '.staging').iterdir():
            if stage.is_dir() and (stage.name != 'slot-000' or 0 not in bank.completed):
                raise ValueError('uncommitted Native staging retained; no automatic replay: ' + str(stage))
        return bank
    except BaseException:
        bank.close()
        raise


def _costs(out):
    attempts = [dict(receipt=binding(path), **json.loads(path.read_bytes()))
                for path in sorted((out / 'attempts').glob('*.json'))]
    cases = []
    for stage in sorted((out / 'gamma' / '.staging').glob('slot-*')):
        path = stage / 'native.json'
        record = {}
        if path.exists():
            try:
                record = json.loads(path.read_bytes())
            except json.JSONDecodeError:
                record = dict(failure='interrupted non-atomic Native receipt write')
        measurements = {}
        timing = stage / 'native.time'
        if timing.exists():
            for line in timing.read_text().splitlines():
                if ': ' in line:
                    key, value = line.strip().rsplit(': ', 1)
                    if key in ('User time (seconds)', 'System time (seconds)',
                               'Maximum resident set size (kbytes)'):
                        measurements[key] = float(value)
        cases.append(dict(bank='gamma', slot=0,
            receipt=binding(path) if path.exists() else None,
            staging=str(stage.resolve()), exit=record.get('exit'),
            failure=record.get('failure'), wall_seconds=record.get('wall_seconds'),
            user_seconds=measurements.get('User time (seconds)'),
            system_seconds=measurements.get('System time (seconds)'),
            maxrss_kib=measurements.get('Maximum resident set size (kbytes)'),
            native_time=binding(timing) if timing.exists() else None))
    unknown = [a['receipt'] for a in attempts if 'elapsed_seconds' not in a]
    incomplete = [c for c in cases if c['wall_seconds'] is None]
    invocation = json.loads((out / 'invocation.json').read_bytes())
    return dict(schema='p254-independent-gamma-cost-v1', native_cases=cases,
        native_cases_recorded=sum(c['receipt'] is not None for c in cases),
        native_cases_successful=sum(c['exit'] == 0 for c in cases),
        native_wall_seconds_recorded=sum(c['wall_seconds'] for c in cases if c['wall_seconds'] is not None),
        native_user_seconds_recorded=sum(c['user_seconds'] for c in cases if c['user_seconds'] is not None),
        native_system_seconds_recorded=sum(c['system_seconds'] for c in cases if c['system_seconds'] is not None),
        native_maxrss_kib=max((c['maxrss_kib'] for c in cases if c['maxrss_kib'] is not None), default=None),
        attempts=attempts,
        owner_elapsed_seconds_recorded=sum(a['elapsed_seconds'] for a in attempts if 'elapsed_seconds' in a),
        incomplete_native_receipts=incomplete, unknown_interrupted_attempt_costs=unknown,
        cumulative_complete=not incomplete and not unknown,
        upstream_input_artifacts=dict(pi=invocation['input_pi'], mean=invocation['input_mean'],
            previous=invocation['input_previous'], binary_verification=invocation['binary_verification']),
        upstream_source_bindings=invocation['source_bindings'],
        cost_scope='Only Gamma continuation costs here. Actual Pi/family/manager, M2, G1 and prior failed/recovery costs remain in bound upstream artifacts; none is reset or asserted zero.',
        accounting='Owner/Native wall times overlap, not additive. Unknown interrupted costs remain unknown. Successful attempt timing ends at raw readout, before final publication.')


def _readout(complete, previous):
    solution = load_solution(Path(complete['path']).parent)
    try:
        ret, less = solution['ret'][0], solution['less'][0]
        old_ret, old_less = previous.gamma_slot()
        covariance = 1j * less[np.arange(31), np.arange(31)]
        n = int(previous.aux['n'])
        return dict(schema='p254-independent-gamma-readout-v1',
            status='complete-gamma-block', full_G2_complete=False,
            times_inverse_ry=solution['times'], gamma= dict(complete=complete),
            probe=dict(name='canonicalGammaNi1/x', index=0, q=0,
                       normalization='1sqrtN', coefficient=float(1 / np.sqrt(n))), n=n,
            intrinsic=gamma_readout(ret, less, n),
            previous=gamma_readout(old_ret, old_less, n),
            original_bare=gamma_readout(previous.aux['source_gamma_ret'],
                previous.aux['source_gamma_less'], n),
            raw_initial_variance=1j * less[0, 0, 0, 0] / n,
            raw_final_variance=1j * less[-1, -1, 0, 0] / n,
            gamma_initial_source_maxabs=float(np.abs(less[0, 0] - previous.aux['source_gamma_less'][0, 0]).max()),
            raw_equal_time_covariance_hermiticity_max=float(np.abs(covariance - covariance.swapaxes(-1, -2).conj()).max()),
            raw_two_time_lesser_antihermiticity_max=max(float(np.abs(less[a, b] + less[b, a].conj().T).max()) for a in range(31) for b in range(31)),
            delta_from_previous=dict(gamma_ret_raw_coordinate_maxabs=raw_difference(ret, old_ret),
                gamma_less_raw_coordinate_maxabs=raw_difference(less, old_less)),
            all_Gamma_position_coordinates_consumed=114,
            interpretation='Unchanged owning A53/history_readout raw probe and normalization. Raw initial/Hermitian defects and iterate differences only; no repair, clipping, calibrated loss, convergence, accuracy, heat or rate. Electrons have not been advanced by this result.')
    finally:
        for key in ('ret', 'less'):
            solution[key]._mmap.close()


def _execute(out, kernels, mean, previous, invocation, binary, order, beta, timeout, resume, started):
    if resume:
        if json.loads((out / 'invocation.json').read_bytes()) != invocation:
            raise ValueError('resume requires identical source/input/binary/order/beta/timeout bindings')
        if (out / 'result.json').exists() and (out / 'handoff.json').exists():
            handoff = json.loads((out / 'handoff.json').read_bytes())
            if handoff['result'] != binding(out / 'result.json'):
                raise ValueError('Gamma handoff/result binding differs')
            with GammaResult(out / 'result.json') as saved:
                saved.validate(kernels, mean, previous, binary=binary,
                    verification=invocation['binary_verification'], order=order, beta=beta)
                return saved.spec
    else:
        atomic_json(out / 'invocation.json', invocation)
    (out / 'attempts').mkdir(exist_ok=True)
    attempt_path = out / 'attempts' / (uuid.uuid4().hex + '.json')
    attempt = dict(schema='p254-independent-gamma-attempt-v1', status='running',
        started_unix=time(), resumed=resume, committed_slots_skipped=0,
        newly_attempted_slots=0, uncommitted_slot=None)
    atomic_json(attempt_path, attempt)
    try:
        with _open_bank(out / 'gamma', mean, invocation['source_bindings'], binary, order, beta) as bank:
            attempt['committed_slots_skipped'] = len(bank.completed)
            atomic_json(attempt_path, attempt)
            if 0 not in bank.completed:
                attempt.update(newly_attempted_slots=1, uncommitted_slot=dict(bank='gamma', slot=0))
                atomic_json(attempt_path, attempt)
                solve_gamma_slot(mean, kernels, bank, binary=binary['path'], order=order,
                    beta=beta, source_bindings=invocation['source_bindings'],
                    previous=previous, timeout=timeout)
                attempt['uncommitted_slot'] = None
                atomic_json(attempt_path, attempt)
            if not (bank.path / 'complete.json').exists():
                bank.finalize()
        complete = binding(out / 'gamma' / 'complete.json')
        readout = _readout(complete, previous)
        atomic_json(out / 'readout.json', json.loads(json.dumps(readout, default=json_value, allow_nan=False)))
        attempt.update(status='complete', elapsed_seconds=perf_counter() - started,
            owner_maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        atomic_json(attempt_path, attempt)
        costs = _costs(out)
        atomic_json(out / 'costs.json', costs)
        result = dict(schema=SCHEMA, status='complete-gamma-block', full_G2_complete=False,
            gamma_slots=1, gamma_coordinates=114, n=64, times_inverse_ry=mean.times.tolist(),
            gamma=dict(complete=complete), readout=binding(out / 'readout.json'),
            invocation=binding(out / 'invocation.json'), costs=binding(out / 'costs.json'),
            input_identity=invocation['input_identity'], input_pi=invocation['input_pi'],
            pi_arrays=invocation['pi_arrays'], input_mean=mean.binding,
            input_previous=previous.binding, preparation_origin=mean.preparation_origin,
            binary=binary, binary_verification=invocation['binary_verification'], order=order, beta=beta,
            source_bindings=invocation['source_bindings'],
            native_cases_recorded=costs['native_cases_recorded'],
            limits=['Only one actual Gamma Jacobi/Picard block, not full G2 or an electronic provider.',
                    'No calibrated loss, accuracy, convergence, conservation, heat or nuclear rate claim.'])
        atomic_json(out / 'result.json', result)
        atomic_json(out / 'handoff.json', dict(schema='p254-independent-gamma-handoff-v1',
            status='complete-gamma-block', full_G2_complete=False,
            result=binding(out / 'result.json'), gamma=result['gamma'], readout=result['readout'],
            costs=result['costs'], input_identity=result['input_identity'],
            input_pi=result['input_pi'], pi_arrays=result['pi_arrays'], input_mean=mean.binding,
            input_previous=previous.binding, binary=binary,
            binary_verification=result['binary_verification'], source_bindings=result['source_bindings'],
            public_reader='0089.gamma_result.GammaResult',
            reuse='Validate exact inputs and all three Pi NPY bytes/axes/dtype against the actual complete fullcatalog before any electronic Native; reuse this external complete binding without copying or recomputing.'))
        return result
    except BaseException as exc:
        attempt.update(status='failed', elapsed_seconds=perf_counter() - started,
            owner_maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            failure_type=type(exc).__name__, failure=str(exc))
        atomic_json(attempt_path, attempt)
        atomic_json(out / 'costs.json', _costs(out))
        raise


def continue_gamma(pi_path, mean_path, previous_path, output, *, binary, order=1,
                   beta=1., timeout=900., resume=False):
    started = perf_counter()
    out = Path(output).resolve()
    with PiKernels(pi_path) as kernels, DressedMean(mean_path) as mean, PreviousHistory(previous_path) as previous:
        identity = validate_inputs(kernels.catalog, mean, previous)
        if order != 1 or beta != 1. or not np.isfinite(timeout) or timeout <= 0:
            raise ValueError('original order1/beta1 and positive finite timeout required')
        pi_arrays = pi_identity(kernels)
        binary_binding, verification = admitted_binary(binary)
        files = [Path(__file__), HERE / 'gamma_result.py', HERE / 'pi_history.py',
            HERE.parent / '0086' / 'dressed_mean.py', HERE.parent / '0063' / 'continue_history.py',
            HERE.parent / '0063' / 'history_provider.py', HERE.parent / '0063' / 'history_storage.py',
            HERE.parent / '0063' / 'saved_mean.py', HERE.parent / '0063' / 'history_readout.py',
            HERE.parent / '0053' / 'consume_material_step.py',
            HERE.parent / '0061' / 'production_transport.py', HERE.parent / '0049' / 'command_receipts.py']
        sources = source_bindings([binding(pi_path), kernels.bindings, pi_arrays,
            mean.bindings, previous.provider.bindings, previous.binding, binary_binding,
            verification, [binding(path) for path in files]])
        for item in sources:
            authenticate(item)
        invocation = dict(schema='p254-independent-gamma-invocation-v1',
            input_identity=identity, input_pi=binding(pi_path), pi_arrays=pi_arrays,
            input_mean=mean.binding, input_previous=previous.binding,
            binary=binary_binding, binary_verification=verification, source_bindings=sources,
            times_inverse_ry=mean.times.tolist(), n=64, order=order, beta=beta, timeout_seconds=timeout,
            layout='A61 Gamma1x114 all31 only; no electronic allocation or fullG2 provider',
            recovery='Explicit resume skips committed slot only; uncommitted staging retained/refused, not an in-Native checkpoint or retry promise.')
        if resume:
            if not (out / 'invocation.json').is_file():
                raise ValueError('resume requires persisted Gamma invocation')
        else:
            out.mkdir(parents=True, exist_ok=False)
        with (out / 'continuation.lock').open('a+b') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            return _execute(out, kernels, mean, previous, invocation, binary_binding,
                            order, beta, timeout, resume, started)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('pi', 'mean', 'previous', 'output', 'binary'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--order', type=int, default=1)
    parser.add_argument('--beta', type=float, default=1.)
    parser.add_argument('--timeout', type=float, default=900.)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    print(json.dumps(continue_gamma(args.pi, args.mean, args.previous, args.output,
        binary=args.binary, order=args.order, beta=args.beta, timeout=args.timeout,
        resume=args.resume)), flush=True)


if __name__ == '__main__':
    main()
