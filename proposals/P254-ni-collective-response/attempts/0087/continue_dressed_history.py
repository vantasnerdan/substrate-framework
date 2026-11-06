"""Advance complete actual dressed kernels with M2, preserving raw history.

This is the A63 slot solve and public SavedHistory/readout path with distinct
actual-provider and original-numerical-reference identities. --resume reopens
A61 durable slot receipts and never reruns committed slots. An uncommitted
staging directory fails closed: no automatic unchanged retry or budget reset.
"""
from pathlib import Path
from contextlib import ExitStack, nullcontext
import argparse
import fcntl
import json
import resource
import sys
from time import perf_counter, time
import uuid

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent / n) for n in ('0089', '0086', '0063', '0062', '0061')]
from dressed_mean import DressedMean
from history_storage import (SOURCE_AUXILIARY, CONTROL_AUXILIARY, binding,
                             atomic_json, authenticate, finite_blocks)
from history_provider import PreviousHistory, SavedHistory, SCHEMA as PROVIDER_SCHEMA
from history_readout import consume, json_value, projected_control_readout
from saved_mean import POSITIONS
from history_banks import OwnedKernels, REFERENCE, PI_CONVENTION
from production_transport import FinalBank
from continue_history import solve_electron_slot, solve_gamma_slot
from gamma_result import GammaResult

INVOCATION_SCHEMA = 'p254-actual-dressed-history-continuation-invocation-v1'
RESULT_SCHEMA = 'p254-actual-dressed-history-continuation-result-v1'


def validate_reference_identity(catalog, mean, previous):
    if previous.provider is None:
        raise ValueError('actual previous dressed provider required, not bare/legacy/control data')
    numerical = previous.provider.spec['input_identity']['numerical_reference']
    if (catalog['input_dataset'] != previous.binding or
            mean.input_identity['actual_dataset'] != previous.binding or
            catalog.get('numerical_reference') != numerical or
            mean.reference_identity != numerical or
            catalog.get('preparation_origin') != mean.preparation_origin):
        raise ValueError('actual previous G1 binding or separate original numerical reference differs')
    authenticate(numerical)
    authenticate(mean.preparation_origin)
    # DressedMean authenticates the complete original-preparation ancestry and
    # actual raw G1 forcing. PreviousHistory.reference_identity is deliberately
    # NOT used here: for a provider it denotes provider bytes, not bundle bytes.
    return dict(kernel_input=catalog['input_dataset'], actual_dataset=previous.binding,
        numerical_reference=numerical, preparation_origin=mean.preparation_origin,
        source_artifact=mean.binding, previous_artifact=previous.binding)


def _binary_binding(binary):
    verification_path = HERE.parent / '0085' / 'verification-001.json'
    verification = json.loads(verification_path.read_bytes())
    actual = binding(binary)
    approved = verification['production_binary']
    if any(actual[key] != approved[key] for key in ('bytes', 'sha256')):
        raise ValueError('the admitted A85 real-axis result binary is required')
    return actual, binding(verification_path)


def _open_bank(path, mean, dimension, slots, statistics, order, beta, sources, binary):
    resume = path.exists()
    if resume and not (path / 'bank.json').exists():
        raise ValueError('incomplete bank allocation retained; inspect before explicit recovery: ' + str(path))
    bank = FinalBank(path, mean.times, dimension, slots, statistics=statistics,
                     order=order, beta=beta, resume=resume)
    try:
        if (not np.array_equal(bank.times, mean.times) or
                bank.meta['shape'] != [slots, 31, 31, dimension, dimension] or
                bank.meta['statistics'] != statistics or bank.meta['order'] != order or
                bank.meta['beta'] != beta):
            raise ValueError('recovered slot bank differs from requested all31 layout')
        for receipt in bank.completed.values():
            provenance = receipt['provenance']
            if provenance['sources'] != sources or provenance['binary'] != binary:
                raise ValueError('committed slot belongs to another actual continuation')
            authenticate(provenance['native'])
        # A61 retains failed/in-flight staging and does not silently overwrite it.
        # Check every uncommitted staging directory before running any fresh slot.
        for stage in (path / '.staging').iterdir():
            if stage.is_dir():
                try:
                    slot = int(stage.name.removeprefix('slot-'))
                except ValueError as exc:
                    raise ValueError('unknown retained staging: ' + str(stage)) from exc
                if slot not in bank.completed:
                    raise ValueError('uncommitted Native staging retained; no automatic replay: ' + str(stage))
        return bank
    except BaseException:
        bank.close()
        raise


def _native_costs(out):
    cases = []
    for bank in ('electron', 'gamma'):
        for path in sorted((out / bank / '.staging').glob('slot-*/native.json')):
            try:
                record = json.loads(path.read_bytes())
            except json.JSONDecodeError:
                record = dict(failure='interrupted non-atomic command-receipt write')
            measurements = {}
            time_path = path.with_name('native.time')
            if time_path.exists():
                for line in time_path.read_text().splitlines():
                    if ': ' in line:
                        key, value = line.strip().rsplit(': ', 1)
                        if key in ('User time (seconds)', 'System time (seconds)',
                                   'Maximum resident set size (kbytes)'):
                            measurements[key] = float(value)
            cases.append(dict(bank=bank, slot=int(path.parent.name.removeprefix('slot-')),
                receipt=binding(path), exit=record.get('exit'),
                failure=record.get('failure'),
                wall_seconds=record.get('wall_seconds'),
                user_seconds=measurements.get('User time (seconds)'),
                system_seconds=measurements.get('System time (seconds)'),
                maxrss_kib=measurements.get('Maximum resident set size (kbytes)'),
                native_time=binding(time_path) if time_path.exists() else None))
    return cases


def _costs(out, precomputed_gamma=None):
    attempts = []
    for path in sorted((out / 'attempts').glob('*.json')):
        record = json.loads(path.read_bytes())
        attempts.append(dict(receipt=binding(path), **record))
    cases = _native_costs(out)
    unfinished = [case for case in cases if case['wall_seconds'] is None]
    unknown_attempts = [a['receipt'] for a in attempts if 'elapsed_seconds' not in a]
    return dict(schema='p254-dressed-continuation-cost-v1', native_cases=cases,
        precomputed_gamma=None if precomputed_gamma is None else dict(
            result=precomputed_gamma.binding, complete=precomputed_gamma.complete,
            costs=precomputed_gamma.spec['costs'],
            original_accounting=json.loads(authenticate(precomputed_gamma.spec['costs']).read_bytes())),
        new_native_scope='All local Native receipts only; precomputed Gamma and its upstream costs are separate, not new electronic calls or zero-cost work.',
        native_cases_recorded=len(cases),
        native_cases_successful=sum(case['exit'] == 0 for case in cases),
        native_wall_seconds_recorded=sum(case['wall_seconds'] for case in cases
                                         if case['wall_seconds'] is not None),
        native_user_seconds_recorded=sum(case['user_seconds'] for case in cases
                                         if case['user_seconds'] is not None),
        native_system_seconds_recorded=sum(case['system_seconds'] for case in cases
                                           if case['system_seconds'] is not None),
        native_maxrss_kib=max((case['maxrss_kib'] for case in cases
                              if case['maxrss_kib'] is not None), default=None),
        attempts=attempts,
        owner_elapsed_seconds_recorded=sum(a['elapsed_seconds'] for a in attempts
                                          if 'elapsed_seconds' in a),
        incomplete_native_receipts=[case['receipt'] for case in unfinished],
        unknown_interrupted_attempt_costs=unknown_attempts,
        cumulative_complete=not unfinished and not unknown_attempts,
        cost_scope='This continuation output only; upstream G1/M2/kernel production costs remain in their bound artifacts, not silently reset or counted as zero.',
        accounting='Owner attempt wall times and Native case wall/CPU are overlapping measures, not additive. Completed-attempt timing ends at the readout boundary; final JSON publication is excluded. Unknown interrupted costs are explicit, never zero.')


def _unique_bindings(items):
    found = {}
    for item in items:
        item = {key: item[key] for key in ('path', 'bytes', 'sha256')}
        path = item['path']
        if path in found and found[path] != item:
            raise ValueError('conflicting source binding: ' + path)
        found[path] = item
    return list(found.values())


def _execute(out, kernels, mean, previous, kernels_path, identities, sources,
             invocation, binary, order, beta, timeout, projected_control, resume, started,
             precomputed_gamma):
    if resume:
        existing = json.loads((out / 'invocation.json').read_bytes())
        if existing != invocation:
            raise ValueError('resume requires identical inputs, binary, sources, order, beta and timeout')
        if (out / 'result.json').exists() and (out / 'handoff.json').exists():
            result = json.loads((out / 'result.json').read_bytes())
            handoff = json.loads((out / 'handoff.json').read_bytes())
            if handoff['result'] != binding(out / 'result.json'):
                raise ValueError('handoff and complete result binding differ')
            for name in ('provider', 'readout', 'invocation', 'costs'):
                authenticate(result[name])
            if result['status'] != 'complete':
                raise ValueError('unexpected noncomplete result artifact')
            # Reopen the actual public provider even for an idempotent call.
            with SavedHistory(result['provider']['path']):
                pass
            return result
    else:
        atomic_json(out / 'invocation.json', invocation)
    (out / 'attempts').mkdir(exist_ok=True)
    attempt_path = out / 'attempts' / (uuid.uuid4().hex + '.json')
    attempt = dict(schema='p254-dressed-continuation-attempt-v1', status='running',
        started_unix=time(), resumed=resume, committed_slots_skipped=0,
        newly_attempted_slots=0, uncommitted_slot=None)
    atomic_json(attempt_path, attempt)
    try:
        # Precomputed Gamma was authenticated before entering this execution.
        # Otherwise validate both local banks before starting electronic work.
        with ExitStack() as banks:
            gamma = None if precomputed_gamma is not None else banks.enter_context(
                _open_bank(out / 'gamma', mean, 114, 1, 1, order, beta, sources, binary))
            electron = banks.enter_context(
                _open_bank(out / 'electron', mean, 12, 128, -1, order, beta, sources, binary))
            attempt['committed_slots_skipped'] = len(electron.completed) + (len(gamma.completed) if gamma is not None else 0)
            attempt['precomputed_gamma_reused'] = precomputed_gamma is not None
            atomic_json(attempt_path, attempt)
            for slot in range(128):
                if slot in electron.completed:
                    continue
                attempt['newly_attempted_slots'] += 1
                attempt['uncommitted_slot'] = dict(bank='electron', slot=slot)
                atomic_json(attempt_path, attempt)
                solve_electron_slot(mean, kernels, electron, slot,
                    binary=binary['path'], order=order, beta=beta,
                    source_bindings=sources, previous=previous, timeout=timeout)
                attempt['uncommitted_slot'] = None
                atomic_json(attempt_path, attempt)
            if not (electron.path / 'complete.json').exists():
                electron.finalize()
            if gamma is not None:
                if 0 not in gamma.completed:
                    attempt['newly_attempted_slots'] += 1
                    attempt['uncommitted_slot'] = dict(bank='gamma', slot=0)
                    atomic_json(attempt_path, attempt)
                    solve_gamma_slot(mean, kernels, gamma, binary=binary['path'],
                        order=order, beta=beta, source_bindings=sources,
                        previous=previous, timeout=timeout)
                    attempt['uncommitted_slot'] = None
                    atomic_json(attempt_path, attempt)
                if not (gamma.path / 'complete.json').exists():
                    gamma.finalize()
        gamma_complete = (precomputed_gamma.complete if precomputed_gamma is not None
                          else binding(out / 'gamma' / 'complete.json'))
        old = previous.provider.spec
        auxiliary = {name: old['auxiliary'][name] for name in SOURCE_AUXILIARY}
        controls = old.get('represented_control')
        if controls is None and previous.control is not None:
            controls = {name: old['auxiliary'][name] for name in CONTROL_AUXILIARY}
        projected = projected_control_readout(previous, projected_control)
        provider = dict(schema=PROVIDER_SCHEMA, status='complete',
            times_inverse_ry=mean.times.tolist(),
            electron=dict(complete=binding(out / 'electron' / 'complete.json')),
            gamma=dict(complete=gamma_complete),
            input_gamma_result=precomputed_gamma.binding if precomputed_gamma is not None else None,
            gamma_origin='precomputed-independent-Gamma-block' if precomputed_gamma is not None else 'new-local-Gamma-slot',
            auxiliary=auxiliary, represented_control=controls,
            projected_control=projected.get('control'), input_identity=identities,
            numerical_reference=identities['numerical_reference'],
            preparation_origin=mean.preparation_origin,
            original_bare_reference=mean.preparation_origin,
            sourcekind='single-physical-Hermitian-contour', integration_order=order,
            beta=beta, reference_source_kind=REFERENCE, pi_convention=PI_CONVENTION,
            input_kernels=binding(kernels_path), input_mean=mean.binding,
            input_previous=previous.binding, source_bindings=sources,
            prescribed_nonGamma_bath=old['prescribed_nonGamma_bath'],
            interpretation='Actual dressed common-functional Picard iterate from previous provider and M2 Gaussian Dyson source, not a refill of G1, fixed-point, temporal accuracy or conservation certificate.')
        atomic_json(out / 'provider.json', provider)
        with SavedHistory(out / 'provider.json') as history:
            readout = consume(history, previous, projected_control=projected_control)
        # Preserve all raw complex observables and defects using the unchanged
        # public readout serializer, with no projection, clipping or Hermitization.
        serialized_readout = json.loads(json.dumps(readout, default=json_value, allow_nan=False))
        atomic_json(out / 'readout.json', serialized_readout)
        attempt.update(status='complete', elapsed_seconds=perf_counter() - started,
            owner_maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        atomic_json(attempt_path, attempt)
        costs = _costs(out, precomputed_gamma)
        atomic_json(out / 'costs.json', costs)
        result = dict(schema=RESULT_SCHEMA, status='complete',
            provider=binding(out / 'provider.json'), readout=binding(out / 'readout.json'),
            invocation=binding(out / 'invocation.json'), costs=binding(out / 'costs.json'),
            input_identity=identities, times_inverse_ry=mean.times.tolist(),
            all128_electron_slots_consumed=True, full114_new_Gamma_consumed=True,
            gamma_complete=gamma_complete,
            precomputed_gamma=precomputed_gamma.binding if precomputed_gamma is not None else None,
            new_native_scope='128 electronic slots only when Gamma is precomputed; otherwise 128 electronic plus one Gamma, less committed slots skipped.',
            native_cases_recorded=costs['native_cases_recorded'],
            case_elapsed_seconds_through_readout=attempt['elapsed_seconds'],
            cumulative_owner_elapsed_seconds_recorded=costs['owner_elapsed_seconds_recorded'],
            cumulative_native_wall_seconds_recorded=costs['native_wall_seconds_recorded'],
            owner_maxrss_kib=attempt['owner_maxrss_kib'],
            delta_from_previous=serialized_readout['delta_from_previous'],
            projected_comparison=serialized_readout['projected_comparison'],
            limits=['Raw defects and iterate differences are retained, not truncation/error/convergence or conservation certificates.',
                    'No clipping, refill, independent reverse-leg inference, control-schema cast, physical calibration, loss, heat or nuclear rate.',
                    'Prescribed nonGamma bath remains the original preparation; only the actual supplied common map and M2 source are advanced.'])
        atomic_json(out / 'result.json', result)
        atomic_json(out / 'handoff.json', dict(
            schema='p254-actual-dressed-history-handoff-v1', status='complete',
            result=binding(out / 'result.json'), provider=result['provider'],
            readout=result['readout'], costs=result['costs'],
            gamma_complete=gamma_complete,
            precomputed_gamma=precomputed_gamma.binding if precomputed_gamma is not None else None,
            input_identity=identities, preparation_origin=mean.preparation_origin,
            actual_next_dataset=result['provider'], public_reader='0063.history_provider.SavedHistory',
            times_inverse_ry=mean.times.tolist(), slots=128, bands=12,
            gamma_coordinates=114, nonclaims=result['limits']))
        return result
    except BaseException as exc:
        attempt.update(status='failed', elapsed_seconds=perf_counter() - started,
            owner_maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            failure_type=type(exc).__name__, failure=str(exc))
        atomic_json(attempt_path, attempt)
        atomic_json(out / 'costs.json', _costs(out, precomputed_gamma))
        raise


def continue_history(kernels_path, mean_path, previous_path, output, *, binary,
                     order, beta=1., timeout=600., resume=False, projected_control=None,
                     gamma_result=None):
    started = perf_counter()
    out = Path(output).resolve()
    with OwnedKernels(kernels_path, require_complete=True) as kernels, \
            DressedMean(mean_path) as mean, PreviousHistory(previous_path) as previous, \
            (GammaResult(gamma_result) if gamma_result is not None else nullcontext()) as precomputed_gamma:
        identities = validate_reference_identity(kernels.catalog, mean, previous)
        if (kernels.times.tolist() != list(range(31)) or
                not np.array_equal(mean.times, kernels.times) or
                not np.array_equal(previous.times, kernels.times) or
                order not in range(1, 6) or not np.isfinite(beta) or beta <= 0 or
                not np.isfinite(timeout) or timeout <= 0):
            raise ValueError('full exact31/all128 source and admitted order/beta/timeout required')
        for descriptor in kernels.catalog.get('source_bindings', []):
            authenticate(descriptor)
        if not all(finite_blocks(array) for bank in (kernels.sigma, kernels.pi)
                   for array in bank.values()):
            raise ValueError('complete finite actual common kernels required')
        if (not np.array_equal(1j * mean.rho0, previous.aux['source_electron_less0']) or
                not np.array_equal(-1j * mean.w0[np.ix_(POSITIONS, POSITIONS)].T,
                                   previous.aux['source_gamma_less'][0, 0]) or
                int(mean.n) != int(previous.aux['n'])):
            raise ValueError('original all128/full114 preparation or exact normalization changed')
        projected = projected_control_readout(previous, projected_control)
        binary_binding, verification = _binary_binding(binary)
        if precomputed_gamma is not None:
            precomputed_gamma.validate(kernels, mean, previous, binary=binary_binding,
                verification=verification, order=order, beta=beta)
        source_files = [Path(__file__), HERE / 'bank_storage.py', HERE / 'functional_bank.py',
            HERE.parent / '0063' / 'continue_history.py',
            HERE.parent / '0063' / 'history_provider.py',
            HERE.parent / '0063' / 'history_readout.py',
            HERE.parent / '0063' / 'history_storage.py',
            HERE.parent / '0063' / 'saved_mean.py',
            HERE.parent / '0062' / 'history_banks.py',
            HERE.parent / '0061' / 'production_transport.py',
            HERE.parent / '0049' / 'command_receipts.py',
            HERE.parent / '0086' / 'dressed_mean.py']
        source_files.append(HERE.parent / '0089' / 'gamma_result.py')
        sources = _unique_bindings([binding(kernels_path), mean.binding,
            previous.binding, binary_binding, verification, *mean.bindings,
            *previous.provider.bindings, *kernels.catalog.get('source_bindings', []),
            *(precomputed_gamma.bindings if precomputed_gamma is not None else []),
            *(binding(path) for path in source_files)])
        invocation = dict(schema=INVOCATION_SCHEMA, source_bindings=sources,
            input_identity=identities, times_inverse_ry=mean.times.tolist(),
            order=order, beta=beta, timeout_seconds=timeout,
            layout='A61 electron128x12 plus Gamma1x114; full31 history',
            sourcekind='single-physical-Hermitian-contour',
            reference_source_kind=REFERENCE, pi_convention=PI_CONVENTION,
            projected_control=projected.get('control'),
            gamma_result=precomputed_gamma.binding if precomputed_gamma is not None else None,
            gamma_complete=precomputed_gamma.complete if precomputed_gamma is not None else None,
            recovery='Committed-slot receipts reopen and skip science; uncommitted staging retained and refused, no automatic retry.')
        if resume:
            if not (out / 'invocation.json').is_file():
                raise ValueError('resume requires a persisted continuation invocation')
        else:
            out.mkdir(parents=True, exist_ok=False)
        with (out / 'continuation.lock').open('a+b') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            return _execute(out, kernels, mean, previous, kernels_path, identities,
                sources, invocation, binary_binding, order, beta, timeout,
                projected_control, resume, started, precomputed_gamma)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('kernels', 'mean', 'previous', 'output', 'binary'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--order', type=int, required=True)
    parser.add_argument('--beta', type=float, default=1.)
    parser.add_argument('--timeout', type=float, default=600.)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--gamma-result', type=Path,
        help='Actual complete independent Gamma result; validate exact inputs before electronic Native and reuse its external bank unchanged.')
    parser.add_argument('--projected-control', type=Path,
        help='Optional actual control for endpoint readout only, never a Dyson source.')
    args = parser.parse_args()
    result = continue_history(args.kernels, args.mean, args.previous, args.output,
        binary=args.binary, order=args.order, beta=args.beta, timeout=args.timeout,
        resume=args.resume, projected_control=args.projected_control,
        gamma_result=args.gamma_result)
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
