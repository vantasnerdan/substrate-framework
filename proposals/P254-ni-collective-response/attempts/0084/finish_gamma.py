"""Recover only the original A84 Gamma solve, then perform the original A63 readout.

The packed input is authenticated and reused in place. No producer, electron
solver, mean integration, control solve, or kernel assembly is called here.
Root supplies the allocation, equivalent binary and finite-oracle provenance.
"""
from pathlib import Path
import argparse
from contextlib import ExitStack
import json
import math
import os
import re
import resource
import sys
from time import perf_counter, time_ns

HERE = Path(__file__).resolve().parent
for directory in ('0061', '0063'):
    sys.path.insert(0, str(HERE.parent / directory))
import numpy as np
from production_transport import (FinalBank, ResultView, HEADER, SOURCEKIND,
    binding, run_recorded, atomic_json, load_solution, plan_layout)
from history_storage import (SOURCE_AUXILIARY, CONTROL_AUXILIARY,
    authenticate, save_array)
from history_provider import PreviousHistory, SavedHistory, SCHEMA, close_solution
from saved_mean import SavedMean, REFERENCE, PI_CONVENTION, POSITIONS
from history_readout import consume, json_value, projected_control_readout

RECOVERY_SCHEMA = 'p254-a84-gamma-only-recovery-v1'
PAID_WALL = 638.905548795592
CASE_BUDGET = 1800.


def read_json(path):
    return json.loads(Path(path).read_bytes())


def positive(value, name, *, zero=False):
    if not math.isfinite(value) or (value < 0 if zero else value <= 0):
        raise ValueError(name + ' must be finite and ' + ('nonnegative' if zero else 'positive'))
    return value


def authenticated_json(item):
    return read_json(authenticate(item))


def allocation(args):
    """Check the actual worker allocation, not just caller labels."""
    if (args.cpu != 6 or args.memory_bytes != 8589934592 or args.swap_bytes != 0
            or args.threads != 1):
        raise ValueError('A85 Gamma allocation requires CPU6, 8GiB, swap0, one thread')
    relative = Path('/proc/self/cgroup').read_text().strip().split('::')[-1]
    group = Path('/sys/fs/cgroup') / relative.lstrip('/')
    current = dict(affinity=sorted(os.sched_getaffinity(0)),
        memory_max=(group/'memory.max').read_text().strip(),
        swap_max=(group/'memory.swap.max').read_text().strip(),
        threads={name: os.environ.get(name) for name in
                 ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
                  'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS')})
    if (current['affinity'] != [args.cpu] or
            current['memory_max'] != str(args.memory_bytes) or
            current['swap_max'] != str(args.swap_bytes) or
            any(value != str(args.threads) for value in current['threads'].values())):
        raise ValueError('actual recovery worker allocation differs from declared allocation')
    return dict(cpu=args.cpu, memory_bytes=args.memory_bytes,
        swap_bytes=args.swap_bytes, threads=args.threads, restart='no', actual=current,
        provenance=binding(args.allocation_record))


def validate_bank(bank, times, dimension, slots, statistics, order, beta):
    if (bank['ret'].shape != (slots, len(times), len(times), dimension, dimension)
            or bank['less'].shape != bank['ret'].shape
            or bank['ret'].dtype != np.dtype('<c16')
            or bank['less'].dtype != np.dtype('<c16')
            or not np.array_equal(bank['times'], times)
            or bank['statistics'] != statistics or bank['order'] != order
            or bank['beta'] != beta):
        raise ValueError('bank shape/order/beta/statistics/times mismatch')


def validate_electrons(out, times, order, beta, expected):
    authenticate(expected)
    bank = load_solution(out/'electron')  # Existing authenticated, read-only reader.
    try:
        validate_bank(bank, times, 12, 128, -1, order, beta)
    finally:
        close_solution(bank)
    complete = read_json(out/'electron'/'complete.json')
    if len(complete['slots']) != 128:
        raise ValueError('exactly 128 completed electron sectors required')
    receipts, wall = [], 0.
    for slot, item in enumerate(complete['slots']):
        receipt = authenticated_json(item)
        if receipt['slot'] != slot:
            raise ValueError('electron committed slot order differs')
        native = authenticated_json(receipt['provenance']['native'])
        if native.get('exit') != 0 or native.get('failure'):
            raise ValueError('committed electron lacks a successful Native receipt')
        wall += positive(native['wall_seconds'], 'electron native wall', zero=True)
        receipts.append(item)
    return dict(complete=expected, committed_slots=receipts, native_wall_seconds=wall)


def packed_header(packed, invocation, original, args):
    with packed.open('rb') as stream:
        header = stream.read(HEADER.size)
        magic, nt, dimension, order, statistics, beta = HEADER.unpack(header)
        if (magic != b'P254VIE1' or nt != 30 or dimension != 114 or statistics != 1
                or order != original['order'] or beta != original['beta']):
            raise ValueError('original full31 Gamma header differs')
        times = np.frombuffer(stream.read((nt+1)*8), dtype='<f8').copy()
    if (not np.array_equal(times, np.asarray(original['times_inverse_ry'], dtype='<f8'))
            or not np.array_equal(times, np.arange(31, dtype='<f8'))):
        raise ValueError('original full31 Gamma times differ')
    layout = plan_layout(len(times), dimension,
        max_input_bytes=args.max_input_bytes, max_result_bytes=args.max_result_bytes,
        max_bank_bytes=args.max_bank_bytes)
    if invocation['layout'] != layout or packed.stat().st_size != layout['input_bytes']:
        raise ValueError('original packed input extent/layout differs')
    return header, times, order, beta, layout


def cost_snapshot(args, stage, started, electrons, old_native):
    # Native wall is a subset of owner wall; never add it a second time.
    earlier = sum(read_json(p)['elapsed_seconds'] for p in stage.glob('execution-*.json'))
    current = perf_counter()-started
    new_native = read_json(stage/'native.json') if (stage/'native.json').exists() else None
    native_wall = 0. if new_native is None else new_native.get('wall_seconds', 0.)
    recovery_wall = earlier + current
    total = args.prior_paid_wall_seconds + args.additional_paid_wall_seconds + recovery_wall
    return dict(prior_case_paid_wall_seconds=args.prior_paid_wall_seconds,
        prior_manager=binding(args.prior_manager_receipt),
        additional_paid_wall_seconds=args.additional_paid_wall_seconds,
        additional_cost_provenance=[binding(p) for p in args.cost_provenance],
        recovery_owner_wall_seconds=recovery_wall,
        cumulative_case_wall_seconds=total, case_budget_seconds=CASE_BUDGET,
        remaining_case_wall_seconds=CASE_BUDGET-total,
        original_electron_native_wall_seconds=electrons['native_wall_seconds'],
        previous_failed_gamma_native_wall_seconds=old_native['wall_seconds'],
        new_gamma_native_wall_seconds=native_wall,
        cumulative_native_wall_seconds=(electrons['native_wall_seconds'] +
            old_native['wall_seconds'] + native_wall),
        accounting='Owner wall is additive across executions; Native is an overlapping subset. Prior 300s Gamma failure remains paid. No previous component/mean/control costs are erased; their existing receipts remain separate.')


def recover(args):
    started = perf_counter()
    out = args.output.resolve()
    if not re.fullmatch(r'recovery-[A-Za-z0-9][A-Za-z0-9_.-]*', args.attempt_id):
        raise ValueError('attempt-id must be a fresh recovery-NAME path component')
    positive(args.timeout, 'timeout')
    if args.timeout > 600:
        raise ValueError('Gamma Native timeout may not exceed 600 seconds')
    positive(args.additional_paid_wall_seconds, 'additional paid wall', zero=True)
    if args.prior_paid_wall_seconds != PAID_WALL:
        raise ValueError('the actual prior 638.905548795592s case cost must not be reset')
    if args.prior_paid_wall_seconds + args.additional_paid_wall_seconds + args.timeout >= CASE_BUDGET:
        raise ValueError('Native cap leaves no owner/readout wall inside the cumulative case budget')
    prereg = read_json(HERE.parent/'0085'/'preregistration.json')
    case = read_json(HERE/'preregistration.json')
    if out != Path(case['output']).resolve():
        raise ValueError('only the original existing A84 output root is accepted')
    for item in prereg['grounding_actual_sources']:
        authenticate(item)
    original_binding = binding(out/'invocation.json')
    original = read_json(out/'invocation.json')
    if (original['schema'] != 'p254-history-continuation-invocation-v1'
            or original['sourcekind'] != SOURCEKIND
            or original['reference_source_kind'] != REFERENCE
            or original['pi_convention'] != PI_CONVENTION
            or original['projected_control'] is not None):
        raise ValueError('literal original A63 invocation conventions required; no control schema substitution')
    for item in original['source_bindings']:
        authenticate(item)
    invocation = authenticated_json(prereg['actual_original_Gamma_invocation'])
    packed = authenticate(invocation['input'])
    if invocation['input'] != prereg['actual_packed_Gamma_input']:
        raise ValueError('Gamma input does not match preregistered original bytes')
    if invocation['sources'] != original['source_bindings']:
        raise ValueError('Gamma and original continuation source bindings differ')
    old_native = authenticated_json(prereg['actual_failed_Gamma_case'])
    if old_native.get('failure') != 'timeout' or old_native.get('exit') is not None:
        raise ValueError('expected original 300s failed Gamma Native receipt')
    manager = read_json(args.prior_manager_receipt)
    if manager.get('exit') != 1 or manager['wall_seconds'] != args.prior_paid_wall_seconds:
        raise ValueError('actual original failed manager receipt/cost differs')
    header, times, order, beta, layout = packed_header(packed, invocation, original, args)
    electrons = validate_electrons(out, times, order, beta,
        prereg['actual_electron_bank_complete'])
    old_stage = packed.parent
    protected = [binding(p) for p in sorted(old_stage.iterdir()) if p.is_file()]
    if any((out/name).exists() for name in ('readout.json', 'result.json')):
        raise ValueError('final readout/result already exists; no overwrite')
    if (out/'provider.json').exists() and not args.finalize_only:
        raise ValueError('bank provider already exists; only explicit finalize-only is allowed')
    runtime = allocation(args)
    binary = binding(args.binary)
    if binary['sha256'] == invocation['binary']['sha256']:
        raise ValueError('unchanged old Native retry is forbidden')
    evidence = [binding(p) for p in args.source_provenance]
    oracle = binding(args.oracle_record)
    reference = binding(args.reference_provenance)
    for path in args.cost_provenance:
        binding(path)
    sources = [*original['source_bindings'], binding(__file__), binary, oracle,
        reference, *evidence, runtime['provenance'],
        binding(HERE.parent/'0085'/'preregistration.json')]
    stage = out/'gamma'/'.staging'/args.attempt_id
    recovery_invocation = dict(schema=RECOVERY_SCHEMA, original_invocation=original_binding,
        original_gamma_invocation=prereg['actual_original_Gamma_invocation'],
        previous_failed_gamma=prereg['actual_failed_Gamma_case'], input=invocation['input'],
        binary=binary, consumer=binding(__file__), source_bindings=sources,
        finite_sensitive_oracle=oracle, reference_provenance=reference,
        allocation=runtime, layout=layout, timeout_seconds=args.timeout,
        caller_limits=dict(max_input_bytes=args.max_input_bytes,
            max_result_bytes=args.max_result_bytes, max_bank_bytes=args.max_bank_bytes),
        prior_paid_wall_seconds=args.prior_paid_wall_seconds,
        additional_paid_wall_seconds=args.additional_paid_wall_seconds,
        additional_cost_provenance=[binding(p) for p in args.cost_provenance],
        recovery='One new equivalent Gamma-only Native invocation; no in-flight Native resume. Original packed input and failed staging retained; all128 electrons read-only reused.')
    if args.finalize_only:
        prior = read_json(stage/'invocation.json')
        # A new finalization worker may have its own allocation record, but may
        # not relabel the executed binary, input, oracle or physical sources.
        for key in ('schema', 'original_invocation', 'original_gamma_invocation',
                    'input', 'binary', 'consumer', 'finite_sensitive_oracle',
                    'reference_provenance', 'layout'):
            if prior[key] != recovery_invocation[key]:
                raise ValueError('finalize-only execution provenance differs: '+key)
        sources = prior['source_bindings']
    else:
        stage.mkdir(exist_ok=False)
        atomic_json(stage/'invocation.json', recovery_invocation)
    execution_path = stage/('execution-%d.json'%time_ns())
    try:
        with ExitStack() as stack:
            mean = stack.enter_context(SavedMean(original['input_identity']['source_artifact']['path']))
            previous = stack.enter_context(PreviousHistory(original['input_identity']['previous_artifact']['path']))
            kernel_item = original['source_bindings'][0]
            catalog = authenticated_json(kernel_item)
            identities = dict(kernel_input=catalog['input_dataset'],
                numerical_reference=previous.reference_identity,
                source_artifact=mean.binding, previous_artifact=previous.binding)
            if (identities != original['input_identity']
                    or mean.reference_identity != previous.reference_identity
                    or catalog['input_dataset'] != previous.reference_identity
                    or not np.array_equal(mean.times, times)
                    or not np.array_equal(previous.times, times)
                    or not np.array_equal(1j*mean.rho0, previous.aux['source_electron_less0'])
                    or not np.array_equal(-1j*mean.w0[np.ix_(POSITIONS, POSITIONS)].T,
                        previous.aux['source_gamma_less'][0, 0])):
                raise ValueError('original reference/source preparation identity changed')
            projected = projected_control_readout(previous, None)
            readout_control = projected.get('control')
            with FinalBank(out/'gamma', resume=True, max_bank_bytes=args.max_bank_bytes) as bank:
                validate_bank(dict(ret=bank.ret, less=bank.less, times=bank.times,
                    statistics=bank.meta['statistics'], order=bank.meta['order'],
                    beta=bank.meta['beta']), times, 114, 1, 1, order, beta)
                if bank.meta['layout'] != layout or bank.meta['sourcekind'] != SOURCEKIND:
                    raise ValueError('Gamma bank layout/sourcekind differs')
                if args.finalize_only:
                    if set(bank.completed) != {0}:
                        raise ValueError('finalize-only requires the previously committed Gamma slot; never runs Native')
                    provenance = bank.completed[0]['provenance']
                    if provenance['recovery_invocation'] != binding(stage/'invocation.json'):
                        raise ValueError('committed Gamma belongs to another recovery attempt')
                    for key in ('input', 'binary', 'native', 'output'):
                        authenticate(provenance[key])
                    native = authenticated_json(provenance['native'])
                    if native.get('exit') != 0 or native.get('failure'):
                        raise ValueError('committed Gamma Native was not successful')
                else:
                    if bank.completed or (out/'gamma'/'complete.json').exists():
                        raise ValueError('Gamma already committed/finalized; no overwrite or Native replay')
                    output = stage/'result.bin'
                    argv = [str(args.binary.resolve()), str(packed), str(output.resolve()),
                        str(args.max_input_bytes)]
                    remaining = CASE_BUDGET-args.prior_paid_wall_seconds-args.additional_paid_wall_seconds-(perf_counter()-started)
                    if remaining <= args.timeout:
                        raise ValueError('preflight consumed the reserved Native/readout budget')
                    run_recorded(['/usr/bin/time', '-v', '-o', str(stage/'native.time'),
                        *argv], stage, 'native', args.timeout)
                    provenance = dict(invocation, binary=binary, consumer=binding(__file__),
                        sources=sources, output=binding(output), native=binding(stage/'native.json'),
                        recovery_invocation=binding(stage/'invocation.json'),
                        previous_failed_gamma=prereg['actual_failed_Gamma_case'],
                        staging_declared=[str(packed), str(output)],
                        recovery=recovery_invocation['recovery'])
                    with ResultView(output, header, times,
                            max_result_bytes=args.max_result_bytes) as result:
                        result.copy_into(bank.ret[0], bank.less[0])
                    bank.commit(0, provenance)
                if not (out/'gamma'/'complete.json').exists():
                    bank.finalize()
            if previous.provider is not None:
                old = previous.provider.spec
                auxiliary = {name: old['auxiliary'][name] for name in SOURCE_AUXILIARY}
                controls = old.get('represented_control')
                if controls is None and previous.control is not None:
                    controls = {name: old['auxiliary'][name] for name in CONTROL_AUXILIARY}
            else:
                auxdir = stage/'auxiliary'
                if args.finalize_only and (stage/'auxiliary.json').exists():
                    saved = read_json(stage/'auxiliary.json')
                    auxiliary, controls = saved['auxiliary'], saved['controls']
                    for item in [*auxiliary.values(), *(controls or {}).values()]:
                        authenticate(item)
                else:
                    auxdir.mkdir(exist_ok=False)
                    auxiliary = {name: save_array(auxdir/(name+'.npy'), previous.aux[name])
                        for name in SOURCE_AUXILIARY}
                    controls = None if previous.control is None else {
                        name: save_array(auxdir/(name+'.npy'), previous.control[name])
                        for name in CONTROL_AUXILIARY}
                    atomic_json(stage/'auxiliary.json', dict(auxiliary=auxiliary, controls=controls))
            provider = dict(schema=SCHEMA, status='complete', times_inverse_ry=mean.times.tolist(),
                electron=dict(complete=binding(out/'electron'/'complete.json')),
                gamma=dict(complete=binding(out/'gamma'/'complete.json')), auxiliary=auxiliary,
                represented_control=controls, projected_control=readout_control, input_identity=identities,
                sourcekind=SOURCEKIND, integration_order=order, beta=beta,
                reference_source_kind=REFERENCE, pi_convention=PI_CONVENTION,
                input_kernels=kernel_item, input_mean=mean.binding,
                input_previous=previous.binding, source_bindings=sources,
                prescribed_nonGamma_bath='literal original bare provider; no dressed other-q update',
                interpretation='One complete supplied-history Picard iterate; not a converged fixed point or a temporal/error/conservation certificate.')
            # Exact A63 sequence: provider-complete means complete banks, not a
            # completed run. Only result.json marks a passed SavedHistory readout.
            if (out/'provider.json').exists():
                if read_json(out/'provider.json') != provider:
                    raise ValueError('existing bank provider differs; no overwrite')
            else:
                atomic_json(out/'provider.json', provider)
            with SavedHistory(out/'provider.json') as history:
                readout = consume(history, previous, projected_control=None)
            # Recheck protected original artifacts, including the entire failed
            # staging inventory. Electron arrays have only ever been mapped read-only.
            for item in [*protected, electrons['complete'], original_binding]:
                authenticate(item)
            if readout['all128_canonical_fermions_consumed'] != 128 or readout['all_Gamma_position_coordinates_consumed'] != 114:
                raise ValueError('actual full original readout incomplete')
            costs = cost_snapshot(args, stage, started, electrons, old_native)
            if costs['remaining_case_wall_seconds'] <= 0:
                raise ValueError('cumulative case budget exceeded; no complete publication')
            (out/'readout.json').write_text(json.dumps(readout, default=json_value, indent=2)+'\n')
            result = dict(schema='p254-common-history-continuation-result-v1', status='complete',
                provider=binding(out/'provider.json'), readout=binding(out/'readout.json'),
                invocation=original_binding, elapsed_seconds=costs['cumulative_case_wall_seconds'],
                owner_maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                native_invocations=130, successful_native_invocations=129,
                reused_electron_sectors=128, new_gamma_sectors=1,
                previous_failed_gamma_invocations=1, recovery_native_invocations=1,
                times_inverse_ry=mean.times.tolist(),
                physical_scope='All128 active12-band sectors, full114 Gamma; original prescribed nonGamma bath.',
                delta_from_previous=readout['delta_from_previous'],
                projected_comparison=readout['projected_comparison'],
                recovery=dict(invocation=binding(stage/'invocation.json'),
                    gamma_receipt=binding(out/'gamma'/'receipts'/'slot-000.json'),
                    previous_failed_gamma=prereg['actual_failed_Gamma_case'],
                    electron_reuse=electrons, costs=costs,
                    finalization_allocation=runtime, source_bindings=sources),
                limits=['Raw iterate differences are not truncation error, convergence or conservation closure.',
                    'No projection/refill, inferred independent companions, physical calibration, loss, heat or nuclear rate.',
                    'No full31 claim unless the actual complete saved31-node mean and common kernels were supplied.'])
            atomic_json(out/'result.json', result)
        atomic_json(execution_path, dict(status='complete', elapsed_seconds=perf_counter()-started,
            result=binding(out/'result.json')))
        return result
    except BaseException as exc:
        atomic_json(execution_path, dict(status='failed', elapsed_seconds=perf_counter()-started,
            error=repr(exc), costs=cost_snapshot(args, stage, started, electrons, old_native),
            preserved='Original failed staging, packed input, electron bank, Gamma bank and new result/receipts retained; no cleanup or automatic Native retry.'))
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='Original existing A84 output root.')
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--attempt-id', required=True, help='Fresh recovery-NAME; never slot-000.')
    parser.add_argument('--timeout', type=float, required=True, help='New Gamma Native cap, at most600s.')
    parser.add_argument('--max-input-bytes', type=int, required=True)
    parser.add_argument('--max-result-bytes', type=int, required=True)
    parser.add_argument('--max-bank-bytes', type=int, required=True)
    parser.add_argument('--cpu', type=int, required=True)
    parser.add_argument('--memory-bytes', type=int, required=True)
    parser.add_argument('--swap-bytes', type=int, required=True)
    parser.add_argument('--threads', type=int, required=True)
    parser.add_argument('--allocation-record', type=Path, required=True)
    parser.add_argument('--source-provenance', type=Path, action='append', required=True,
        help='Repeat for new native sources/build record and any additional Root source bindings.')
    parser.add_argument('--oracle-record', type=Path, required=True,
        help='Root actual passed finite-sensitive oracle report, bound without inventing its schema.')
    parser.add_argument('--reference-provenance', type=Path, required=True)
    parser.add_argument('--prior-manager-receipt', type=Path,
        default=HERE/'receipts-001'/'manager.json')
    parser.add_argument('--prior-paid-wall-seconds', type=float, default=PAID_WALL)
    parser.add_argument('--additional-paid-wall-seconds', type=float, required=True,
        help='Root actual additional compilation/oracle/other case wall not already in the638.9055s.')
    parser.add_argument('--cost-provenance', type=Path, action='append', required=True)
    parser.add_argument('--finalize-only', action='store_true',
        help='Never run Native; validate this attempt already committed Gamma, then finish readout/publication.')
    print(json.dumps(recover(parser.parse_args())), flush=True)


if __name__ == '__main__':
    main()
