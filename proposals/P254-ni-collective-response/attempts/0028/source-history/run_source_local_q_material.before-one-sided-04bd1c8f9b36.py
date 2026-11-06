#!/usr/bin/env python3
"""Selected actual material-noise work using source-local discrete pullbacks.

plan is a metadata forecast; probe measures one original full-PW t0 source;
pair alone evaluates the ordered two-time covariance on the saved actual mean
history. None certifies whole-64q reconstruction, accuracy, heat or coherence.
No mode authorizes production allocation or automatically restarts itself.
"""
import argparse
import json
import os
from pathlib import Path
import sys
from time import time
from uuid import uuid4

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0025'), str(HERE.parent/'0024')]
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'

from mixed_material_driver_support import (
    canonical, check_loaded_sources, owned_directory, source_closure)
from source_local_material_owner import (
    SCHEMA, CAP_NAMES, context_bindings, encode, error_record, execution_caps,
    finish_receipt, isolated_output, material_budget, material_owner,
    reserve_terminal_and_ledger, terminal_signals)
from source_local_material_runtime import (
    build_runtime, material_inputs, material_sources, metadata_plan)
from source_local_material_probe import measure_edge
from substrate_framework.joint_field_probes import OBSERVABLES
from substrate_framework.source_local_q_covariance import SourceLocalTwoSidedCovariance
from substrate_framework.streamed_material_state import file_binding
from run_mixed_q_material import watchdog


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument('--mode', choices=('plan', 'probe', 'pair'), required=True)
    result.add_argument('--case', required=True)
    for name in ('initialization', 'history', 'physical-qualification', 'source-local-oracle'):
        result.add_argument('--'+name, type=Path, required=True)
        result.add_argument('--'+name+'-sha256', required=True)
    result.add_argument('--output', type=Path, required=True)
    result.add_argument('--resume', action='store_true',
        help='Reuse this exact owner, original deadline, action ledger and scientific namespace')
    result.add_argument('--q', nargs=3, type=float, required=True)
    result.add_argument('--cpus', nargs='+', type=int, choices=range(6, 11), required=True)
    for name in CAP_NAMES:
        result.add_argument('--'+name.replace('_', '-'),
                            type=float if name == 'max_seconds' else int, required=True)
    result.add_argument('--time', type=float, required=True)
    result.add_argument('--second-time', type=float, required=True)
    result.add_argument('--subdivisions', type=int, required=True)
    result.add_argument('--atom', type=int, default=37)
    result.add_argument('--axis', choices=('x', 'y', 'z'), default='x')
    result.add_argument('--observable', choices=OBSERVABLES, default='total-force')
    result.add_argument('--field-axis', choices=('x', 'y', 'z'),
                        help='Cartesian global field readout axis, independent of the ionic seed axis')
    result.add_argument('--dynamics', choices=('screened-joint', 'bare-electronic'),
                        default='screened-joint')
    result.add_argument('--block-size', type=int, default=8)
    result.add_argument('--frame-admission', nargs=2, action='append', default=[],
                        metavar=('PATH', 'SHA256'), help='Original validated frame admission; repeatable')
    result.add_argument('--edge-index', type=int, default=0,
                        help='Original source edge for the full-PW t0 probe (default: 0)')
    for side in ('left', 'right'):
        result.add_argument('--'+side+'-context', type=Path,
                            help='Owned complete source-local pullback to authenticate on --resume')
        result.add_argument('--'+side+'-context-sha256')
    return result


def scientific_request(args, sources, caps):
    selectors = {name: getattr(args, name) for name in (
        'mode', 'case', 'q', 'time', 'second_time', 'subdivisions', 'atom', 'axis',
        'observable', 'field_axis', 'dynamics', 'block_size', 'edge_index')}
    explicit = {}
    for name in ('initialization', 'history', 'physical_qualification', 'source_local_oracle'):
        explicit[name] = dict(path=str(getattr(args, name).resolve()),
                              sha256=getattr(args, name+'_sha256'))
    explicit['frame_admissions'] = [dict(path=str(Path(path).resolve()), sha256=sha)
                                   for path, sha in args.frame_admission]
    return canonical(encode(dict(schema=SCHEMA+'-request', selectors=selectors,
        explicit_inputs=explicit,
        source_bindings=sources, caps=caps, blas_threads=1,
        scientific_namespace='scientific-work', original_time_origin=0.,
        continuation='owned source-local contexts only; explicit bindings are authenticated restart hints')))


def execute(args, inputs, sources, budget, directory, contexts, result):
    if args.mode == 'plan':
        result['forecast'] = metadata_plan(args, inputs)
        result['scientific_scope'] = 'metadata-only source-local forecast, not a measured covariance'
        return
    runtime = build_runtime(args, inputs, budget=budget,
                            directory=directory/'scientific-work', source_bindings=sources)
    result.update(provenance=runtime.provenance, resource_plan=runtime.resource_plan,
                  time_grid=runtime.times, second_time_grid=runtime.second_times)
    if args.mode == 'probe':
        result['one_source_diagnostic'] = measure_edge(runtime, index=args.edge_index, time=0.0)
        result['scientific_scope'] = 'one original full-PW t0 source measurement; not covariance'
        return
    result['covariance'] = SourceLocalTwoSidedCovariance(
        runtime.preparation, runtime.executor).pair(runtime.probe, runtime.probe,
        runtime.times, runtime.second_times, directory/'scientific-work'/'pair', **contexts)
    result['covariance_executed'] = True
    result['scientific_scope'] = 'actual selected-q ordered two-time covariance on original saved history'


def run_owned(args, directory, invocation, started_unix, result):
    budget = None
    try:
        caps = execution_caps(args)
        sources = source_closure(material_sources()+[
            file_binding(HERE/'run_source_local_q_material.py'),
            file_binding(HERE/'source_local_material_owner.py')])
        request = scientific_request(args, sources, caps)
        owner, binding, identity = material_owner(directory, request,
            resume=args.resume, started_unix=started_unix)
        result.update(owner=binding, material_identity=identity,
                      original_started_unix=owner['started_unix'],
                      source_bindings=sources, explicit_inputs=request['explicit_inputs'])
        budget = material_budget(directory, owner, caps)
        reserve_terminal_and_ledger(budget)
        contexts = context_bindings(args, directory)
        result['explicit_source_local_contexts'] = contexts
        os.sched_setaffinity(0, {args.cpus[0]})
        with watchdog(budget):
            inputs = material_inputs(args)
            execute(args, inputs, sources, budget, directory, contexts, result)
            check_loaded_sources(sources)
        result['status'] = 'complete'
    except BaseException as error:
        result['status'] = 'interrupted' if isinstance(error, (InterruptedError, KeyboardInterrupt)) else 'failed'
        result['error'] = error_record(error)
    finally:
        result['finished_unix'] = time()
        binding = finish_receipt(directory, invocation, result, budget)
    return binding


def main(arguments=None):
    started_unix = time()
    try:
        args = parser().parse_args(arguments)
    except SystemExit as error:
        if error.code:
            print(json.dumps(dict(schema=SCHEMA+'-exit', status='failed',
                error=error_record(error), receipt_location='stderr; no owned output admitted')),
                file=sys.stderr)
        return int(error.code or 0)
    invocation = uuid4().hex
    result = dict(schema=SCHEMA+'-result', status='failed', mode=args.mode,
        invocation=invocation, invocation_started_unix=started_unix, resume=args.resume,
        cpus=args.cpus, blas_threads=1, covariance_executed=False,
        whole_spatial_reconstruction=False, required_q_for_whole_local=64,
        numerical_accuracy_certified=False, heat_claim=False, coherence_claim=False,
        production_authorized=False, automatic_restart=False,
        terminal_accounting='reserved metadata; pre-budget failures are recorded without scientific admission')
    try:
        with terminal_signals():
            directory = isolated_output(args, ROOT)
            with owned_directory(directory, args.resume) as directory:
                binding = run_owned(args, directory, invocation, started_unix, result)
        print(json.dumps(dict(status=result['status'], receipt=binding), allow_nan=False))
    except BaseException as error:
        # A busy/unsafe output or unavailable filesystem must never be mutated to
        # manufacture an on-disk receipt. Retain the catchable exit on stderr.
        result['status'] = 'failed'
        result['receipt_error'] = error_record(error)
        result['receipt_location'] = 'stderr; owned output unavailable or terminal publication failed'
        print(json.dumps(encode(result), allow_nan=False), file=sys.stderr)
    return 0 if result['status'] == 'complete' else 2


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        raise SystemExit(main())
