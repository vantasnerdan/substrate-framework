#!/usr/bin/env python3
"""Actual37-ion finite ionic scalar/111-force pilot; Root grant REQUIRED.

CLI: --admission JSON --output DIR [--context-only]. Existing A24/A25
OutputRun/Publication/compound grant controls remain authoritative. Inputs are
xml/upfs{Ni,H}/native_stdout original bindings. Explicit settings are
material_ionic_recipe's five A27 authority bindings, numerical Ewald contract,
Native evaluation_index=0 and original coordinate identity. No library_binding.

One immutable whole-ionic-baseline.npz/.json boundary retains scalar energy,
real/reciprocal/self/background terms, full force/gradient, geometry/valences,
original Native ion-ion decimal/stage identity and all signed differences.
Resume requires an original bound receipt at resume_receipts[whole-ionic-baseline]
and the identical recipe plus a fresh remaining-time grant. Unfinished sums
are lost, never resumed. No automatic retry, cutoff tightening or force repair.
"""
import time
BOOT = time.monotonic()
import argparse
import json
from pathlib import Path
import signal
import sys
import traceback
import uuid

ROOT = Path(__file__).resolve().parents[4]
ATTEMPTS = Path(__file__).resolve().parents[1]
for path in (ROOT/'src', ATTEMPTS/'0024', ATTEMPTS/'0025'):
    sys.path.insert(0, str(path))


def produce(run, grant, summary, recipe):
    from material_run_contract import bound_json, implementation, verify_consumers
    from substrate_framework.material_ionic_recipe import (
        KEY, COORDINATES, block_metadata, load_ionic_inputs, sum_contract, validate_grant)
    from substrate_framework.material_ionic_arrays import baseline_arrays, validate_arrays
    from substrate_framework.material_ionic_output import load_saved_ionic_block
    from substrate_framework.qe_ionic_ewald import periodic_ionic_ewald
    validate_grant(grant)
    if not grant['resume'] and any((run.target/(KEY+suffix)).exists() for suffix in ('.npz', '.json')):
        raise ValueError('fresh execution refuses any existing immutable ionic boundary')
    run.check()
    verify_consumers(grant, implementation((Path(__file__),)))
    start = time.monotonic()
    data = load_ionic_inputs(grant['inputs'], grant['settings'])
    summary['timing_seconds']['original_xml_upf_native_decode'] = time.monotonic()-start
    run.check()
    start = time.monotonic()
    if grant['resume'] and KEY in grant['resume_receipts']:
        expected = bound_json(grant['resume_receipts'][KEY])
        saved = load_saved_ionic_block(run.target, run.identity, expected, recipe, data=data)
        run.record(KEY, dict(saved.receipt), 'completed-cache-reuse')
        summary['completion_mode'] = 'completed-whole-baseline-reuse; no Ewald evaluation'
    else:
        term = periodic_ionic_ewald(
            data.reference.cell_bohr, data.reference.positions_bohr.T, data.valence_charges,
            contract=sum_contract(grant['settings']), coordinate_identity=COORDINATES,
            source_authority='Original bound saved XML geometry and ordinary NCUPF VALENCE; no nuclear-charge substitution, fitted force or support',
            sources=data.sources)
        run.check()
        arrays = baseline_arrays(term, data)
        metadata = block_metadata(recipe, data)
        metadata['comparison'] = validate_arrays(arrays, data)
        receipt = run.commit(KEY, arrays, metadata)
        # Consume the real committed archive, not the still-resident candidate.
        saved = load_saved_ionic_block(run.target, run.identity, receipt, recipe, data=data)
        summary['completion_mode'] = 'new-whole-baseline; committed archive independently restored'
    summary['timing_seconds']['whole_baseline_or_restore_including_archive_readback'] = time.monotonic()-start
    summary['ionic_result'] = dict(energy_ry=saved.term.energy_ry,
                                  energy_terms_ry=dict(saved.term.energy_terms_ry),
                                  sum_metadata=dict(saved.term.sum_metadata),
                                  comparison=dict(saved.comparison),
                                  saved_readback='authenticated all111 coordinates and typed full archive; Ewald not rerun')
    summary['input_sources'] = list(data.sources)+[dict(data.native_ledger.source)]
    summary['recovery_receipt'] = dict(saved.receipt)
    run.check()
    sources = implementation((Path(__file__),))
    verify_consumers(grant, sources)
    summary['final_implementation_sources'] = sources
    summary['source_coverage_after_lazy_consumers'] = 'actual producer and archive-reader inventory authenticated'
    run.check()


def main():
    run = None
    code = 1
    summary = dict(status='failed', generation='0027', run_id=uuid.uuid4().hex,
                   outputs=[], first_useful_block_seconds=None, timing_seconds={},
                   scope='full37-ion finite ionic scalar/111-force only; no prepared total-force or equilibrium claim',
                   recovery_scope='completed whole baseline only; all unfinished sums lost; setup charged on reuse')
    try:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('--admission', required=True, type=Path)
        parser.add_argument('--output', required=True, type=Path)
        parser.add_argument('--context-only', action='store_true')
        args = parser.parse_args()
        target = args.output.resolve()
        target.mkdir(parents=True, exist_ok=True)
        from material_run_contract import admit, implementation, runtime_context, binding
        from material_reference_output import OutputRun
        from material_publication import Publication
        # Preload the real producer, XML/UPF/stdout and saved-output consumers
        # before source inventory. Context-only does not decode material.
        from substrate_framework import material_ionic_arrays, material_ionic_output
        from substrate_framework.material_ionic_recipe import KIND, preload_parsers
        preload_parsers()
        run = OutputRun(target, summary, BOOT)

        def interrupted(signum, frame):
            raise InterruptedError(f'signal {signum}; only completed whole ionic baseline recovers')

        for signum in (signal.SIGTERM, signal.SIGINT, signal.SIGALRM):
            signal.signal(signum, interrupted)
        summary['context'] = runtime_context()
        summary['implementation_sources'] = implementation((Path(__file__),))
        context = dict(generation='0027', run_id=summary['run_id'], context=summary['context'],
                       implementation_sources=summary['implementation_sources'])
        grant, recipe, identity = admit(args.admission, KIND, summary['implementation_sources'],
                                       summary['context'], args.context_only)
        run.publisher = Publication(target, summary, BOOT, grant, context_only=args.context_only)
        run.publisher.arm()
        if not args.context_only:
            summary.update(recipe=recipe, recipe_sha256=identity, admission=binding(args.admission),
                           remaining_seconds_granted=grant['remaining_seconds'], resume=grant['resume'])
            run.configure(grant, identity)
        context_path = target/('context-'+summary['run_id']+'.json')
        run.publisher.write(context_path, context)
        if args.context_only:
            print(json.dumps(dict(event='material_ionic_context', **context), allow_nan=False), flush=True)
            summary['status'] = 'context-only'
        else:
            run.event('material_ionic_context', path=str(context_path), context=summary['context'])
            produce(run, grant, summary, recipe)
            summary['status'] = 'complete'
        code = 0
    except BaseException as error:
        summary['failure'] = dict(type=type(error).__name__, message=str(error), traceback=traceback.format_exc())
        if isinstance(error, SystemExit) and error.code == 0:
            code = 0
    finally:
        if run is not None:
            code = run.finish(code)
        elif code:
            print(json.dumps(dict(event='pre_output_failure', run_id=summary['run_id'], exit_code=code,
                                  error=str(summary.get('failure'))[:240])), flush=True)
    return code


if __name__ == '__main__':
    raise SystemExit(main())
