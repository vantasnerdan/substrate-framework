#!/usr/bin/env python3
"""Bound conditional neutral preparation and one ORIGINAL Gamma raw-D consumer.

CLI: --admission JSON --output DIR [--context-only]. Shared registration and
source_access remain the compound A25 authority. Kind=material-neutral-preparation,
cap120s cumulative including production/readback; Root supplies reduced grants.
inputs: bloch_terminal (ACTUAL completed authenticated A27), factor_terminal
(original admitted A25). settings: continuation_registration, implementation_registration,
preparation_oracle bindings; preparation_recipe='neutral-free-M-saved-label-exact-block-empty-complement-v1';
selected_spin=1. No library_binding. Fresh resume_receipts={}; recovery maps
neutral-populations and/or spin1-gamma-preparation to exact prior receipt bindings,
with identical recipe/resume_identity and separately reduced remaining_seconds.
No material input is read in --context-only. No actual execution is authorized
by this source file. First useful output is the complete neutral ledger; a
complete run also requires actual all-row/all-column raw-D and probe output.
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

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
for directory in (ROOT / 'src', HERE.parent / '0025', HERE.parent / '0024'):
    sys.path.insert(0, str(directory))


def preload():
    """Load every lazy scientific/readback/parser consumer before source inventory."""
    from substrate_framework import (material_preparation_inputs, material_preparation_gamma,
        material_preparation_output, material_preparation_consumer, material_preparation_restore)
    material_preparation_inputs.preload_parsers()
    return (Path(__file__), *(Path(module.__file__) for module in (
        material_preparation_inputs, material_preparation_gamma, material_preparation_output,
        material_preparation_consumer, material_preparation_restore)))


def main():
    run, code = None, 1
    summary = dict(status='failed', generation='0027', run_id=uuid.uuid4().hex,
        outputs=[], timing_seconds={}, first_useful_block_seconds=None,
        scope='complete mean-neutral population ledger plus separately selected spin1 Gamma raw preparation',
        ensemble='conditional independent spin/full-k quasifree fixed-mean N; free M/common multiplier',
        stationarity='not imposed or established; no material propagation', physical_temperature=None,
        ionic=None, cross=None, held_state=None, total=None,
        recovery_scope='completed neutral ledger / whole selected Gamma result only; unfinished phase lost; no automatic retry',
        time_accounting='boot-to-publication pre-exit snapshot; actual post-exit manager outcome/readback charged separately')
    try:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('--admission', type=Path, required=True)
        parser.add_argument('--output', type=Path, required=True)
        parser.add_argument('--context-only', action='store_true')
        args = parser.parse_args()
        def interrupted(signum, frame):
            raise InterruptedError(f'signal {signum}; only completed preparation boundaries reusable')
        for signum in (signal.SIGTERM, signal.SIGINT, signal.SIGALRM):
            signal.signal(signum, interrupted)
        paths = preload()
        from material_run_contract import admit, implementation, runtime_context, binding
        from material_reference_output import OutputRun
        from material_publication import Publication
        from substrate_framework.material_preparation_consumer import produce
        target = args.output.resolve()
        target.mkdir(parents=True, exist_ok=True)
        run = OutputRun(target, summary, BOOT)
        context, sources = runtime_context(), implementation(paths)
        summary.update(context=context, implementation_sources=sources)
        grant, recipe, identity = admit(args.admission, 'material-neutral-preparation', sources,
                                        context, args.context_only)
        run.publisher = Publication(target, summary, BOOT, grant, context_only=args.context_only)
        run.publisher.arm()
        if not args.context_only:
            summary.update(recipe=recipe, recipe_sha256=identity, admission=binding(args.admission),
                           remaining_seconds_granted=grant['remaining_seconds'], resume=grant['resume'])
            run.configure(grant, identity)
        context_path = target / ('context-' + summary['run_id'] + '.json')
        run.publisher.write(context_path, dict(context=context, implementation_sources=sources,
                                               generation='0027', run_id=summary['run_id']))
        run.event('material_preparation_context', context=binding(context_path))
        if args.context_only:
            summary['status'] = 'context-only'
        else:
            produce(run, grant, paths)
            summary['status'] = 'complete'
        code = 0
    except BaseException as error:
        summary['status'] = 'failed'
        summary['failure'] = dict(type=type(error).__name__, message=str(error), traceback=traceback.format_exc())
        if isinstance(error, SystemExit) and error.code == 0:
            code = 0
    finally:
        if run is not None:
            code = run.finish(code)
        elif code:
            print(json.dumps(dict(event='pre_output_failure', run_id=summary['run_id'],
                                  exit_code=code, error=str(summary.get('failure'))[:240])), flush=True)
    return code


if __name__ == '__main__':
    sys.exit(main())
