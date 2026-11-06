#!/usr/bin/env python3
"""Bound actual prepared-source consumer using the established material runtime.

Kind material-reference; A28 settings identify the common-state construction.
--context-only opens no material. Existing immutable block/readback, source
inventory, own-manager context and publication APIs remain unchanged.
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
for directory in (ROOT/'src', HERE.parent/'0025', HERE.parent/'0024'):
    sys.path.insert(0, str(directory))


def main():
    run, code = None, 1
    summary = dict(status='failed', generation='0028-common-reference', run_id=uuid.uuid4().hex,
        outputs=[], first_useful_block_seconds=None,
        scope='actual full-column conditional source preparation and selected unweighted material moments',
        recovery_scope='completed source/slot archives only; unfinished work is not recoverable')
    try:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('--admission', type=Path, required=True)
        parser.add_argument('--output', type=Path, required=True)
        parser.add_argument('--context-only', action='store_true')
        args = parser.parse_args()
        def interrupted(signum, frame):
            raise InterruptedError('signal '+str(signum)+'; committed blocks retained')
        for signum in (signal.SIGTERM, signal.SIGINT, signal.SIGALRM):
            signal.signal(signum, interrupted)
        import material_reference_producer
        from substrate_framework.material_bloch_output import preload_parsers
        from material_run_contract import admit, binding, implementation, runtime_context
        from material_reference_output import OutputRun
        from material_publication import Publication
        preload_parsers()
        paths = (Path(__file__), Path(material_reference_producer.__file__))
        context, sources = runtime_context(), implementation(paths)
        summary.update(context=context, implementation_sources=sources)
        target = args.output.resolve()
        target.mkdir(parents=True, exist_ok=True)
        run = OutputRun(target, summary, BOOT)
        grant, recipe, identity = admit(args.admission, 'material-reference', sources,
                                        context, args.context_only)
        run.publisher = Publication(target, summary, BOOT, grant, context_only=args.context_only)
        run.publisher.arm()
        if not args.context_only:
            summary.update(recipe=recipe, recipe_sha256=identity, admission=binding(args.admission))
            run.configure(grant, identity)
        context_path = target/('context-'+summary['run_id']+'.json')
        run.publisher.write(context_path, dict(context=context, implementation_sources=sources,
                                               generation=summary['generation'], run_id=summary['run_id']))
        run.event('common_reference_context', context=binding(context_path))
        if args.context_only:
            summary['status'] = 'context-only'
        else:
            material_reference_producer.produce(run, recipe, paths)
        code = 0
    except BaseException as error:
        summary['failure'] = dict(type=type(error).__name__, message=str(error), traceback=traceback.format_exc())
    finally:
        if run is not None:
            code = run.finish(code)
        else:
            print(json.dumps(dict(event='pre_output_failure', failure=summary.get('failure'))), flush=True)
        print('[dft-done] A28 common reference' if code == 0 else '[dft-fail] A28 common reference', flush=True)
    return code


if __name__ == '__main__':
    raise SystemExit(main())
