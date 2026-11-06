"""Actual A27 Bloch admission/probe: --admission JSON --output DIR [--context-only].

Uses compound A25 runtime authority plus bound A27 continuation/implementation/
passed-oracle settings. Inputs: xml,density,upfs{Ni,H},wavefunctions{1:{1..10},
2:{1..10}},reference_terminal (original accepted A25). Exact settings are in
material_bloch_admission.validate_recipe. Resume additionally binds each
completed source-N-K or coverage-probe receipt; unchanged recipe, reduced clock.
No old producer, PBE backend, Native action or finite-q placeholder is invoked.
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
for directory in (HERE.parent/'0024',HERE.parent/'0025',HERE.parents[3]/'src'):
    sys.path.insert(0,str(directory))


def main():
    run = publisher = None
    summary = dict(status='failed',generation='0027',run_id=uuid.uuid4().hex,outputs=[],
        first_useful_block_seconds=None,first_useful_target_seconds=30,
        scope='all20 original sources/all128 lazy spin-grid slots; selected nonGamma raw q0 block only',
        state=None,D=None,temperature=None,finite_q=None,
        recovery_boundary='complete original-source archives and final complete coverage/probe only',
        source_boundary_hard_seconds=30,probe_boundary_hard_seconds=90,saved_readback_hard_seconds=45,
        time_accounting='All setup/admission/readback/publication charged; terminal is pre-publication snapshot; manager post-exit success separate')
    code = 1
    try:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('--admission',required=True,type=Path)
        parser.add_argument('--output',required=True,type=Path)
        parser.add_argument('--context-only',action='store_true')
        args = parser.parse_args()
        target = args.output.resolve()
        target.mkdir(parents=True,exist_ok=True)
        def interrupted(signum, frame):
            raise InterruptedError(f'signal {signum}; only complete immutable boundaries recover')
        for signum in (signal.SIGTERM,signal.SIGINT,signal.SIGALRM):
            signal.signal(signum,interrupted)
        from material_run_contract import runtime_context,implementation,admit,binding
        from material_reference_output import OutputRun
        from material_publication import Publication
        from substrate_framework.material_bloch_admission import validate_grant
        from substrate_framework.material_bloch_production import produce
        from substrate_framework.material_bloch_output import preload_parsers
        preload_parsers()
        paths = (Path(__file__),)
        context,sources = runtime_context(),implementation(paths)
        summary.update(context=context,implementation_sources=sources)
        grant,recipe,identity = admit(args.admission,'material-bloch-grid',sources,context,args.context_only)
        if args.context_only:
            publisher = Publication(target,summary,BOOT,grant,context_only=True)
            publisher.arm()
        else:
            validate_grant(grant)
            run = OutputRun(target,summary,BOOT)
            run.configure(grant,identity)
            publisher = run.publisher
            summary.update(recipe=recipe,recipe_sha256=identity,admission=binding(args.admission),
                           remaining_seconds_granted=grant['remaining_seconds'],
                           first_block_seconds_granted=grant['first_block_seconds'],
                           mode='explicit-resume' if grant['resume'] else 'new')
        context_path = target/('context-'+summary['run_id']+'.json')
        publisher.write(context_path,dict(generation='0027',run_id=summary['run_id'],
                                         context=context,implementation_sources=sources))
        publisher.event('material_bloch_context',context=binding(context_path))
        if args.context_only:
            summary['status'] = 'context-only'
        else:
            produce(run,recipe,paths)
        code = 0
    except BaseException as error:
        summary.update(status='failed',failure=dict(type=type(error).__name__,message=str(error),
                                                   traceback=traceback.format_exc()))
        if isinstance(error,SystemExit):
            code = int(error.code or 0)
    finally:
        if publisher is not None:
            code = publisher.finish(code)
        elif code:
            print(json.dumps(dict(event='publication_not_admitted',run_id=summary['run_id'],
                                  exit_code=code,error=str(summary.get('failure'))[:240])),flush=True)
    return code


if __name__ == '__main__':
    sys.exit(main())
