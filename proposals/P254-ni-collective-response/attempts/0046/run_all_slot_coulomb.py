"""Future full128/all64q/all16-pair bare second-Born CLI; Root owns run grants.

This module can consume an actual SavedTransitions bundle, or an imported
reference provider factory(coefficients,times) exposing electron(a,b) and
nonempty original-byte bindings. No U is inferred from rho. --bounds-only
loads no numerical material. The original feature-pair route was measured too
slow; this shared144x144 backend has finite-oracle/synthetic-kernel evidence,
NOT an actual complete material throughput measurement or completion claim.
"""
import argparse
import importlib
import inspect
import json
from pathlib import Path
import sys
from time import perf_counter
import resource
import numpy as np

HERE = Path(__file__).resolve().parent
A44 = HERE.parent/'0044'
sys.path.insert(0,str(A44))
from coulomb_material_inputs import load_coulomb_inputs, DEFAULT_MANIFEST, FrozenReference
from bare_coulomb_factors import BareCoulombFactors
from material_electron_ion_inputs import SavedTransitions
from run_material_bare_coulomb import MATCHING as ORIGINAL_MATCHING
from all_slot_coulomb_contractions import source_bounds as original_source_bounds
from all_slot_coulomb_shared import evaluate_time_grid,bounds as shared_bounds
from all_slot_coulomb_checkpoint import (AllSlotPrefix,binding,authenticate,atomic_json,
                                          publish_result,consume_result)


def source_bounds(*,batch=4):
    result = original_source_bounds(batch=batch)
    result.update(algorithm='shared(q,c,timepair) Greens/static144x144 superoperator',
        conservative_algorithm_ram_bytes_upper=3584*1024**2,
        exact_shared_numeric_bounds=shared_bounds(external_batch=batch,time_batch=8),
        original_material_route='actual q0 ordered calls131-133s; no completedq',
        shared_material_throughput='UNMEASURED; synthetic full-width forecast exceeds current grant',
        material_run_authorized='Root WaveII grant conditional on forecast; currently paused')
    return result


def specification(c,reference,times,batch,extra_sources=()):
    shape = [c.ns,len(times),len(times),c.m,c.m]
    sources = [binding(p) for p in sorted(HERE.glob('*.py'))]
    sources += [binding(A44/name) for name in ('bare_coulomb_factors.py','second_born_functional.py',
        'material_second_born.py','coulomb_material_inputs.py','run_material_bare_coulomb.py')]
    sources += [binding(HERE.parent/'0041'/'material_electron_ion_inputs.py'),*extra_sources]
    history = list(reference.bindings)
    inputs = [c.bindings['manifest'],c.bindings['original_boundary'],
        c.bindings['background_receipt'],c.bindings['background_archive'],
        *c.bindings['authenticated_archives'],*c.bindings['producer_implementation'],*history]
    for item in sources+inputs:
        authenticate(item)
    matching = dict(ORIGINAL_MATCHING,reference=reference.label,
        limitations=['not selfconsistent KB dynamics','no passive/complement collision loops',
                     'not full-G or all-electron','no screened-W','no rate/passivity/whole-E claim'])
    return dict(schema='p254-a46-full-slot-specification-v1',shape=shape,nq=c.nq,
        external_slots=list(range(c.ns)),times_inverse_ry=list(map(float,times)),
        external_batches=[list(range(i,min(i+batch,c.ns))) for i in range(0,c.ns,batch)],
        matching=matching,reference_history=history,input_bindings=c.bindings,
        authenticated_sources_and_inputs=sources+inputs,bounds=source_bounds(batch=batch))


def run_reference(c,reference,times,output,*,external_batch=4,extra_sources=()):
    """Numerical provider entrypoint; independent ordered pairs are all gathered."""
    spec = specification(c,reference,times,external_batch,extra_sources)
    complete = consume_result(output,spec)
    if complete is not None:
        return dict(status='complete-reused',archive=complete['record']['archive'],
            decoded_norms=complete['decoded_norms'],counts=complete['record']['counts'])
    electrons = {(a,b):reference.electron(a,b) for a in range(len(times)) for b in range(len(times))}
    if any(len(legs) != 2 or any(x.shape != (c.ns,c.m,c.m) or not np.isfinite(x).all()
        for x in legs) for legs in electrons.values()):
        raise ValueError('provider ordered electron legs differ')
    checkpoint = AllSlotPrefix(Path(output)/'numeric-prefix',spec,spec['shape'],c.nq)
    started = perf_counter()
    arrays,counts = evaluate_time_grid(BareCoulombFactors(c),electrons,times,checkpoint,
        external_batch=external_batch,after_commit=lambda q:print(json.dumps(dict(
            completed_q=q,external_slots=c.ns,elapsed_seconds=perf_counter()-started)),flush=True))
    result = publish_result(output,spec,arrays,counts,checkpoint.recoveries)
    return dict(status='complete',archive=result['record']['archive'],
        decoded_norms=result['decoded_norms'],counts=counts,seconds=perf_counter()-started,
        peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',default=DEFAULT_MANIFEST)
    p.add_argument('--output',type=Path)
    p.add_argument('--resume',action='store_true')
    p.add_argument('--bounds-only',action='store_true')
    group = p.add_mutually_exclusive_group()
    group.add_argument('--transitions',type=Path,help='actual p254-electron-ion-transitions-v1 bundle')
    group.add_argument('--provider',help='module:factory returning actual source-bound reference')
    group.add_argument('--frozen',action='store_true',help='explicit conditional frozen comparison only')
    p.add_argument('--external-batch',type=int,choices=range(1,5),default=4)
    p.add_argument('--density-cap-mib',type=int,default=900)
    args = p.parse_args()
    if args.bounds_only:
        print(json.dumps(source_bounds(batch=args.external_batch),indent=2)); return
    if args.output is None or not (args.transitions or args.provider or args.frozen):
        p.error('--output and an explicit actual reference choice are required')
    if args.output.exists() and not args.resume:
        raise FileExistsError('output exists; exact-source --resume required')
    args.output.mkdir(parents=True,exist_ok=True)
    atomic_json(args.output/'invocation.json',dict(argv=sys.argv,start_unix=__import__('time').time(),
        note='operator/Root must separately grant actual material execution'))
    c = load_coulomb_inputs(args.manifest,density_cap_bytes=args.density_cap_mib*1024**2)
    times = np.array([0.,1.,10.,30.])
    extra = []
    if args.transitions:
        reference = SavedTransitions(c,times,args.transitions,c.bindings['manifest'])
    elif args.provider:
        module,name = args.provider.split(':',1)
        factory = getattr(importlib.import_module(module),name)
        extra = [binding(inspect.getsourcefile(factory))]
        reference = factory(c,times)
        if not reference.bindings:
            raise ValueError('actual reference provider must bind its physical history')
    else:
        reference = FrozenReference(c,times)
    try:
        print(json.dumps(run_reference(c,reference,times,args.output,
            external_batch=args.external_batch,extra_sources=extra)),flush=True)
    finally:
        if hasattr(reference,'close'):
            reference.close()


if __name__ == '__main__':
    main()
