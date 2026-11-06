"""Evaluate correlation-only Ni/H bare-Coulomb 2B on all original time pairs.

External spin/k slots are declared, never a truncation of internal sums. Root
owns execution authorization. --bounds-only consumes no numerical material.
"""
import argparse
import json
from pathlib import Path
from time import perf_counter
import resource
import numpy as np
from coulomb_material_inputs import load_coulomb_inputs, FrozenReference, DEFAULT_MANIFEST, file_binding, HERE
from bare_coulomb_factors import BareCoulombFactors
from material_second_born import evaluate_material_time_grid, source_bounds
from reduced_projection_cache import ProjectionCache, atomic_json
from material_coulomb_result import consume_complete_result

MATCHING = dict(
    model='finite-feature all64q/full12 active NC valence pseudo-space bare Coulomb',
    interaction='H=1/4 sum A_ij;kl cdag_i cdag_j c_l c_k; A=V_ij;kl-V_ij;lk',
    multiplier='8*pi/|G+q|^2 Ry; admitted |G+q|^2/tpiba^2>1e-8; lnolr=False',
    measure='explicit128spin-k, canonical G without weights; each vertex1/N, N=64; D carries1/sqrt(Omega)',
    source='original density members only; NOT saved PBE K, ionic moments, fitted U or screened W',
    reference='literal frozen original h0/rho0 coherent Gaussian; no occupation refill',
    residual='reference-normal-ordered quartic residual defines correlation-only second Born',
    retained_mean='original DFT/Hartree/PBE mean and original conservative K density feedback remain separate and unchanged',
    excluded_mean='NO fresh Hartree/Fock insertion; normal ordering is not an ab-initio XC counterterm',
    ordered='greater=>,>,reverse<; lesser=<,<,reverse>; retarded=theta(t>=s)*(greater-lesser),theta(0)=1',
    limitations=['not moving/selfconsistent KB dynamics','no passive/complement collision loops',
                 'not full-G or all-electron','no screened-W','no rate/passivity/whole-E claim'])


def external_indices(text):
    result = list(range(128)) if text == 'all' else [int(x) for x in text.split(',')]
    if not result or len(set(result)) != len(result) or min(result)<0 or max(result)>=128:
        raise ValueError('declare distinct external spin-k slots0..127')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',default=DEFAULT_MANIFEST)
    p.add_argument('--external',required=True,help='comma-separated explicit spin-k slots, or all')
    p.add_argument('--bounds-only',action='store_true')
    p.add_argument('--output',type=Path)
    p.add_argument('--resume',action='store_true',help='reuse only byte-bound committed all-time transfer-q chunks')
    p.add_argument('--density-cap-mib',type=int,default=900)
    args = p.parse_args()
    slots = external_indices(args.external)
    bounds = source_bounds(len(slots),16)
    if args.bounds_only:
        print(json.dumps(dict(bounds=bounds,matching=MATCHING),indent=2)); return
    if args.output is None:
        p.error('--output required for an actual material evaluation')
    if args.output.exists() and not args.resume:
        raise FileExistsError('output exists; exact-source --resume required')
    args.output.mkdir(parents=True,exist_ok=True)
    started = perf_counter()
    if args.resume:
        complete = consume_complete_result(args.output,args.manifest,slots,MATCHING,bounds)
        if complete is not None:
            print(json.dumps(complete),flush=True)
            return
    c = load_coulomb_inputs(args.manifest,density_cap_bytes=args.density_cap_mib*1024**2)
    factors = BareCoulombFactors(c)
    times = np.array([0.,1.,10.,30.])
    reference = FrozenReference(c,times)
    electrons = {(a,b):reference.electron(a,b) for a in range(4) for b in range(4)}
    receipt = dict(schema='p254-a44-original-bare-coulomb-material-v1',status='running',
        matching=MATCHING,source_bounds=bounds,input_bindings=c.bindings,
        implementation=[file_binding(path) for path in sorted(HERE.glob('*.py'))],
        times_inverse_ry=times.tolist(),external_slots=slots,items=[],
        external_slot_bindings=[dict(slot=i,spin=int(c.slot_spin[i]),k=int(c.slot_k[i])) for i in slots],
        feature_counts=[len(g) for g in c.feature_g],
        head_admitted=[int(x.sum()) for x in factors.admitted],omega_bohr3=c.omega,tpiba_bohr_inverse=c.tpiba)
    specification = {k:receipt[k] for k in ('schema','matching','source_bounds','input_bindings',
                                           'implementation','times_inverse_ry','external_slots')}
    cache = ProjectionCache(args.output/'transfer-q-chunks',specification)
    record = args.output/'result.json'
    atomic_json(record,receipt)
    for external in slots:
        unit_start = perf_counter()
        arrays,count = evaluate_material_time_grid(factors,external,electrons,times,cache)
        path = args.output/f'slot-{external:03d}-all-original-pairs.npz'
        np.savez(path,**arrays)
        with np.load(path,allow_pickle=False) as saved:
            norms = {key:[[float(np.linalg.norm(saved[key][a,b])) for b in range(4)]
                          for a in range(4)] for key in arrays}
        receipt['items'].append(dict(external=external,archive=file_binding(path),
            count=count,norms=norms,seconds=perf_counter()-unit_start))
        atomic_json(record,receipt)
        print(json.dumps(receipt['items'][-1]),flush=True)
    cache.manifest['status'] = 'complete'
    cache.publish()
    receipt.update(status='complete',seconds=perf_counter()-started,
        peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    atomic_json(record,receipt)
    observed = consume_complete_result(args.output,args.manifest,slots,MATCHING,bounds)
    print(json.dumps(dict(**observed,seconds=receipt['seconds'])),flush=True)


if __name__ == '__main__':
    main()
