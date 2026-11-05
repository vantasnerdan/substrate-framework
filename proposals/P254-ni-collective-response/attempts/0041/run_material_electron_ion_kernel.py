"""Evaluate actual all-internal-sum Ni/H Sigma/Pi, not a lifetime or KB solver.

All original 0/1/10/30 inverse-Ry time pairs are retained. External electron
slots and external Pi q are declared by the caller; these do NOT truncate any
internal momentum/spin sum. Output includes full12x12 Sigma and111x111 Pi.
Numerical payloads stay outside Git; output creation is exclusive.
"""
import argparse
from dataclasses import fields
import json
from pathlib import Path
from time import perf_counter
import resource

import numpy as np

from electron_ion_functional import evaluate_components
from material_electron_ion_inputs import (HERE, ROOT, DEFAULT_MANIFEST,
    FrozenReference, SavedTransitions, load_material, edge_pair, reference_diagnostics,
    file_binding)
from substrate_framework.partitioned_q_storage import release_mapped_pages
from reduced_projection_cache import ProjectionCache


MATCHING = dict(
    interaction='H_ei=sum u_q cdag_(k+q) g(q,k)c_k/sqrt(N); g=dH/du, force=-g',
    measure='global canonical rho; explicit two spins and64 k; no slot weights; two vertices give1/N',
    boson='connected original ordered Cartesian u fluctuations; no mean dyad, FDT or oscillator rescale',
    curvature='A28 material_matched_reference.assemble_frozen: C_II(q)+absolute frozen <d2H_e/du2>; contact is NOT a bubble',
    insertion='these raw kernels are NOT fed back here; replace the explicit-electron tangent polarization if closing D, never stack it',
    mean='existing h0 includes restored K(original density-saved background); no new Hartree/tadpole or equilibrium refill',
    dft='second-Born residual electron-ion correction about this fixed selected DFT mean; no claim of ab-initio XC double-counting removal',
    initial='uncorrelated electron/ion product Gaussian IVP at0; connected electronic bilinear force noise comes from SAME Pi lesser',
    missing='selfconsistent KB evolution, moving-mean U generation, active/passive/complement loops, cubic-ion/e-e collisions and physical rate/ranking remain open')


def indices(text, count):
    values = list(range(count)) if text == 'all' else [int(x) for x in text.split(',')]
    if not values or len(set(values)) != len(values) or min(values) < 0 or max(values) >= count:
        raise ValueError('distinct external indices within original axes required')
    return values


def _compute_q(c, reference, q, sigma_slots, want_pi, electrons):
    """One recoverable q contribution, without discarded unrequested diagrams."""
    nt = len(reference.times)
    sigma = {key: np.zeros((nt,nt,len(sigma_slots),c.m,c.m),complex)
             for key in ('greater','lesser','retarded')}
    pi = ({key: np.zeros((nt,nt,c.d,c.d),complex)
           for key in ('greater','lesser','retarded')} if want_pi else {})
    slot_positions = {slot: j for j,slot in enumerate(sigma_slots)}
    defects = dict(vertex_paired_adjoint_max=0., boson_ordered_ccr_max=0.,
                   boson_paired_symplectic_max=0.)
    data = c.q(q)
    handles = [getattr(data,name) for name in ('ionic','ionic_reverse','density','potential','kernel')]
    sources = range(c.ns) if want_pi else [
        int(c.target[int(c.minus[q]),slot]) for slot in sigma_slots]
    count = 0
    try:
        defects['vertex_paired_adjoint_max'] = float(
            np.linalg.norm(data.ionic_reverse-data.ionic.conj()))
        sq, sm = reference.transitions(q)
        jq = c.mechanical.j[q]
        qm = int(c.minus[q])
        for a in range(nt):
            defects['boson_paired_symplectic_max'] = max(
                defects['boson_paired_symplectic_max'],
                float(np.linalg.norm(sq[a]@jq@sm[a].T-jq)))
            ordered = sq[a]@c.mechanical.w0[q]@sm[a].T
            reverse = sm[a]@c.mechanical.w0[qm]@sq[a].T
            defects['boson_ordered_ccr_max'] = max(
                defects['boson_ordered_ccr_max'],
                float(np.linalg.norm(ordered-reverse.T-1j*jq)))
            for b in range(nt):
                boson = reference.boson(q,a,b)
                for source in sources:
                    target = int(c.target[q,source])
                    pair = edge_pair(data,source,target,electrons[a,b],electrons[b,a],boson)
                    wanted = tuple(name for name, keep in (
                        ('sigma', target in slot_positions), ('pi', want_pi)) if keep)
                    kernel = evaluate_components(pair,n=c.n,
                        t=float(reference.times[a]),s=float(reference.times[b]),components=wanted)
                    count += 1
                    if target in slot_positions:
                        for component in sigma:
                            sigma[component][a,b,slot_positions[target]] += kernel['sigma_'+component]
                    for component in pi:
                        pi[component][a,b] += kernel['pi_'+component]
    finally:
        for handle in handles:
            release_mapped_pages(handle)
    arrays = {'sigma_'+key: value for key,value in sigma.items()}
    arrays.update({'pi_'+key: value for key,value in pi.items()})
    if not all(np.isfinite(value).all() for value in arrays.values()):
        raise FloatingPointError('nonfinite q kernels; unstable reference was NOT clipped')
    return arrays, count, defects


def compute(c, reference, sigma_slots, pi_q, cache):
    """Compose authenticated immutable q results in their original sum order."""
    nt = len(reference.times)
    sigma = {key: np.zeros((nt,nt,len(sigma_slots),c.m,c.m),complex)
             for key in ('greater','lesser','retarded')}
    pi = {key: np.zeros((nt,nt,len(pi_q),c.d,c.d),complex)
          for key in ('greater','lesser','retarded')}
    q_positions = {q: j for j,q in enumerate(pi_q)}
    electrons = {(a,b): reference.electron(a,b) for a in range(nt) for b in range(nt)}
    count, reused_count = 0, 0
    defects = dict(vertex_paired_adjoint_max=0., boson_ordered_ccr_max=0.,
                   boson_paired_symplectic_max=0.)
    for q in range(c.nq):
        key = 'q/%02d' % q
        if cache.has(key):
            reused_count += cache.manifest['items'][key]['metadata']['edge_pair_evaluations']
        else:
            started = perf_counter()
            arrays, q_count, q_defects = _compute_q(
                c,reference,q,sigma_slots,q in q_positions,electrons)
            cache.path.with_suffix('.json.prev').write_bytes(cache.path.read_bytes())
            cache.put(key,arrays,metadata=dict(q_index=q,
                edge_pair_evaluations=q_count,diagnostics=q_defects),
                timings=dict(seconds=perf_counter()-started))
            count += q_count
            del arrays
        # The actual final consumer reads committed q outputs, not RAM echoes.
        saved = cache.get(key)
        metadata = cache.manifest['items'][key]['metadata']
        for component in sigma:
            sigma[component] += saved['sigma_'+component]
        if q in q_positions:
            for component in pi:
                pi[component][:,:,q_positions[q]] = saved['pi_'+component]
        for name in defects:
            defects[name] = max(defects[name],metadata['diagnostics'][name])
        print(json.dumps(dict(completed_q=q,archive=cache.manifest['items'][key]['archive'],
            newly_computed_edge_pairs=count,reused_edge_pairs=reused_count)),flush=True)
        del saved
    arrays = {'sigma_'+key: value for key,value in sigma.items()}
    arrays.update({'pi_'+key: value for key,value in pi.items()})
    arrays.update(times_inverse_ry=reference.times,sigma_slots=np.asarray(sigma_slots,int),
        pi_q=np.asarray(pi_q,int),slot_spin=c.slot_spin,slot_k=c.slot_k,
        target=c.target,minus=c.minus)
    cache.manifest['status'] = 'complete'
    cache.publish()
    return arrays, count, reused_count, defects


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',type=Path,default=Path(DEFAULT_MANIFEST))
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--resume',action='store_true',
                        help='reuse only exact-input committed q kernels in this output')
    parser.add_argument('--sigma-slots',default='0,64',help='external slots, or all; internal q always all64')
    parser.add_argument('--pi-q',default='0,1,2',help='external q, or all; internal slots always all128')
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--frozen-reference',action='store_true')
    mode.add_argument('--saved-transitions',type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    if output == ROOT or ROOT in output.parents:
        raise ValueError('material numerical output must remain outside repository')
    sigma_slots, pi_q = indices(args.sigma_slots,128), indices(args.pi_q,64)
    output.mkdir(parents=True,exist_ok=args.resume)
    started = perf_counter()
    c, manifest, binding = load_material(args.manifest)
    times = [0.,1.,10.,30.]
    reference = (FrozenReference(c,times) if args.frozen_reference else
                 SavedTransitions(c,times,args.saved_transitions,binding))
    source_bindings = [file_binding(HERE/name) for name in ('electron_ion_functional.py',
        'material_electron_ion_inputs.py','run_material_electron_ion_kernel.py')]
    cache = ProjectionCache(output/'q-kernels',dict(
        kind='electron-ion-two-time-kernel-components',coefficient_manifest=binding,
        reference=reference.label,transition_bindings=reference.bindings,
        times_inverse_ry=times,sigma_slots=sigma_slots,pi_q=pi_q,
        source_bindings=source_bindings))
    try:
        arrays,count,reused_count,defects = compute(c,reference,sigma_slots,pi_q,cache)
        path = output/'kernels.npz'
        with path.open('xb') as stream:
            np.savez(stream,**arrays)
        record = dict(schema='p254-material-electron-ion-kernel-v1',status='complete',
            reference=reference.label,coefficient_manifest=binding,
            original_boundary=manifest['specification']['original_boundary'],
            transition_bindings=reference.bindings,matching=MATCHING,
            dimensions=dict(bands=12,feature_family=27,Cartesian=111,explicit_spin_k=128,q=64),
            times_inverse_ry=times,requested_physical_times_inverse_ry=[1.,10.,30.],
            sigma_external_slots=sigma_slots,pi_external_q=pi_q,
            internal_sums=dict(sigma_all64q=True,pi_all128spin_k=True,implicit_spin_factor=False),
            units=dict(G='dimensionless',D='bohr^2',g='Ry/bohr',Sigma='Ry^2',Pi='Ry^2/bohr^2',time='inverse Ry'),
            output=file_binding(path),actual_edge_pair_evaluations=count,
            reused_edge_pair_evaluations=reused_count,
            completed_q_checkpoint=file_binding(cache.path),
            analytic_edge_pair_evaluations=16*(128*len(pi_q)+len(sigma_slots)*(64-len(pi_q))),
            diagnostics=dict(reference=reference_diagnostics(reference),**defects),
            source_bindings=source_bindings,
            elapsed_seconds=float(perf_counter()-started),
            maxrss_kib=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
            selfconsistent=False,moving_mean_reconstructed=False,lifetime_or_rate_computed=False,
            full_material_goal_complete=False)
        with (output/'result.json').open('x') as stream:
            json.dump(record,stream,indent=2,sort_keys=True,allow_nan=False)
            stream.write('\n')
        print(json.dumps(dict(result=str(output/'result.json'),output=record['output'])))
    finally:
        if isinstance(reference,SavedTransitions):
            reference.close()


if __name__ == '__main__':
    main()
