"""Produce genuine moving transitions; consume through unchanged A41 public reader."""
import argparse
import json
from pathlib import Path
from time import perf_counter
import zipfile
import numpy as np
from moving_reference_inputs import (inputs, own_sources, TIMES, CONTRACT,
    atomic_json, file_binding, authenticate, SavedTransitions)
from moving_reference_propagation import AffineSavedHamiltonian, integrate_u, checkpoint_u, bare_sequence


def produce(output):
    output=Path(output);output.mkdir(exist_ok=False)
    started=perf_counter()
    sources=own_sources()
    c,baseline,coefficient_binding,history=inputs()
    history_bindings=[CONTRACT,history,*sources]
    h=AffineSavedHamiltonian(baseline)
    print(json.dumps(dict(phase='actual_affine_h_ready',seconds=perf_counter()-started,
        segments=len(h.coefficients))),flush=True)
    u=np.broadcast_to(np.eye(c.m,dtype=complex),(c.ns,c.m,c.m)).copy()
    us=[u.copy()];work=[];checkpoints=[];diagnostics=[]
    u,checkpoint=checkpoint_u(output,0,u,history_bindings,[]);checkpoints.append(checkpoint)
    for left,right in zip(TIMES[:-1],TIMES[1:]):
        before=perf_counter()
        u,w=integrate_u(h,left,right,u)
        w.update(left=left,right=right,seconds=perf_counter()-before)
        work.append(w)
        # Every next integration ACTUALLY consumes authenticated persisted U.
        u,checkpoint=checkpoint_u(output,right,u,history_bindings,work.copy())
        checkpoints.append(checkpoint);us.append(u.copy())
        rho,direct_h=baseline(right)
        reconstructed=u@c.rho0@u.conj().swapaxes(-1,-2)
        diagnostics.append(dict(time_inverse_ry=right,
            U_unitarity_frobenius=float(np.linalg.norm(u@u.conj().swapaxes(-1,-2)-np.eye(c.m))),
            rho_from_U_minus_original_saved_rho_frobenius=float(np.linalg.norm(reconstructed-rho)),
            rho_from_U_minus_original_saved_rho_maxabs=float(np.max(np.abs(reconstructed-rho))),
            saved_rho_hermiticity_frobenius=float(np.linalg.norm(rho-rho.conj().swapaxes(-1,-2))),
            h_hermiticity_frobenius=float(np.linalg.norm(direct_h-direct_h.conj().swapaxes(-1,-2))),
            affine_h_interpolation_direct_maxabs=float(np.max(np.abs(h(right)-direct_h)))))
        print(json.dumps(dict(phase='actual_U_boundary',work=w,diagnostics=diagnostics[-1],
            checkpoint=checkpoint,elapsed_seconds=perf_counter()-started)),flush=True)
        if right==1. and perf_counter()-started>360:
            atomic_json(output/'allocation-trigger.json',dict(decision='stop at authenticated U1 boundary; change approach required',
                elapsed_seconds=perf_counter()-started,trigger='3x expected120s first checkpoint'))
            raise RuntimeError('first checkpoint exceeded3x expected cost; no unchanged continuation')
    archive=output/'transitions.npz';s_diagnostics=[]
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_STORED,allowZip64=True) as stream:
        with stream.open('U.npy','w',force_zip64=True) as member:
            np.lib.format.write_array(member,np.asarray(us),allow_pickle=False)
        for q in range(c.nq):
            sq=bare_sequence(c.mechanical.generators[q],TIMES)
            if not np.isfinite(sq).all():
                raise FloatingPointError('nonfinite bare original S; no clipping permitted')
            with stream.open('S_%02d.npy'%q,'w',force_zip64=True) as member:
                np.lib.format.write_array(member,sq,allow_pickle=False)
            # Preserve raw generator/invariant diagnostics, including unstable directions.
            qm=int(c.minus[q]);sm=sq if q==qm else bare_sequence(c.mechanical.generators[qm],TIMES)
            jq=c.mechanical.j[q]
            s_diagnostics.append(dict(q=q,phase_space=len(jq),
                S_maxabs=[float(np.max(np.abs(s))) for s in sq],
                paired_symplectic_frobenius=[float(np.linalg.norm(s@jq@sm[a].T-jq)) for a,s in enumerate(sq)],
                ordered_CCR_frobenius=[float(np.linalg.norm(s@c.mechanical.w0[q]@sm[a].T-
                    (sm[a]@c.mechanical.w0[qm]@s.T).T-1j*jq)) for a,s in enumerate(sq)]))
            print(json.dumps(dict(phase='actual_bare_S',q=q,elapsed_seconds=perf_counter()-started)),flush=True)
    bundle=dict(schema='p254-electron-ion-transitions-v1',coefficient_manifest=coefficient_binding,
        times_inverse_ry=TIMES,archive=file_binding(archive),
        physical_history_bindings=history_bindings+checkpoints,
        approximation='timeordered U from original saved moving mean h_baseline(t); original bare Gaussian S from constant original mechanical generators; NO electron/cubic/collision feedback')
    atomic_json(output/'bundle.json',bundle)
    reference=SavedTransitions(c,TIMES,output/'bundle.json',coefficient_binding)
    try:
        validation=[]
        for q in range(c.nq):
            sq,sm=reference.transitions(q)
            greater,lesser=reference.boson(q,1,2)
            validation.append(dict(q=q,S_shape=list(sq.shape),boson_greater_norm=float(np.linalg.norm(greater)),
                boson_lesser_norm=float(np.linalg.norm(lesser))))
        gg,gl=reference.electron(1,2)
        consumer=dict(reader='unchanged A41 SavedTransitions',all64q=validation,
            electron_greater_norm=float(np.linalg.norm(gg)),electron_lesser_norm=float(np.linalg.norm(gl)))
    finally:reference.close()
    record=dict(schema='p254-a45-moving-reference-producer-v1',status='complete',
        bundle=file_binding(output/'bundle.json'),source_bindings=sources,
        baseline=history,coefficient_manifest=coefficient_binding,checkpoints=checkpoints,
        actual_work=work,hamiltonian_construction_seconds=h.construction_seconds,
        raw_U_rho_diagnostics=diagnostics,raw_S_diagnostics=s_diagnostics,
        actual_persisted_consumer=consumer,elapsed_seconds=perf_counter()-started,
        settings=dict(method='DOP853',rtol=2e-8,atol=2e-10),
        recovery='U time0/1/10/30 boundaries actually reloaded before dependent intervals; adaptive DOP853 history lost at boundary; unfinished interval replays',
        accuracy='raw diagnostics only; no repair/projection/clipping, no rigorous accuracy claim',
        scope='original bare Gaussian moving-reference transitions, not selfconsistent KB or collision/cubic feedback')
    atomic_json(output/'result.json',record)
    return file_binding(output/'result.json')


def recover(kernel_output,bundle,output):
    # Existing A41 compute is the real immutable q-result consumer. No old CLI
    # replacement and no kernel implementation changes are needed here.
    from material_electron_ion_inputs import load_material
    from run_material_electron_ion_kernel import compute
    from reduced_projection_cache import ProjectionCache
    kernel_output=Path(kernel_output)
    original=json.loads((kernel_output/'result.json').read_bytes())
    c,_,binding=load_material()
    reference=SavedTransitions(c,TIMES,bundle,binding)
    manifest=json.loads((kernel_output/'q-kernels/manifest.json').read_bytes())
    cache=ProjectionCache(kernel_output/'q-kernels',manifest['specification'])
    if cache.manifest['specification']['transition_bindings']!=reference.bindings:
        raise ValueError('kernel recovery transition identity differs')
    before=file_binding(cache.path)
    try:
        arrays,count,reused,defects=compute(c,reference,[0,64],[0,1,2],cache)
        maximum={}
        with np.load(authenticate(original['output']),allow_pickle=False) as saved:
            for key,value in arrays.items():
                if not np.array_equal(value,saved[key]):
                    raise AssertionError('persisted q recovery differs: '+key)
                maximum[key]=float(np.max(np.abs(value)))
        if count!=0 or reused!=original['analytic_edge_pair_evaluations']:
            raise AssertionError('complete exact-input q recovery was not exercised')
        record=dict(schema='p254-a45-moving-kernel-recovery-v1',passed=True,
            original_result=file_binding(kernel_output/'result.json'),
            recovered_new_pairs=count,recovered_saved_pairs=reused,
            all_arrays_byte_equal=True,computed_maxabs=maximum,defects=defects,
            original_manifest=before,consumed_manifest=file_binding(cache.path),
            source_bindings=own_sources(),
            recovery='complete committed q boundaries only; unfinished q is not retained; original A41 compute unchanged')
        atomic_json(output,record)
    finally:reference.close()
    return file_binding(output)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['produce','recover'])
    parser.add_argument('--output',required=True)
    parser.add_argument('--kernel-output');parser.add_argument('--bundle')
    args=parser.parse_args()
    result=produce(args.output) if args.mode=='produce' else recover(args.kernel_output,args.bundle,args.output)
    print(json.dumps(result),flush=True)
