"""Genuine saved-U numerical consumer; never calls the U producer/integrator."""
import sys,json
from pathlib import Path
from time import perf_counter
import numpy as np
from general_moving_reference import (HERE,inputs,SavedMovingReference,atomic_json,file_binding,authenticate)
sys.path.insert(0,str(HERE.parent/'0050'))
from joint_gamma_adapter import ordered_coordinates

def main(bundle,output):
    began=perf_counter();c,baseline,binding,history=inputs()
    reference=SavedMovingReference(c,bundle,binding)
    electron=[]
    for a,t in enumerate(reference.times):
        gg,gl=reference.electron(a,0);rg,rl=reference.electron(0,a)
        equal=reference.electron(a,a)[1]
        electron.append(dict(time=float(t),all128_shape=list(gg.shape),greater_norm=float(np.linalg.norm(gg)),lesser_norm=float(np.linalg.norm(gl)),reverse_greater_norm=float(np.linalg.norm(rg)),reverse_lesser_norm=float(np.linalg.norm(rl)),density_trace_real=float(np.trace(-1j*equal,axis1=-2,axis2=-1).sum().real)))
    old_handoff=json.loads((HERE.parent/'0045/actual-complete-original-handoff-001.json').read_bytes())
    selected=[];comparisons=[]
    with np.load(authenticate(old_handoff['outputs']['raw_transitions']),allow_pickle=False) as old:
        for q in (0,1,3):
            qm=int(c.minus[q]);sq,sm=reference.transitions(q)
            assert sq.shape==(31,len(c.mechanical.generators[q]),len(c.mechanical.generators[q]))
            old_s=old['S_%02d'%q]
            for k,t in enumerate((0,1,10,30)):
                delta=sq[t]-old_s[k]
                comparisons.append(dict(q=q,time=t,maxabs=float(np.max(np.abs(delta))),frobenius=float(np.linalg.norm(delta))))
            records=[]
            for a,t in enumerate(reference.times):
                b=max(0,a-1)
                dg,dl=reference.boson(q,a,b)
                pg,pl=reference.ordered(q,a,b,'positions')
                records.append(dict(time=float(t),s=float(reference.times[b]),ionic_greater_norm=float(np.linalg.norm(dg)),ionic_lesser_norm=float(np.linalg.norm(dl)),positions_greater_norm=float(np.linalg.norm(pg)),positions_lesser_norm=float(np.linalg.norm(pl)),positions_shape=list(pg.shape)))
            selected.append(dict(q=q,minus=qm,S_shape=list(sq.shape),S_maxabs=float(np.max(np.abs(sq))),position_consumption=records))
    # Existing A50 full114-coordinate callback consumes the new provider unchanged.
    g,l=ordered_coordinates(reference,30,17)
    direct=reference.ordered(0,30,17,'positions')
    gamma_gap=max(float(np.max(np.abs(g-direct[0]))),float(np.max(np.abs(l-direct[1]))))
    assert g.shape==(114,114) and gamma_gap<1e-11
    # Saved callback selection preserves exact requested fractional semantics:
    # absent fractions reject rather than silently interpolating four/31 nodes.
    rejected=False
    try:SavedMovingReference(c,bundle,binding,times=[.5])
    except ValueError as error:rejected='requested time absent' in str(error)
    assert rejected and reference.u_integration_calls==0
    result=dict(schema='p254-a52-genuine-saved-reference-consumer-v1',status='complete',bundle=file_binding(bundle),
        electron_all128_all31=electron,selected_actual_q=selected,old_four_node_S_comparison=comparisons,
        old_node_policy='same direct scipy.linalg.expm(t*Aq), original full mechanical generators/preparation; raw differences only, no forced matches',
        A50_existing_ordered_coordinates=dict(shape=list(g.shape),time=30,s=17,maxabs_difference=gamma_gap,greater_norm=float(np.linalg.norm(g)),lesser_norm=float(np.linalg.norm(l))),
        saved_U_integration_calls=reference.u_integration_calls,missing_fraction_rejected=rejected,
        S_expm_calls=reference.expm_calls,S_cache_hits=reference.cache_hits,S_cached_bytes=reference.cache_size,S_cache_cap=reference.cache_bytes,
        elapsed_seconds=perf_counter()-began,reference_bindings=reference.bindings,
        scope='numerical consumption of saved uniform reference only; no uniform selfenergy, physicalaccuracy, feedback or closure claim')
    atomic_json(output,result);reference.close();print(json.dumps(result),flush=True)
if __name__=='__main__':main(sys.argv[1],sys.argv[2])
