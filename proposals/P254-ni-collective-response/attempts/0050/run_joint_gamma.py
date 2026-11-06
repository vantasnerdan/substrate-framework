"""Compute actual all114 Gamma component and extract same-sector channel."""
import json
from pathlib import Path
import sys
import time
import numpy as np
from joint_gamma_adapter import common_pair,vertices,coordinate_indices
from joint_material_inputs import load_inputs,file_binding
from moving_reference_inputs import atomic_json
from represented_electron_channel import RepresentedElectronChannel,CUTOVER


def main(output,mode,common_archive=None):
    output=Path(output)
    if mode=='common':
        output.mkdir(exist_ok=False)
    start=time.monotonic()
    c,baseline,reference,bindings=load_inputs()
    if mode=='represented':
        if common_archive is not None:
            output.mkdir(exist_ok=False)
            atomic_json(output/'actual-inputs.json',bindings)
        return extract(output,c,baseline,reference,file_binding(common_archive or output/'common.npz'),start)
    atomic_json(output/'actual-inputs.json',bindings)
    atomic_json(output/'input-loaded.json',dict(elapsed=time.monotonic()-start,
        gamma_state_size=c.ns*c.m*c.m,coordinates=c.d+c.f,full_original_mechanical_size=len(c.mechanical.generators[0])))
    print('A50 full original Gamma inputs loaded',flush=True)
    common={}
    for a,t in enumerate(reference.times):
        for b,s in enumerate(reference.times):
            joint=common_pair(c,reference,a,b)
            ionic=common_pair(c,reference,a,b,ionic_only=True)
            for prefix,values in (('joint_',joint),('ionic_',ionic)):
                for key,value in values.items():
                    name=prefix+key
                    if name not in common:
                        common[name]=np.empty((4,4,*value.shape),complex)
                    common[name][a,b]=value
            atomic_json(output/'common-boundary.json',dict(completed_pair=[a,b],elapsed=time.monotonic()-start))
            print(json.dumps(dict(stage='common',pair=[float(t),float(s)],elapsed=time.monotonic()-start)),flush=True)
    np.savez(output/'common.npz',**common)
    common_binding=file_binding(output/'common.npz')
    del common
    reference.close()
    return common_binding


def extract(output,c,baseline,reference,common_binding,start):
    channel=RepresentedElectronChannel(c,baseline)
    zero_duals=[];responses=np.zeros((4,4,c.d+c.f,c.d+c.f),complex);work=[]
    for a,t in enumerate(reference.times):
        before=time.monotonic()
        zero=np.empty((*channel.shape,c.d+c.f),complex)
        batches=[]
        for first in range(0,c.d+c.f,19):
            indices=np.arange(first,min(first+19,c.d+c.f))
            pulled,receipt=channel.pullback(float(t),reference.times[:a+1],indices)
            zero[:,:,:,indices]=pulled[0.]
            for b in range(a+1):
                responses[a,b,indices,:]=channel.response_from_dual(
                    pulled[float(reference.times[b])],float(reference.times[b]))
            receipt.update(elapsed=time.monotonic()-before,cumulative_elapsed=time.monotonic()-start)
            batches.append(receipt)
            atomic_json(output/('represented-boundary-%d-batch-%d.json'%(a,first)),receipt)
            print(json.dumps(dict(stage='represented',**receipt)),flush=True)
            del pulled
        zero_duals.append(zero)
        work.append(dict(terminal=float(t),batches=batches,elapsed=time.monotonic()-before))
    ab=np.empty_like(responses)
    for a in range(4):
        for b in range(4):
            ab[a,b]=channel.ordered_noise(zero_duals[a],zero_duals[b])
    ba=ab.swapaxes(0,1).swapaxes(-1,-2).copy()
    np.savez(output/'represented.npz',pi_retarded=responses,pi_greater=-1j*ab,pi_lesser=-1j*ba,
             ordered_ab=ab,ordered_ba=ba)
    forward,reverse=vertices(c)
    diagnostics=dict(stored_reverse_minus_forward_transpose_norm=float(np.linalg.norm(reverse-forward.swapaxes(-1,-2))),
        rho0_hermiticity_defect=float(np.linalg.norm(c.rho0-c.rho0.conj().swapaxes(-1,-2))),
        raw_mean_rho_h_commutator_norms=[float(np.linalg.norm(h@rho-rho@h)) for rho,h in [baseline(float(t)) for t in reference.times]],
        mechanical=c.mechanical.coefficient_diagnostics(),
        retained_coordinate_indices=coordinate_indices(c).tolist())
    record=dict(schema='p254-a50-original-joint-gamma-component-v1',status='complete',
        times=reference.times.tolist(),common_archive=common_binding,represented_archive=file_binding(output/'represented.npz'),
        inputs=file_binding(output/'actual-inputs.json'),work=work,elapsed_seconds=time.monotonic()-start,
        diagnostics=diagnostics,cutover=CUTOVER,
        scope='original all111 ionic plus3 canonicalfield Gamma; all128 CAR slots/all12 bands; all16 original timepairs; COMPONENT ONLY',
        represented_scope='density/K-dressed linearized electron channel about original moving mean with prescribed mechanical forcing; NOT undressed A41 bubble',
        approximation=reference.label,limits=['no allq closure','no selfconsistent loop','no XC double-count prescription','no physical rates/fullE','no passive/complement dynamic response/noise'],
        findings=[])
    atomic_json(output/'result.json',record)
    reference.close()
    print(json.dumps(dict(status='complete',result=str(output/'result.json'),elapsed=record['elapsed_seconds'])),flush=True)


if __name__=='__main__':main(sys.argv[1],sys.argv[2],sys.argv[3] if len(sys.argv)>3 else None)
