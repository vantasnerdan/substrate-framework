#!/usr/bin/env python3
"""Actual spin1 ik2 -> ik3 full-row finite-q operator consumer.

q_fractional=(0,0,1/4), fold=(0,0,-1). Original raw columns are used ONLY
as operator probes, never as a populated dynamics approximation. Restores the
accepted A25 background verbatim and exercises physical Kq, not a fitted
kernel. No H0 equilibrium, full material mean, or trajectory is asserted.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
sys.path[:0]=[str(ROOT/'src'),str(HERE.parent/'0025')]

import numpy as np
from substrate_framework.qe_saved_native import read_bound_bytes
from substrate_framework.qe_material_reference import load_material_inputs
from substrate_framework.qe_bloch_grid import build_bloch_grid, UNITARY_POLICY
from substrate_framework.qe_bloch_edge import BlochEdge
from substrate_framework.qe_bloch_finiteq_ionic import BlochFiniteQIonic
from material_reference_output import restore_background
from material_run_contract import runtime_context, implementation


def binding(path):
    raw=path.read_bytes()
    return dict(path=str(path.resolve()),bytes=str(len(raw)),sha256=hashlib.sha256(raw).hexdigest())


def consume(inputs, background, source_frame, target_frame, *, block_size=2):
    """All PW rows, source columns 0..2 and target columns 0..1, no populations.

    Frames must come from bound original-byte loaders. Returns raw complex
    actions, densities, gradients, Kq fields and pairing residuals, not a
    pass/fail convergence certificate. Root owns numerical execution/allocation.
    """
    s,t=source_frame.basis.slot,target_frame.basis.slot
    if (s.ispin!=1 or t.ispin!=1 or s.source_ik!=2 or t.source_ik!=3
            or tuple(s.mesh_index)!=(0,0,1) or tuple(t.mesh_index)!=(0,0,2)):
        raise ValueError('required original spin1 ik2 -> ik3 edge not supplied')
    if background.inputs is not inputs:
        raise ValueError('accepted background must share exact material inputs')
    edge=BlochEdge.create(inputs.reference,source_frame.basis,target_frame.basis,(0.,0.,.25))
    if not np.array_equal(edge.reciprocal_fold,(0,0,-1)):
        raise ValueError('known original edge fold changed')
    source=source_frame.coefficients[:,:3]
    target=target_frame.coefficients[:,:2]
    if source.shape[1]!=3 or target.shape[1]!=2:
        raise ValueError('three source and two target saved probe columns required')
    coupling=BlochFiniteQIonic(inputs,edge,block_size=block_size)
    bridge=coupling.bridge
    shape=bridge.source_fft.physical
    # Separable physical profile, flattened in the original F grid convention.
    profile=(np.cos(2*np.pi*np.arange(shape[0])[:,None,None]/shape[0])
             +1j*np.sin(2*np.pi*np.arange(shape[1])[None,:,None]/shape[1]))
    profile=np.broadcast_to(profile,shape).reshape(-1,order='F')
    direction=np.zeros((inputs.reference.nat,3),complex)
    direction[31,0]=1.; direction[0,1]=1j
    coefficient=np.array([[.7+.2j,-.3j,.11],[-.19,.4+.1j,.23j]])
    other=np.array([[.13j,.22,-.17+.04j],[.29,.12j,-.31]])
    weight=s.weight
    local=bridge.local_apply(profile,source)
    local_adj=bridge.local_adjoint_apply(profile,target)
    ionic=coupling.prepare(direction)
    gsource,gadj=ionic.apply(source),ionic.adjoint_apply(target)
    density,gradient=coupling.moments(target,source,coefficient,weight=weight)
    second=bridge.density_from_dyads(target,source,other,weight=weight)
    density_gradient=bridge.density_gradient(density)
    field=background.apply_kernel(density,q_bohr_inv=edge.q_cart_bohr_inv)
    second_field=background.apply_kernel(second,q_bohr_inv=edge.q_cart_bohr_inv)
    feedback=bridge.local_apply(field[:,0],source)
    feedback_adj=bridge.local_adjoint_apply(field[:,0],target)
    local_overlap=target.conj().T@local
    ionic_overlap=target.conj().T@gsource
    feedback_overlap=target.conj().T@feedback
    # Independent small factor contractions; no target-by-source density matrix.
    pairs=dict(
        local_density=np.asarray([bridge.density_pair(density,profile),weight*np.vdot(coefficient,local_overlap)]),
        ionic_gradient=np.asarray([np.vdot(direction.ravel(),gradient),np.conj(weight*np.vdot(coefficient,ionic_overlap))]),
        feedback_density=np.asarray([bridge.density_pair(density,field[:,0]),weight*np.vdot(coefficient,feedback_overlap)]),
        kernel_hermitian=np.asarray([bridge.grid_measure_bohr3*np.vdot(density,second_field),bridge.grid_measure_bohr3*np.vdot(field,second)]))
    # Exact finite-support integration by parts with i(q+G), retaining all
    # complex values. The arbitrary test vector uses both explicit channels.
    scalar=np.column_stack((profile,.3j*profile.conj()))
    left=0j; right=0j
    for axis,factor in enumerate((.7+.2j,-.4j,.11)):
        test=factor*scalar
        left+=bridge.grid_measure_bohr3*np.vdot(density_gradient[axis],test)
        divergence=bridge.density_fft.from_g(1j*bridge.density_fft.wavevectors[axis,:,None]*bridge.density_fft.to_g(test))
        right-=bridge.grid_measure_bohr3*np.vdot(density,divergence)
    pairs['density_gradient_adjoint']=np.asarray([left,right])
    arrays=dict(source_probe_columns=source,target_probe_columns=target,coefficient=coefficient,
        second_coefficient=other,profile=profile,direction=direction,
        local_action=local,local_adjoint_action=local_adj,ionic_action=gsource,ionic_adjoint_action=gadj,
        density_bohr3=density,positive_gradient_ry_bohr=gradient,density_gradient_bohr4=density_gradient,
        kernel_field_ry=field,second_density_bohr3=second,second_kernel_field_ry=second_field,
        feedback_action=feedback,feedback_adjoint_action=feedback_adj,
        q_fractional=edge.q_fractional,q_bohr_inv=edge.q_cart_bohr_inv,reciprocal_fold=edge.reciprocal_fold,
        **{name+'_pair':value for name,value in pairs.items()})
    report=dict(schema='P254-finiteq-material-operator-v1',
        scope='actual-material selected-column operator validation, all physical PW rows; NOT populated dynamics',
        source_basis_identity=source_frame.basis.identity,target_basis_identity=target_frame.basis.identity,
        source_original_binding=dict(s.source_binding),target_original_binding=dict(t.source_binding),
        source_dimensions=list(source_frame.coefficients.shape),target_dimensions=list(target_frame.coefficients.shape),
        selected_source_columns_zero_based=[0,1,2],selected_target_columns_zero_based=[0,1],
        q_fractional=edge.q_fractional.tolist(),q_bohr_inv=edge.q_cart_bohr_inv.tolist(),fold=edge.reciprocal_fold.tolist(),
        source_star_shift=s.reciprocal_shift.tolist(),target_star_shift=t.reciprocal_shift.tolist(),
        weight=weight,grid_measure_bohr3=bridge.grid_measure_bohr3,ispin=1,
        pairing_absolute_residuals={name:float(abs(value[0]-value[1])) for name,value in pairs.items()},
        action_adjoint_residuals=dict(local=float(np.linalg.norm(local_overlap-local_adj.conj().T@source)),
            ionic=float(np.linalg.norm(ionic_overlap-gadj.conj().T@source)),
            feedback=float(np.linalg.norm(feedback_overlap-feedback_adj.conj().T@source))),
        norms={name:float(np.linalg.norm(value)) for name,value in arrays.items()},
        radial_ranges={name:dict(original=inputs.tables[name].qmax_bohr_inv,used=table.qmax_bohr_inv) for name,table in coupling.tables.items()},
        limitations=['No selected-column population interpretation','Original source-state and selected density/gradient checkpoints do not establish full material equilibrium',
            'H0 source-state result remains pending; no H0 stationarity assertion','No actual full-grid trajectory or equilibrium/FDT claim',
            'No automatic Hermitization, clipping, tolerance certification or q-support repair','Cq and populated full-grid mean must be supplied for response trajectories'])
    return arrays,report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registration',type=Path,default=HERE/'finiteq-response-registration.json')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    context=runtime_context()
    with (args.output/'runtime-context.json').open('x') as stream:
        json.dump(context,stream,indent=2)
    registration_raw=args.registration.read_bytes()
    registration=json.loads(registration_raw)
    authorities={Path(entry['path']).name:entry for entry in registration['source_bindings']+registration['additional_original_bytes']}
    def authenticated_json(name):
        data,_=read_bound_bytes(authorities[name]); return json.loads(data)
    # Authenticate the receipts before using any nested binding or saved bytes.
    optical=authenticated_json('optical-current-sources.json')['accepted_material']
    source_metadata=authenticated_json('source-1-2.json')['metadata']
    source_state=authenticated_json('source-state.json')
    checkpoint=authenticated_json('slot-1-density-gradient.json')
    inputs_spec=source_metadata['inputs']
    if inputs_spec['xml']!=optical['inputs']['xml'] or inputs_spec['density']!=optical['inputs']['density'] or inputs_spec['upfs']!=optical['inputs']['upfs']:
        raise ValueError('accepted background and Bloch source input authorities differ')
    # All source-code authorities are checked before scientific consumption.
    for entry in registration['source_bindings']+registration['additional_original_bytes']:
        read_bound_bytes(entry)
    inputs=load_material_inputs(inputs_spec['xml'],inputs_spec['density'],inputs_spec['upfs'])
    raw,_=read_bound_bytes(optical['background_metadata']); background_metadata=json.loads(raw)
    raw,_=read_bound_bytes(background_metadata['archive'])
    with np.load(io.BytesIO(raw),allow_pickle=False) as saved:
        background=restore_background(inputs,{name:saved[name] for name in saved.files})
    del raw
    bindings={(int(spin),int(ik)):entry for spin,items in inputs_spec['wavefunctions'].items() for ik,entry in items.items()}
    grid=build_bloch_grid(inputs.reference,bindings,symmetry_policy=UNITARY_POLICY)
    ss=next(slot for slot in grid.slots if slot.ispin==1 and slot.mesh_index==(0,0,1))
    ts=next(slot for slot in grid.slots if slot.ispin==1 and slot.mesh_index==(0,0,2))
    if dict(ss.source_binding)!=source_state['metadata']['source'] or dict(ss.source_binding)!=checkpoint['metadata']['source']:
        raise ValueError('completed state/moment component receipts bind a different original source')
    # Existing loader authenticates each original WFC decode buffer, preserving
    # all raw rows/columns. No state receipt is treated as a replacement frame.
    source,target=grid.load_frame(ss),grid.load_frame(ts)
    if source.basis.identity!=source_state['metadata']['original_basis_identity'] or source.basis.identity!=checkpoint['metadata']['basis_identity']:
        raise ValueError('original state/moment basis identity differs')
    arrays,report=consume(inputs,background,source,target)
    report.update(registration=dict(path=str(args.registration.resolve()),bytes=str(len(registration_raw)),sha256=hashlib.sha256(registration_raw).hexdigest()),
        accepted_background=optical['background_metadata'],background_archive=background_metadata['archive'],
        completed_component_receipts={name:authorities[name] for name in ('source-state.json','slot-1-density-gradient.json')})
    report['runtime_context']=context
    report['implementation_sources']=implementation((Path(__file__),))
    archive=args.output/'finiteq-material.npz'
    with archive.open('xb') as stream: np.savez(stream,**arrays)
    report['archive']=binding(archive)
    with (args.output/'result.json').open('x') as stream: json.dump(report,stream,indent=2,allow_nan=False); stream.write('\n')
    print(json.dumps(dict(output=str(args.output),pairing_absolute_residuals=report['pairing_absolute_residuals'],action_adjoint_residuals=report['action_adjoint_residuals']),indent=2))


if __name__=='__main__':
    main()
