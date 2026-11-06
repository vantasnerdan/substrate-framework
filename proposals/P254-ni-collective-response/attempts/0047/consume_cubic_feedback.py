"""Actual full37 source consumer and only the new generalized-boundary oracle.

Old physical Ewald/NC/Gaussian/Fock oracles are carried, not replayed. One
transient111^3 reference expands the already-qualified donor for an independent
entrywise contraction check; no Phi cube or external/time kernel bank is saved.
"""
import argparse
import json
from pathlib import Path
import sys
from time import perf_counter
import numpy as np

HERE=Path(__file__).resolve().parent;A43=HERE.parent/'0043'
sys.path.insert(0,str(A43))
from cubic_feedback import CompactCubic
from reduced_projection_cache import ProjectionCache,atomic_json
from substrate_framework.streamed_material_state import authenticate,file_binding


def relative(value,expected):
    return float(np.linalg.norm(np.asarray(value)-np.asarray(expected))/max(1.,float(np.linalg.norm(expected))))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();started=perf_counter()
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=False)
    model=CompactCubic.from_source(args.source)
    if (model.n,model.nat,model.d)!=(64,37,111):
        raise ValueError('actual full37 original material source required')
    specification=dict(kind='source-complete-compact-cubic-feedback-boundary-consumer',source=model.source_binding,
        source_bindings=[file_binding(HERE/'cubic_feedback.py'),file_binding(__file__)],
        all_atoms=37,all_coordinates=111,all_internal_q=64,
        no_persistent_phi_tensor_bank=True,no_external_time_kernel_bank=True)
    cache=ProjectionCache(output/'feedback',specification)
    rng=np.random.default_rng(4702026)
    def random(shape,scale=1.):return scale*(rng.normal(size=shape)+1j*rng.normal(size=shape))
    means=np.zeros((64,111),complex);dmeans=np.zeros_like(means);cot=np.zeros_like(means)
    for q in (0,1,3):means[q]=random((111,),.02)
    for q in (0,1,2,3):
        dmeans[q]=random((111,),.01);cot[q]=random((111,))
    covariances=random((64,111,111),.002)
    dcovariances=random((64,111,111),.001)
    x,y,z=random((3,111),.02);dx,dy,dz=random((3,111),.01)
    point_cot=random((111,));scalar_cot=.7+.6j
    probe_inputs=dict(means=means,dmeans=dmeans,cotangent=cot,covariances=covariances,
        dcovariances=dcovariances,x=x,y=y,z=z,dx=dx,dy=dy,dz=dz,point_cotangent=point_cot,
        scalar_cotangent=np.asarray(scalar_cot))
    raw_input_bytes=sum(value.nbytes for value in probe_inputs.values())
    projected_output_bytes=64*111*111*16+6*64*111*16+8*111*16
    if raw_input_bytes+projected_output_bytes+(8<<20)>64<<20:
        raise ValueError('declared finite input/action payload would exceed64MiB')
    cache.put('probe-inputs',probe_inputs,metadata=dict(purpose='complex linearization probes, not a modified electronic source state',
        active_mean_q=[0,1,3],active_mean_tangent_q=[0,1,2,3],all_covariance_q=list(range(64)),all_coordinates=111))
    actual=cache.get('probe-inputs')
    means=actual['means'];dmeans=actual['dmeans'];cot=actual['cotangent']
    covariances=actual['covariances'];dcovariances=actual['dcovariances']
    x,y,z=actual['x'],actual['y'],actual['z'];dx,dy,dz=actual['dx'],actual['dy'],actual['dz']
    point_cot=actual['point_cotangent']
    qids=[1,4,15];qs=model.indices[qids]
    phi=model.ionic.block(qs)
    onsite=(model.local+model.nonlocal_term)/np.sqrt(model.n)
    for atom in range(37):
        sl=slice(3*atom,3*atom+3);phi[sl,sl,sl]+=onsite[atom]
    point_w=covariances[0];point_dw=dcovariances[0]
    expected_force=-.5*np.einsum('abc,bc->a',phi,point_w+y[:,None]*z[None,:],optimize=True)
    expected_scalar=np.einsum('abc,a,b,c->',phi,x,y,z,optimize=True)/6
    expected_jvp=-.5*np.einsum('abc,bc->a',phi,point_dw+dy[:,None]*z[None,:]+y[:,None]*dz[None,:],optimize=True)
    point_force=model.mean_force(qs,y,z,point_w)
    point_scalar=model.ordered_scalar(qs,x,y,z)
    point_jvp=model.mean_force_jvp(qs,y,z,dy,dz,point_dw)
    point_dual=model.mean_force_transpose(qs,y,z,point_cot)
    expected_dual=(-.5*np.einsum('abc,a,c->b',phi,point_cot,z,optimize=True),
        -.5*np.einsum('abc,a,b->c',phi,point_cot,y,optimize=True),
        -.5*np.einsum('abc,a->bc',phi,point_cot,optimize=True))
    errors=dict(point_force_vs_qualified_dense_reference=relative(point_force,expected_force),
        ordered_scalar_vs_qualified_dense_reference=relative(point_scalar,expected_scalar),
        point_jvp_vs_independent_entrywise_polynomial=relative(point_jvp,expected_jvp),
        point_true_transpose_vs_entrywise=max(relative(value,expected) for value,expected in zip(point_dual,expected_dual)))
    lhs=np.sum(point_cot*point_jvp)
    rhs=np.sum(point_dual[0]*dy)+np.sum(point_dual[1]*dz)+np.sum(point_dual[2]*point_dw)
    errors['point_bilinear_transpose_identity']=relative(lhs,rhs)
    h=.125
    finite_point=(model.mean_force(qs,y+h*dy,z+h*dz,point_w+h*point_dw)
        -model.mean_force(qs,y-h*dy,z-h*dz,point_w-h*point_dw))/(2*h)
    errors['point_central_finite_direction']=relative(point_jvp,finite_point)
    scalar_direction=model.scalar_jvp(qs,x,y,z,dx,dy,dz)
    scalar_dual=model.scalar_transpose(qs,x,y,z,scalar_cot)
    errors['scalar_bilinear_transpose_identity']=relative(scalar_cot*scalar_direction,
        np.sum(scalar_dual[0]*dx)+np.sum(scalar_dual[1]*dy)+np.sum(scalar_dual[2]*dz))
    mutations=dict(point_wrong_covariance_transpose=float(np.linalg.norm(point_force-model.mean_force(qs,y,z,point_w.T))/max(1e-30,np.linalg.norm(point_force))),
        point_hilbert_adjoint=float(abs(lhs-(np.sum(point_dual[0].conj()*dy)+np.sum(point_dual[1].conj()*dz)+np.sum(point_dual[2].conj()*point_dw)))/max(1e-30,abs(lhs))),
        missing_canonical_one_over_sqrtN=float(np.linalg.norm(point_force*8-point_force)/max(1e-30,np.linalg.norm(point_force))))
    del phi
    print(json.dumps(dict(checkpoint='actual111-coordinate general-q scalar/force/JVP/true-transpose',elapsed_seconds=perf_counter()-started,errors=errors,mutations=mutations)),flush=True)
    if max(errors.values())>=2e-12 or any(value<=1e-3 for value in mutations.values()):
        atomic_json(output/'failed-boundary-findings.json',dict(errors=errors,mutations=mutations,source=model.source_binding))
        raise AssertionError(dict(errors=errors,mutations=mutations))
    manifest_path=Path(model.source_record['coefficient_manifest']['path'])
    manifest=json.loads(Path(authenticate(model.source_record['coefficient_manifest'])).read_bytes())
    mechanical_binding=manifest['items']['mechanical']['archive']
    with np.load(authenticate(mechanical_binding),allow_pickle=False) as saved_mechanical:
        original_ionic_ordered=saved_mechanical['ionic_ordered']
    original_uu=original_ionic_ordered[:111,:111]
    physical_mean=np.zeros((64,111),complex)
    for k in range(64):
        physical_qs=model.indices[[0,k,int(model.minus[k])]]
        contribution=-.5*model.contract(physical_qs,covariance=original_uu)
        physical_mean[0]+=contribution
        print(json.dumps(dict(checkpoint='actual saved-Gaussian full111 mean-force route',
            completed_internal_q=k,contribution_norm=float(np.linalg.norm(contribution)),
            elapsed_seconds=perf_counter()-started)),flush=True)
    old_result_path=Path('/var/tmp/p254-a43-material-cubic-kernel-002/result.json')
    old_result=json.loads(old_result_path.read_bytes())
    with np.load(authenticate(old_result['output']),allow_pickle=False) as saved_old_kernel:
        old_mean=saved_old_kernel['mean_force'][0,0]
    rows=np.asarray([0,1,2,93,94,95])
    errors['actual_saved_original_Gaussian_t0_selected_readout_bridge']=relative(physical_mean[0,rows],old_mean)
    mean=model.field_mean_force(means,covariances)
    tangent=model.field_mean_jvp(means,dmeans,dcovariances)
    mean_dual,covariance_dual=model.field_mean_transpose(means,cot)
    lhs=np.sum(cot*tangent);rhs=np.sum(mean_dual*dmeans)+np.sum(covariance_dual*dcovariances)
    errors['all64_internal_q_field_bilinear_transpose_identity']=relative(lhs,rhs)
    finite=(model.field_mean_force(means+h*dmeans,covariances+h*dcovariances)
        -model.field_mean_force(means-h*dmeans,covariances-h*dcovariances))/(2*h)
    errors['all64_internal_q_mean_central_direction']=relative(tangent,finite)
    scalar=model.field_mean_scalar(means,covariances)
    scalar_jvp=model.field_mean_scalar_jvp(means,covariances,dmeans,dcovariances)
    f=lambda t:model.field_mean_scalar(means+t*dmeans,covariances+t*dcovariances)
    scalar_finite=(f(-2*h)-8*f(-h)+8*f(h)-f(2*h))/(12*h)
    errors['common_Gaussian_cubic_scalar_finite_direction']=relative(scalar_jvp,scalar_finite)
    scalar_mu_dual,scalar_w_dual=model.field_mean_scalar_transpose(means,covariances,scalar_cot)
    errors['common_Gaussian_scalar_true_transpose_identity']=relative(scalar_cot*scalar_jvp,
        np.sum(scalar_mu_dual*dmeans)+np.sum(scalar_w_dual*dcovariances))
    errors['common_scalar_coordinate_gradient_is_negative_mean_force']=relative(scalar_mu_dual,-scalar_cot*mean[model.minus])
    mutations['field_hilbert_adjoint']=float(abs(lhs-(np.sum(mean_dual.conj()*dmeans)+np.sum(covariance_dual.conj()*dcovariances)))/max(1e-30,abs(lhs)))
    if max(errors.values())>=2e-12 or any(value<=1e-3 for value in mutations.values()):
        atomic_json(output/'failed-boundary-findings.json',dict(errors=errors,mutations=mutations,source=model.source_binding))
        raise AssertionError(dict(errors=errors,mutations=mutations))
    outputs=dict(mean_force=mean,directional_action=tangent,mean_transpose=mean_dual,covariance_transpose=covariance_dual,
        scalar=np.asarray(scalar),scalar_direction=np.asarray(scalar_jvp),point_force=point_force,
        point_direction=point_jvp,original_Gaussian_t0_mean_force=physical_mean)
    cache.put('actual-actions',outputs,metadata=dict(source=model.source_binding,probe_inputs=cache.manifest['items']['probe-inputs']['archive'],
        all_coordinates=111,all_internal_covariance_q=64,cubic_scalar_and_force_share_same_source=True,
        covariance_class='explicit translation-invariant ordered W[k,-k]; general covariance(k,l) route API also exposed',
        no_mass_scaling_no_equilibrium_no_rate=True))
    saved_outputs=cache.get('actual-actions')
    if any(not np.isfinite(value).all() for value in saved_outputs.values()):
        raise ValueError('genuine saved generalized action output nonfinite')
    cache.manifest['status']='complete';cache.publish()
    result=dict(schema='p254-a47-actual-source-complete-cubic-feedback-consumer-v1',status='complete',
        source=model.source_binding,source_bindings=specification['source_bindings'],
        probe_inputs=cache.manifest['items']['probe-inputs']['archive'],saved_actions=cache.manifest['items']['actual-actions']['archive'],
        completed_cache=file_binding(cache.path),errors=errors,algebraic_and_finite_tolerance=2e-12,mutations=mutations,
        original_Gaussian_input=dict(coefficient_manifest=model.source_record['coefficient_manifest'],mechanical_archive=mechanical_binding,
            selected_saved_member='ionic_ordered[:111,:111]',all64_original_covariance_blocks_share_original_product_preparation=True),
        original_A43_selected_readout_bridge=dict(result=file_binding(old_result_path),arrays=old_result['output'],rows=rows.tolist(),time_inverse_ry=0),
        actual_saved_quantities={name:dict(shape=list(value.shape),norm=float(np.linalg.norm(value)),max_abs=float(np.max(np.abs(value)))) for name,value in saved_outputs.items()},
        actual_original_Gaussian_t0_all111_mean_norm=float(np.linalg.norm(physical_mean[0])),
        compact_pair_cache=dict(entries=len(model.ionic._reciprocal_cache),array_bytes=sum(value.nbytes for pair in model.ionic._reciprocal_cache.values() for value in pair),maximum_entries=512),
        raw_probe_input_bytes=raw_input_bytes,raw_action_output_bytes=sum(value.nbytes for value in outputs.values()),
        temporary_reference_phi_cubes=1,persisted_phi_cubes=0,dense_allq_111_cubed_bank=False,external_time_kernel_bank=False,
        all128_electronic_sectors=True,all37_electronic_atoms=True,all111_coordinate_inputs=True,all64_internal_covariance_q=True,
        active_finite_mean_q=[0,1,3],active_finite_direction_q=[0,1,2,3],
        scope='Actual compact common Gaussian cubic scalar/mean-force directional map and true complex-bilinear transpose. Arbitrary111-coordinate legs; all64 internal covariance q. Electronic source state unchanged.',
        units=dict(scalar='Ry',mean='bohr',covariance='bohr^2',mean_force='Ry/bohr',
            Phi='Ry/bohr^3',time='inverse Ry; this boundary is instantaneous',
            transpose='complex-bilinear coordinate dual; never a Hilbert adjoint'),
        old_physical_oracles_replayed=False,whole_physical_goal_complete=False,lifetime_or_rate_computed=False,
        elapsed_seconds=perf_counter()-started)
    atomic_json(output/'feedback-result.json',result)
    print(json.dumps(result),flush=True)

if __name__=='__main__':
    main()
