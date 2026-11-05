"""Original source0/Gamma1038 physical density/force/current same-source join.

Consumes ONE explicitly reconstructed durable original-model baseline. Reuses
paid A37 basis/state/spectrum/force matrices and A38 current matrices. NEW work:
backward dual, full-PW density actions, and current-generation joint spectral
convolutions/CAR. No H/F/J/basis/eigen/DFT/preparation reconstruction. All retained
pieces are subtracted BEFORE same-state contraction, never added variances.
"""
from pathlib import Path
import argparse
import json
import resource
import sys
from time import perf_counter
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
sys.path[:0]=[str(ROOT/'src')]+[str(HERE.parent/n) for n in ('0041','0038','0037','0031','0030','0029','0028')]
import numpy as np
from canonical_force_filter import chebyshev_lobatto_transform
from canonical_density_filter import build_canonical_density_filter,load_canonical_density_filter
from durable_active_baseline import SavedMovingBaseline
from density_projection import BorrowedBathFrame,FourierDensityProjector
from density_integration import integrate_density_residual
from density_joint_readout import compose_same_source,substitution_readout
from consume_joint_influence import authenticate_handoff,read_bound
from current_projection_inputs import restore_inputs,load_projection_bank
from export_frozen_bath_comparator import car_moments,projected_factors
from expanded_influence_readout import save_json,save_npz
from integrate_saved_spectral import load_saved_spectral
from material_electron_ion_inputs import load_material
from reduced_electronic_correlations import ActiveBandCorrelations
from reduced_projection_cache import ProjectionCache, atomic_json
from substrate_framework.streamed_material_state import authenticate,file_binding

CONTRACT=dict(path=str(HERE/'implementation-contract-001.json'),bytes=7438,
    sha256='6cb22f829d8fd3450d0dba2606a57f50695b5382d213d4d06757182312695755')
SCOPE=('Original source0/Gamma expanded1038, Ni1/x t30 inverseRy. Same-state '
    'original moving-active output plus frozen-H full-minus-retained two-mode '
    'ionic, uniform-current and physical-density residual forcing, all cross '
    'orders and preparation retained. NOT selfconsistent A32 hybrid covariance, '
    'moving fullPW mean, allsource/allq result, collision/loss/heat/rate or whole E.')


def joint_coefficients(density_filter,directions,amplitude):
    """Put three physical kinds on common knots/degrees and separate probe legs.

    Column0 ionic,1 field,2 density. This permits ONE shared gap-moment pass
    while retaining each contribution and all ordered CAR cross pairs.
    """
    force=density_filter.force_filter
    nodes,transform=chebyshev_lobatto_transform(14)
    density=density_filter.density_chebyshev
    out=np.zeros((len(density),15,5+density.shape[2],3),complex)
    out[:,:,5:,2]=density[:,:,:,0]
    for segment,(left,right) in enumerate(density_filter.density_segment_bounds):
        samples=[]
        for j,node in enumerate(nodes):
            time=right if j==0 else left if j==14 else left+(right-left)*(node+1)/2
            ionic,field=force.evaluate(time)
            samples.append(np.r_[(directions.conj().T@ionic)[:,0],
                                  (amplitude.T@field)[:,0]])
        coeff=transform@np.asarray(samples)
        out[segment,:,:2,0]=coeff[:,:2]
        out[segment,:,2:5,1]=coeff[:,2:]
    return out


def consume(baseline_receipt,current_projection_complete,output,*,max_bytes=256*1024**2,
            original_source_snapshots=(),density_vertex_cache=None,
            vertex_cache_max_bytes=1024**3,density_vertex_donor=None,
            projection_equivalence=None,filter_receipt=None):
    started=perf_counter()
    contract=json.loads(authenticate(CONTRACT).read_bytes())
    producer,basis,temporal,amplitude,bindings=authenticate_handoff(
        current_projection_complete,original_source_snapshots=original_source_snapshots)
    baseline_binding=file_binding(baseline_receipt)
    baseline_meta=read_bound(baseline_binding)
    original=read_bound(baseline_meta['source']['original_result'])
    specification=original['configuration']['coefficient_specification']
    manifest_binding=contract['settled_material']['manifest']
    c,manifest,actual_manifest_binding=load_material(authenticate(manifest_binding))
    if actual_manifest_binding!=manifest_binding:raise ValueError('settled material binding differs')
    if c.m!=12 or c.d!=111 or c.n!=64:raise ValueError('original loaded12/27 layout required')
    baseline=SavedMovingBaseline(baseline_receipt,c)
    if baseline_meta['coefficient_specification']!=specification:
        raise ValueError('baseline original coefficient specification differs')
    old_filter=read_bound(bindings['original_filter_packet'])
    if baseline_meta['source']['original_result']!=old_filter['source']['original_result']:
        raise ValueError('baseline and original A30 observable source differ')
    settings=baseline_meta['settings']
    corr=ActiveBandCorrelations(c,baseline,rtol=settings['rtol'],atol=settings['atol'])
    reused_filter=file_binding(filter_receipt) if filter_receipt else None
    filter_packet=read_bound(reused_filter) if reused_filter else None
    if filter_packet:
        if (filter_packet['baseline']!=baseline_binding
                or filter_packet['filter']['q_index']!=0
                or filter_packet['filter']['force_filter']['coordinates']!=[0]
                or filter_packet['filter']['force_filter']['final_time_inverse_ry']!=30.):
            raise ValueError('saved filter is not the same original baseline/output')
        density_filter=load_canonical_density_filter(filter_packet)
    else:
        density_filter=build_canonical_density_filter(corr,0,30.,(0,))
    directory=Path(output).resolve();directory.mkdir(parents=True,exist_ok=False)
    filter_bytes=sum(a.nbytes for a in density_filter.arrays.values())
    output_matrix_bytes=6*1038**2*16
    baseline_bytes=sum(int(b['bytes']) for b in baseline_meta['arrays'].values())
    if baseline_bytes+filter_bytes+output_matrix_bytes>max_bytes:
        raise ValueError('complete baseline/filter/joint payload exceeds shared256MiB cap; nothing truncated')
    filter_binding=(filter_packet['archive'] if filter_packet else
                    save_npz(directory/'canonical-density-filter.npz',density_filter.arrays))
    save_json(directory/'density-filter-complete.json',dict(status='complete',
        baseline=baseline_binding,archive=filter_binding,filter=density_filter.metadata,
        original_completed_filter=reused_filter))
    progress=lambda event:print('[physics] '+json.dumps(event,allow_nan=False),flush=True)
    reference=read_bound(producer['current_reference'])
    inputs=restore_inputs(basis,reference,progress)
    bank=load_projection_bank(producer['basis'],basis,inputs)
    with np.load(authenticate(producer['basis']),allow_pickle=False) as archive:
        q=archive['projected_active_coordinates'];projector=archive['active_projector']
    with np.load(authenticate(producer['current_bank']),allow_pickle=False) as archive:
        currents=[archive['current_'+axis] for axis in 'xyz']
    saved=load_saved_spectral(producer['a37_filtered_operators'])
    frame=BorrowedBathFrame(inputs.context.psi,bank['inverse_sqrt_gram'],bank['extra_basis'])
    density_projector=FourierDensityProjector(
        inputs.context._coupling,frame,frame,ncells=64,block_size=128)
    block=c.q(0)
    # Source record0 is one explicit spin; other-spin vertices are exact zeros,
    # not absent channels or an implicit factor2.
    spin=specification['indexing']['slot_spin'][0]
    physical=[int(i) for i in np.flatnonzero(np.asarray(block.feature_spin)==spin)]
    ng=len(block.feature_g)
    vertex_directory=(Path(density_vertex_cache).resolve() if density_vertex_cache
                      else directory/'physical-density-vertices')
    if vertex_directory==ROOT or ROOT in vertex_directory.parents:
        raise ValueError('physical density vertex cache must remain outside Git')
    physical_spec=dict(kind='full-PW-physical-density-vertices',basis=producer['basis'],
        physical_reference=producer['current_reference'],
        original_pw_basis_identity=producer['original_pw_basis_identity'],
        ncells=c.n,explicit_spin=spin,feature_g=block.feature_g.tolist())
    current_projection=[file_binding(HERE/'density_projection.py')]
    donor_binding=file_binding(density_vertex_donor) if density_vertex_donor else None
    equivalence_binding=file_binding(projection_equivalence) if projection_equivalence else None
    vertices=ProjectionCache(vertex_directory,dict(**physical_spec,
        current_projection_implementation=current_projection,
        donor_manifest=donor_binding,projection_equivalence=equivalence_binding))
    if donor_binding:
        donor=read_bound(donor_binding)
        original_projection=donor['specification']['projection_implementation']
        if {k:v for k,v in donor['specification'].items()
                if k!='projection_implementation'}!=physical_spec:
            raise ValueError('donor density matrices have different physical inputs')
        if equivalence_binding is None:
            raise ValueError('changed-representation donor requires actual oracle admission')
        equivalence=read_bound(equivalence_binding)
        if (equivalence['passed'] is not True
                or equivalence['new_projection']!=current_projection[0]
                or equivalence['original_projection']!=original_projection[0]):
            raise ValueError('donor/current projection implementation authority differs')
        for key,row in donor['items'].items():
            if not vertices.has(key):
                authenticate(row['archive'])
                vertices.manifest['items'][key]=dict(archive=row['archive'],
                    metadata=dict(**row['metadata'],producer_manifest=donor_binding,
                        producer_implementation=original_projection))
        vertices.publish()
    calls=[]
    def provider(index):
        if index<2:
            name='mode'+str(index)
            return bank['force_'+name]/8,bank['adjoint_force_'+name]/8
        if index<5:
            # Existing qualified physical Hermitian current legs, no new J action.
            return currents[index-2],currents[index-2]
        feature=index-5;start=perf_counter()
        record=dict(feature=feature,spin=spin,columns=1038,column_block=128,status='started',
            proposed_forward_actions=(1038+127)//128,proposed_adjoint_actions=(1038+127)//128,
            representation='exact finite physical FFT slot shift; not new field/cutoff')
        calls.append(record)
        atomic_json(directory/'density-actions.json',dict(calls=calls))
        try:
            key='density/%03d'%feature
            if vertices.has(key):
                record.update(status='reused',physical_forward_actions=0,
                              physical_adjoint_actions=0)
            else:
                value=density_projector.project(block.feature_g[feature%ng])
                existing=sum(int(row['archive']['bytes'])
                             for row in vertices.manifest['items'].values())
                if existing+value[0].nbytes+value[1].nbytes+4096>vertex_cache_max_bytes:
                    raise ValueError('new physical vertex cache exceeds its authorized cap')
                vertices.put(key,dict(density=value[0],adjoint_density=value[1]),
                    metadata=dict(feature=feature,spin=spin,ncells=c.n,
                                  producer_implementation=current_projection))
                del value
                record.update(status='complete',physical_forward_actions=(1038+127)//128,
                              physical_adjoint_actions=(1038+127)//128)
            # The spectral consumer uses the actual authenticated persisted pair.
            saved_vertex=vertices.get(key)
            value=saved_vertex['density'],saved_vertex['adjoint_density']
            record.update(archive=vertices.manifest['items'][key]['archive'],
                active_density_vertex_defect_fro=float(np.linalg.norm(
                    q.conj().T@value[0]@q-block.density[0,feature].T/8)))
            if len(vertices.manifest['items'])==len(physical):
                vertices.manifest['status']='complete'
                vertices.publish()
            progress(dict(density_feature=feature,status=record['status'],archive=record['archive']))
            return value
        except BaseException as error:
            record.update(status='failed',error=repr(error))
            raise
        finally:
            record['seconds']=perf_counter()-start
            atomic_json(directory/'density-actions.json',dict(calls=calls))
    coefficients=joint_coefficients(density_filter,inputs.directions,amplitude)
    convolution=ProjectionCache(vertex_directory/'spectral-convolution-v1',dict(
        kind='full-density-spectral-row-convolution',row_block_size=128,
        baseline=baseline_binding,filter_archive=filter_binding,
        physical_vertex_manifest=file_binding(vertices.path),physical_reference=producer['current_reference'],
        basis=producer['basis'],current_bank=producer['current_bank'],
        spectral_archive=producer['a37_filtered_operators'],
        implementation=[file_binding(HERE/name) for name in
            ('density_integration.py','density_moments.py','consume_density_influence.py')]))
    integration=integrate_density_residual(saved,provider,[*range(5),*(5+i for i in physical)],
        projector,density_filter.density_segment_bounds,coefficients,progress=progress,checkpoint=convolution)
    initial=density_filter.force_filter.initial_dual
    _,zs,_=corr._split(0,initial)
    prep=q@zs[0,:,:,0].T@q.conj().T
    # Original active compression must be matched, never assumed from labels.
    rho_full=bank['projected_frame']@bank['density_coefficients']@bank['projected_frame'].conj().T
    compression_defect=float(np.linalg.norm(q.conj().T@rho_full@q-c.rho0[0]))
    matrices,adjoints=compose_same_source(prep,prep.conj().T,
        integration['operators'][0],integration['adjoint_operators'][0],
        integration['operators'][1],integration['adjoint_operators'][1],
        integration['operators'][2],integration['adjoint_operators'][2])
    moments=car_moments(bank['projected_frame'],bank['gram'],bank['density_coefficients'],
        projected_factors(matrices,adjoints,bank['projected_frame'],basis['basis_identity'],
            producer['original_state_authority'],'bohr'),producer['original_factor_identity'],
        basis['basis_identity'],producer['original_state_authority'])
    original_ordered=corr.initial_ordered(0,initial,initial)[0,0]
    from reduced_correlation_transpose import wick_edge_pair
    represented=wick_edge_pair(zs[:1],zs[:1],c.rho0[:1],c.rho0[:1])[0,0]
    substitution=substitution_readout(moments,original_ordered,represented)
    operators=save_npz(directory/'density-joint-operators.npz',dict(
        density_residual=matrices['density'],density_residual_adjoint=adjoints['density'],
        matched=matrices['matched'],matched_adjoint=adjoints['matched']))
    result=dict(schema='P254-A42-source0-density-joint-v1',status='complete',scope=SCOPE,
        baseline=baseline_binding,baseline_generation=baseline_meta['generation'],
        filter=filter_binding,filter_metadata=density_filter.metadata,operators=operators,
        inputs=bindings,source=producer['original_source'],basis=producer['basis'],
        spectral_reference=producer['a37_filtered_operators'],current_bank=producer['current_bank'],
        original_state_authority=producer['original_state_authority'],
        original_factor_identity=producer['original_factor_identity'],
        original_pw_basis_identity=producer['original_pw_basis_identity'],
        moments=moments,substitution=substitution,source_active_compression_defect=compression_defect,
        normalization=dict(ncells=64,explicit_spin=spin,density_profile='exp(2pi i G.r)/sqrt(Omega)',
            density_vertex='physical fullPW O_G/sqrt(N)',density_filter='K.T P.T C_rho.T z/sqrt(N)',
            ionic='paid F/sqrt(N)',field='paid J, coefficient sqrt(N*w)*M.T*zPi',
            extra_spin_factor=1,extra_source_weight=1,spatial_inverse_fourier_applied=False),
        work=dict(total_seconds=perf_counter()-started,physical_density_calls=calls,
            integration=integration['diagnostics'],H_F_J_replays=0,basis_eigensolves=0,
            mean_replays_in_consumer=0,backward_solves=0 if reused_filter else 1,
            reused_completed_filter=reused_filter,
            joint_convolutions='NEW current-generation convolutions of paid physical matrices; old integrated operators not silently reused',
            baseline_filter_operator_payload_bytes=baseline_bytes+filter_bytes+output_matrix_bytes,
            physical_density_vertex_cache=file_binding(vertices.path),
            vertex_cache_max_bytes=vertex_cache_max_bytes,
            spectral_convolution_cache=file_binding(convolution.path),
            recovery='committed physical vertices and completed spectral row blocks reusable; only current block and eigenvertex RAM work uncommitted',
            peak_process_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024),
        implementation=[file_binding(p) for p in [
            *(HERE/name for name in ('consume_density_influence.py','canonical_density_filter.py',
                'density_projection.py','density_integration.py','density_moments.py','density_joint_readout.py',
                'durable_active_baseline.py')),
            HERE.parent/'0041/material_electron_ion_inputs.py',
            HERE.parent/'0030/canonical_force_filter.py',
            HERE.parent/'0038/consume_joint_influence.py']])
    return save_json(directory/'result.json',result)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-receipt',required=True)
    parser.add_argument('--current-projection-complete',required=True)
    parser.add_argument('--output',required=True)
    parser.add_argument('--max-bytes',type=int,default=256*1024**2)
    parser.add_argument('--original-source-snapshot',action='append',default=[])
    parser.add_argument('--density-vertex-cache',
                        help='reuse only exact-input independently projected physical density matrices')
    parser.add_argument('--vertex-cache-max-bytes',type=int,default=1024**3)
    parser.add_argument('--density-vertex-donor',help='original immutable vertex manifest; no bank copy')
    parser.add_argument('--projection-equivalence',help='actual changed-representation oracle admission')
    parser.add_argument('--filter-receipt',help='reuse the actual same-baseline/output completed dual filter')
    args=parser.parse_args()
    print(json.dumps(consume(args.baseline_receipt,args.current_projection_complete,args.output,
        max_bytes=args.max_bytes,original_source_snapshots=args.original_source_snapshot,
        density_vertex_cache=args.density_vertex_cache,
        vertex_cache_max_bytes=args.vertex_cache_max_bytes,
        density_vertex_donor=args.density_vertex_donor,
        projection_equivalence=args.projection_equivalence,
        filter_receipt=args.filter_receipt),allow_nan=False))
