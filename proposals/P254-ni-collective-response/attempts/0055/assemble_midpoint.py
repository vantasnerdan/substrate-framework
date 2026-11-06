"""Merge disjoint new-pair EE sums, bind complete catalog, consume real bytes."""
from refinement_common import *
from refine_ee import KEYS

def _load(path):return json.loads(Path(path).read_bytes())

def merge_ee():
    prior=HERE/'material/ee-consumer.json'
    if prior.exists():
        complete=_load(prior)
        if complete['status']=='complete':
            authenticate(complete['archive'])
            return complete
    records=[_load(HERE/'material'/name/'scope.json') for name in ('ee-a-001','ee-b-001','ee-c-001')]
    owned=[q for r in records for q in r['completed_outer_q']]
    if sorted(owned)!=list(range(64)) or len(set(owned))!=64 or any(r['status']!='complete' for r in records):
        return dict(status='partial',completed_outer_q=owned,required_outer_q=list(range(64)),shards=records)
    sums={key:np.zeros((128,5,12,12),complex) for key in KEYS}
    for record in records:
        with np.load(authenticate(record['archive']),allow_pickle=False) as z:
            for key in KEYS:sums[key]+=z[key]
    old_binding=catalog()['archives']['ee'];path=HERE/'material/ee.npz'
    with np.load(authenticate(old_binding),allow_pickle=False) as old:
        arrays={key:old_merge(old[key],sums[key],True) for key in KEYS}
    np.savez(path,**arrays)
    with np.load(path,allow_pickle=False) as saved,np.load(old_binding['path'],allow_pickle=False) as old:
        unchanged=all(saved[key][:,a,b].tobytes()==old[key][:,oa,ob].tobytes() for key in KEYS for a,oa in OLD_MAP.items() for b,ob in OLD_MAP.items())
        vector=np.linspace(.3,1.1,12)*(1+.17j)
        actions={key:float(np.linalg.norm(saved[key][:,1,0]@vector)) for key in KEYS}
        finite=all(np.isfinite(saved[key]).all() for key in KEYS)
    if not unchanged or not finite:raise AssertionError('EE final merge failed literal old/finite consumer')
    receipt=dict(status='complete',archive=binding(path),completed_outer_q=owned,shards=[r['archive'] for r in records],literal_oldfour_byte_identical=unchanged,actual_midpoint_pair_actions=actions,new_pair_count=5,old_pairs_recomputed=0,external_slots=128,inner_h=64,bands=12)
    save(HERE/'material/ee-consumer.json',receipt);return receipt

def run():
    ee=merge_ee();old=catalog();pieces={'ee':ee}
    for kind in ('ei','joint','represented'):
        p=HERE/'material'/kind/'consumer.json'
        pieces[kind]=_load(p) if p.exists() else dict(status='missing')
    p=HERE/'material/cubic-complete/result.json'
    if not p.exists():p=HERE/'material/cubic/result.json'
    pieces['cubic']=_load(p) if p.exists() else dict(status='missing')
    complete=all(r.get('status')=='complete' for r in pieces.values())
    old_equal={}
    for kind,piece in pieces.items():
        if piece.get('status')!='complete':continue
        with np.load(authenticate(piece['archive']),allow_pickle=False) as new,np.load(authenticate(old['archives'][kind]),allow_pickle=False) as source:
            names=KEYS if kind=='ee' else (('greater','lesser','retarded') if kind in ('ei','cubic') else tuple(k for k in new.files if k in source.files))
            same=all((new[k][:,a,b] if kind=='ee' else new[k][a,b]).tobytes()==(source[k][:,oa,ob] if kind=='ee' else source[k][oa,ob]).tobytes() for k in names for a,oa in OLD_MAP.items() for b,ob in OLD_MAP.items())
            if not same:raise AssertionError('literal original four changed: '+kind)
            old_equal[kind]=dict(byte_identical=True,array_keys=list(names))
    if pieces['ei'].get('status')=='complete':
        with np.load(pieces['ei']['archive']['path'],allow_pickle=False) as new,np.load(authenticate(old['archives']['joint']),allow_pickle=False) as source:
            same=all(new['pi_'+channel][a,b].tobytes()==source['ionic_pi_'+channel][oa,ob].tobytes() for channel in ('greater','lesser','retarded') for a,oa in OLD_MAP.items() for b,ob in OLD_MAP.items())
            if not same:raise AssertionError('literal ionic Gamma Pi four changed')
            old_equal['ei_Gamma_ionic_Pi']=dict(byte_identical=True,original_archive=old['archives']['joint'])
    combined={}
    if all(pieces[k].get('status')=='complete' for k in ('ei','joint','represented','cubic')):
        with np.load(pieces['ei']['archive']['path'],allow_pickle=False) as ei,np.load(pieces['joint']['archive']['path'],allow_pickle=False) as joint,np.load(pieces['represented']['archive']['path'],allow_pickle=False) as represented,np.load(pieces['cubic']['archive']['path'],allow_pickle=False) as cubic:
            for channel in ('greater','lesser','retarded'):
                sig=ei[channel]+joint['joint_sigma_'+channel]-joint['ionic_sigma_'+channel]
                pi=joint['joint_pi_'+channel]-represented['pi_'+channel]
                pi[:,:,:111,:111]+=cubic[channel]
                x=np.linspace(.1,.9,114)*(1+.37j);v=np.linspace(.3,1.1,12)*(1+.17j)
                if not np.isfinite(sig).all() or not np.isfinite(pi).all():raise FloatingPointError('nonfinite combined actual component')
                combined[channel]=dict(electron_Sigma_midpoint_action_norm=float(np.linalg.norm(sig[1,0]@v)),Gamma_Pi_midpoint_action_norm=float(np.linalg.norm(pi[1,0]@x)),Gamma_ordered_noise_or_response_norm=float(np.linalg.norm(pi)),directfield_ionicGamma_replacement_count=1,represented_Pi_subtraction_count=1,cubic_ionic111_insertion_count=1)
    runs=[_load(p) for p in sorted((HERE/'invocations').glob('*/exit.json'))]
    new_total=sum(r['elapsed_seconds'] for r in runs)
    new_ee=sum(r['elapsed_seconds'] for r in runs if r['mode']=='ee')
    new_cubic=sum(r['elapsed_seconds'] for r in runs if r['mode'].startswith('cubic') or r['mode']=='merge')
    grant_path=HERE/'direct-cubic-continuation-grant.json'
    cubic_cap=_load(grant_path)['amended_combined_cubic_cap_seconds'] if grant_path.exists() else 9000
    ledgers=dict(all_new_native_seconds=new_total,all_new_native_cap_seconds=6000,all_new_ee_native_seconds=new_ee,A46_original_material_seconds=6633.830,A46_material_cumulative_seconds=6633.830+new_ee,A46_material_cap_seconds=24000,A46_WaveI_separate_seconds=6.722,cubic_new_native_seconds=new_cubic,cubic_combined_original_seconds=8611.287,cubic_combined_cumulative_seconds=8611.287+new_cubic,cubic_combined_cap_seconds=cubic_cap,cubic_reader_accounting='all component merge/consumer invocation elapsed conservatively included')
    payload=sum(p.stat().st_size for p in HERE.rglob('*') if p.is_file())
    if payload>256*1024**2 or new_total>6000 or 8611.287+new_cubic>cubic_cap:raise RuntimeError('actual allocation exceeded')
    physical_sources={
        '0046':('all_slot_coulomb_shared.py','all_slot_coulomb_blas.py','all_slot_coulomb_contractions.py'),
        '0044':('bare_coulomb_factors.py','coulomb_material_inputs.py'),
        '0041':('electron_ion_functional.py','material_electron_ion_inputs.py'),
        '0050':('joint_gamma_adapter.py','represented_electron_channel.py','joint_material_inputs.py'),
        '0051':('cubic_contour.py',),
        '0047':('cubic_feedback.py',),
        '0052':('general_moving_reference.py',),
        '0045':('moving_reference_propagation.py',)}
    imported_sources=[binding(HERE.parent/attempt/name) for attempt,names in physical_sources.items() for name in names]
    result=dict(schema='p254-a55-genuine-three-node-handoff-v1',status='complete' if complete else 'partial',times_inverse_ry=TIMES.tolist(),pieces=pieces,literal_oldfour_verification=old_equal,actual_saved_numeric_consumer=combined,all_new_invocations=runs,ledgers=ledgers,new_owned_payload_bytes=payload,reference=_load(HERE/'reference/result.json'),source_bindings=[binding(p) for p in sorted(HERE.glob('*.py'))],original_catalog=binding(HERE.parent/'0053/actual-input-catalog-001.json'),findings=['Genuine midpoint U integrated from saved Hamiltonian; literal original U0/U1 retained','Only five new ordered pairs evaluated; complete q sums use original physical backends','Old functional proofs carried unchanged, not replayed'],limits=['Fixed components only, not selfconsistent coupled evolution','Midpoint alone does not license T30/global temporal accuracy/physical calibration/rates/fullE','No symmetrization/clipping/tolerance change or missing piece zero fill'])
    result['physical_backend_source_bindings']=imported_sources
    if complete:
        new_catalog=dict(old,schema='p254-actual-fixed-kernel-component-catalog-v1',component_times_inverse_ry=TIMES.tolist(),selected_times_inverse_ry=TIMES.tolist(),archives={k:r['archive'] for k,r in pieces.items()},reference_bundle=binding(HERE/'reference/bundle.json'),refinement_identity=dict(case='P254-original-full128-all64-12band-Gamma114',new_pairs=NEW_PAIRS,old_pairs=OLD_MAP,original_catalog=result['original_catalog']))
        save(HERE/'actual-input-catalog-001.json',new_catalog);result['catalog']=binding(HERE/'actual-input-catalog-001.json')
        sys.path.insert(0,str(HERE.parent/'0053'))
        from material_step_inputs import load_step
        c,ref,consumed_catalog,sigma,pi,represented=load_step(HERE/'actual-input-catalog-001.json')
        try:
            v=np.linspace(.3,1.1,12)*(1+.17j)
            x=np.linspace(.1,.9,114)*(1+.37j)
            numeric=dict(status='passed',reader_source=binding(HERE.parent/'0053/material_step_inputs.py'),catalog=result['catalog'],times_inverse_ry=ref.times.tolist(),fermion_actions={k:float(np.linalg.norm(value[1,0]@v)) for k,value in sigma.items()},Gamma_Pi_actions={k:float(np.linalg.norm(value[1,0]@x)) for k,value in pi.items()},Gamma_represented_actions={k:float(np.linalg.norm(value[1,0]@x)) for k,value in represented.items()},fermion_shapes={k:list(value.shape) for k,value in sigma.items()},Gamma_shapes={k:list(value.shape) for k,value in pi.items()},uniform_history=True,canonical_raw_component_assembly=True,physical_or_evolved_success_claim=False)
            if not all(np.isfinite(value).all() for value in [*sigma.values(),*pi.values(),*represented.values()]):raise FloatingPointError('canonical actual consumer nonfinite')
            save(HERE/'actual-canonical-saved-component-consumer-001.json',numeric)
            result['canonical_saved_component_consumer']=binding(HERE/'actual-canonical-saved-component-consumer-001.json')
        finally:ref.close()
    save(HERE/'actual-genuine-midpoint-handoff-001.json',result);return result
