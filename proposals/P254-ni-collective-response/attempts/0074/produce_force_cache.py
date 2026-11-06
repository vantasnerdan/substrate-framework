"""Produce the genuine31-node/all64 bare covariance force, without a mean solve."""
from pathlib import Path
import json
import time
import numpy as np
from fullwindow_force_cache import (HERE, SCHEMA, SOURCE, CompactCubic,
    FullWindowForceCache, bare_nodes, gamma_force, binding, authenticate, SavedMean)
from history_storage import atomic_json
from paid_mean_feedback_cache import PaidMeanFeedbackCache
from general_moving_reference import GaussianTransitions
from material_electron_ion_inputs import load_material


def produce(output):
    started = time.monotonic()
    out = Path(output)
    out.mkdir(exist_ok=False)
    frozen = json.loads((HERE/'preregistration.json').read_bytes())
    for name in ('genuine_reference','source_descriptor','input_old_cache','existing_force_cache'):
        authenticate(frozen[name])
    c,manifest,coefficient_binding = load_material()
    source = CompactCubic.from_source(SOURCE)
    old_path = authenticate(frozen['input_old_cache'])
    old_spec = json.loads(old_path.read_bytes())
    old_result_path = authenticate(old_spec['original_result'])
    old_result = json.loads(old_result_path.read_bytes())
    old_catalog_path = authenticate(old_result['baseline_catalog'])
    old_catalog = json.loads(old_catalog_path.read_bytes())
    old_reference_path = authenticate(old_catalog['reference_bundle'])
    old_reference = json.loads(old_reference_path.read_bytes())
    new_reference = json.loads(authenticate(frozen['genuine_reference']).read_bytes())
    # The old three-node all64 force contains dressed Gamma. Only its unchanged
    # Gaussian nonGamma bath at the exact shared nodes0,1 can be carried forward.
    if (old_reference['coefficient_manifest'] != coefficient_binding or
            new_reference['coefficient_manifest'] != coefficient_binding or
            old_reference['physical_history_bindings'][:4] != new_reference['physical_history_bindings'] or
            old_spec['compact_source'] != source.source_binding or
            old_spec['compact_source_arrays'] != source.source_record['arrays']):
        raise ValueError('source-compatible old nonGamma/Phi000 reuse not established')
    inputs = [binding(HERE/'preregistration.json'),frozen['source_descriptor'],
        frozen['genuine_reference'],binding(old_path),binding(old_catalog_path),
        binding(old_reference_path),old_result['baseline_catalog'],
        *old_result['physical_history_bindings'],*old_result['source_bindings']]
    for item in inputs:
        authenticate(item)
    old = PaidMeanFeedbackCache(old_path,coefficient_binding)
    try:
        with SavedMean(frozen['source_descriptor']['path']) as bare:
            times = bare.times.copy()
            if (not bare.is_bare or bare.reference_identity != frozen['genuine_reference'] or
                    not np.array_equal(times,np.arange(31,dtype=float)) or
                    not np.array_equal(bare.rho0,c.rho0) or
                    not np.array_equal(bare.w0,c.mechanical.w0[0]) or
                    int(bare.archive['n']) != c.n or
                    (c.n,c.nq,c.ns,c.m,c.d,c.f) != (64,64,128,12,111,3)):
                raise ValueError('actual original all31 bare source/preparation differs')
            rho,w = bare_nodes(bare)
            inputs += bare.bindings+old.bindings
        gamma = gamma_force(source,w)
        force = np.zeros((31,111),complex)
        shared = []
        for a,t in enumerate(times):
            found = np.flatnonzero(old.times==t)
            if len(found):
                b = int(found[0])
                force[a] = old.full_force[b]+(gamma[a]-old.old_gamma_force[b])
                shared.append(dict(new_index=a,old_index=b,time=float(t),
                    bare_gamma_norm=float(np.linalg.norm(gamma[a])),
                    old_dressed_gamma_norm=float(np.linalg.norm(old.old_gamma_force[b])),
                    gamma_delta_norm=float(np.linalg.norm(gamma[a]-old.old_gamma_force[b]))))
        if [x['new_index'] for x in shared] != [0,1]:
            raise ValueError('unexpected exact old/new node intersection')
        phi = np.asarray(old.phi000).copy()
    finally:
        old.close()
    pending = list(range(2,31))
    force[pending] = gamma[pending]
    bath = GaussianTransitions(c,times,cache_bytes=64*1024**2)
    qcost = []
    (out/'boundaries').mkdir()
    # q0 is the full114 raw source's true ionic block, not a field-cubic term.
    qcost.append(dict(internal_q=0,new_node_contractions=31,seconds=0.,
        source='genuine saved228-phase source -> raw114 greater -> ionic111'))
    for k in range(1,64):
        step = time.monotonic()
        qs = source.indices[[0,k,int(source.minus[k])]]
        for a in pending:
            # Original A57 force contraction and independently evaluated q/-q
            # raw Gaussian legs. No conjugate/opposite-q reconstruction.
            greater,_ = bath.boson(k,a,a)
            force[a] -= .5*source.contract(qs,covariance=1j*greater)
        entry = dict(internal_q=k,new_node_contractions=len(pending),
            seconds=time.monotonic()-step,expm_calls=bath.expm_calls)
        qcost.append(entry)
        boundary = out/'boundaries'/('through-q%02d.npz'%k)
        np.savez(boundary,times=times,partial_force=force,bare_gamma_force=gamma,
            completed_internal_q=np.arange(k+1),complete_nodes=np.array([0,1]),
            newly_computed_nodes=np.array(pending))
        atomic_json(out/'checkpoint.json',dict(boundary=binding(boundary),
            completed_internal_q=list(range(k+1)),shared_node_reuse=shared,
            complete=k==63,timings=qcost,elapsed_seconds=time.monotonic()-started,
            recovery='Complete numeric q boundaries retained; automatic/warm continuation unsupported.'))
        print(json.dumps(dict(stage='A74_BARE_FORCE_Q',**entry)),flush=True)
    archive = out/'bare-covariance-force.npz'
    np.savez(archive,times=times,bare_full_force=force,bare_gamma_force=gamma,
        Phi000=phi,node_rho=rho,node_ordered_gamma_covariance=w)
    implementation = [binding(HERE/n) for n in
        ('fullwindow_force_cache.py','produce_force_cache.py','native_force_cache.py','launch_force_cache.py')]
    implementation += [binding(HERE.parent/d/n) for d,n in
        (('0057','mean_feedback_functional.py'),('0057','paid_mean_feedback_cache.py'),
         ('0057','produce_mean_feedback.py'),('0047','cubic_feedback.py'),
         ('0052','general_moving_reference.py'),('0063','saved_mean.py'),
         ('0063','history_storage.py'),('0063','history_provider.py'),
         ('0041','material_electron_ion_inputs.py'))]
    spec = dict(schema=SCHEMA,status='complete',archive=binding(archive),
        coefficient_manifest=coefficient_binding,original_input_dataset=frozen['source_descriptor'],
        compact_source=source.source_binding,compact_source_arrays=source.source_record['arrays'],
        original_result=binding(out/'checkpoint.json'),input_bindings=inputs,
        implementation_bindings=implementation,completed_original_internal_q=list(range(64)),
        times_inverse_ry=times.tolist(),shared_node_reuse=shared,
        reuse=dict(Phi000=old_spec['archive'],old_cache=binding(old_path),
            unchanged_nonGamma_nodes=[0.,1.],new_all64_nodes=list(range(2,31)),
            arithmetic='Fbare_shared=Fpaid_full+(Fbare_Gamma-Fpaid_dressed_Gamma); new nodes accumulated q0..q63. Future Fnew=Fbare_full+(Fnew_Gamma-Fbare_Gamma). Ordinary binary64 delta arithmetic.'),
        scalar='Original A57 CommonMeanFunctional/CompactCubic: force=-sumq Phi[0,q,-q]:Wq/2; Phi000 retained from authenticated paid source; no field cubic.',
        transitions=dict(cache_bytes=bath.cache_bytes,expm_calls=bath.expm_calls,
            cache_hits=bath.cache_hits,all64_phase_history_allocated=False),
        recovery='Per-q complete numeric boundaries saved. Automatic/warm resume and interrupted source restoration unsupported; partial boundaries are not a complete cache.',
        elapsed_seconds=time.monotonic()-started,
        limits=['Prescribed raw bare covariance force only; not a common mean or intrinsic trajectory.',
            'No interpolation, refill, projection, stationarity, conservation, temporal-accuracy, ranking, rate or heat claim.'])
    atomic_json(out/'result.json',spec)
    # The one affected smoke consumes the actual persisted all31 source/cache,
    # computes rho/W/force through the public API, and uses stored Phi in the
    # owning CommonMeanFunctional scalar, with no electronic/mean evolution.
    from mean_feedback_functional import CommonMeanFunctional
    with FullWindowForceCache(out/'result.json',coefficient_binding) as cache:
        actual = cache.forcing_for(frozen['source_descriptor']['path'],source)
        model = CommonMeanFunctional(c,source,phi=cache.phi000)
        probe = np.linspace(.01,.03,model.nx).astype(complex)*(1+.17j)
        consumer = dict(result=binding(out/'result.json'),archive=binding(archive),
            input_dataset=actual['input_dataset'],times=actual['times'].tolist(),
            force_byte_identical=actual['cubic_covariance_force'].tobytes()==force.tobytes(),
            gamma_byte_identical=actual['gamma_force'].tobytes()==gamma.tobytes(),
            rho_byte_identical=actual['rho'].tobytes()==rho.tobytes(),
            covariance_byte_identical=actual['ordered_gamma_covariance'].tobytes()==w.tobytes(),
            force_norm=[float(np.linalg.norm(f)) for f in actual['cubic_covariance_force']],
            gamma_norm=[float(np.linalg.norm(f)) for f in actual['gamma_force']],
            Phi000_norm=float(np.linalg.norm(cache.phi000)),
            cubic_scalar_probe=[dict(real=float(model.cubic_scalar(probe,f).real),
                imag=float(model.cubic_scalar(probe,f).imag)) for f in actual['cubic_covariance_force']],
            nonGamma_contractions=actual['nonGamma_contractions'],
            Gamma_contractions=actual['Gamma_contractions'],native_history_solves=0,
            mean_evolution_calls=0,completed_internal_q=list(range(64)))
        if not all(consumer[k] for k in ('force_byte_identical','gamma_byte_identical',
                'rho_byte_identical','covariance_byte_identical')):
            raise AssertionError('actual saved all31 public cache consumer differs')
        atomic_json(out/'saved-consumer.json',consumer)
    bath.cache.clear()
    return dict(result=binding(out/'result.json'),archive=binding(archive),
        consumer=binding(out/'saved-consumer.json'),elapsed_seconds=time.monotonic()-started)
