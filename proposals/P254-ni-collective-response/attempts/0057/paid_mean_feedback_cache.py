"""Post-terminal export of paid all64 mean forcing and reusable Gamma vertex.

The all64 sum is not replayed. Store the paid full sum plus its old q0 term;
replace only that term for a later dressed dataset. Arithmetic cancellation is
ordinary binary64, not an assertion of bit-identical separate nonGamma sums.
"""
from pathlib import Path
import json,numpy as np
from mean_feedback_functional import SOURCE,CompactCubic,binding,authenticate,save
from cubic_contour import vertex

class PaidMeanFeedbackCache:
    def __init__(self,result_path,coefficient_binding):
        self.spec=json.loads(Path(result_path).read_bytes())
        if self.spec['schema']!='p254-paid-mean-feedback-cache-v1' or self.spec['coefficient_manifest']!=coefficient_binding:
            raise ValueError('paid cache/material binding mismatch')
        for entry in [self.spec['original_result'],self.spec['original_input_dataset'],self.spec['compact_source'],self.spec['compact_source_arrays'],*self.spec['implementation_bindings']]:authenticate(entry)
        self._archive=np.load(authenticate(self.spec['archive']),allow_pickle=False)
        self.times=self._archive['times'];self.phi000=self._archive['Phi000'];self.full_force=self._archive['paid_full_force'];self.old_gamma_force=self._archive['paid_gamma_force']
        if self.times.tolist()!=[0.,.5,1.] or self.phi000.shape!=(111,111,111) or self.full_force.shape!=(3,111) or self.old_gamma_force.shape!=(3,111):raise ValueError('complete paid force/Gamma vertex axes differ')
        self.bindings=[binding(result_path),self.spec['archive'],self.spec['original_result'],self.spec['original_input_dataset'],self.spec['compact_source'],self.spec['compact_source_arrays'],*self.spec['implementation_bindings']]

    def forcing_for(self,dataset,source):
        """Return actual raw node rho/W and full force, reusing all paid q!=0.

        Caller supplies the same CompactCubic source already needed for the mean
        scalar. No material bank, bare bath, or nonGamma contraction is rebuilt.
        """
        if source.source_binding!=self.spec['compact_source'] or source.source_record['arrays']!=self.spec['compact_source_arrays']:raise ValueError('cubic source changed: paid bath not reusable')
        with np.load(dataset,allow_pickle=False) as data:
            times=data['times']
            if not np.array_equal(times,self.times):raise ValueError('paid nonGamma bath nodes changed')
            rho=-1j*np.asarray([data['electron_less'][i,i] for i in range(3)])
            w=1j*np.asarray([data['gamma_less'][i,i]+data['gamma_ret'][i,i] for i in range(3)])
        if rho.shape!=(3,128,12,12) or w.shape!=(3,114,114):raise ValueError('actual full128/Gamma114 dressed nodes required')
        new_gamma=np.asarray([-.5*source.contract(source.indices[[0,0,0]],covariance=value[:111,:111]) for value in w])
        full=self.full_force+(new_gamma-self.old_gamma_force)
        return dict(times=times,rho=rho,ordered_gamma_covariance=w,cubic_covariance_force=full,gamma_force=new_gamma,input_dataset=binding(dataset),cache_bindings=self.bindings,nonGamma_contractions=0,Gamma_contractions=3)

    def close(self):self._archive.close()


def export_paid_mean_cache(result_path,output):
    result=json.loads(Path(result_path).read_bytes());out=Path(output);out.mkdir(parents=True,exist_ok=False)
    authenticate(result['archive']);authenticate(result['input_dataset'])
    source=CompactCubic.from_source(SOURCE)
    with np.load(result['archive']['path'],allow_pickle=False) as original:
        times=original['times'];full=original['node_cubic_covariance_force'];w=original['node_ordered_gamma_covariance']
    oldgamma=np.asarray([-.5*source.contract(source.indices[[0,0,0]],covariance=value[:111,:111]) for value in w])
    phi=vertex(source,source.indices[[0,0,0]])
    archive=out/'paid-mean-forcing.npz'
    np.savez(archive,times=times,paid_full_force=full,paid_gamma_force=oldgamma,Phi000=phi)
    spec=dict(schema='p254-paid-mean-feedback-cache-v1',archive=binding(archive),coefficient_manifest=result['coefficient_manifest'],original_result=binding(result_path),original_input_dataset=result['input_dataset'],compact_source=source.source_binding,compact_source_arrays=source.source_record['arrays'],implementation_bindings=[binding(__file__)],
        completed_original_internal_q=list(range(64)),reused_nonGamma_internal_q=list(range(1,64)),new_export_contractions_internal_q=[0],
        update='Fnew=Fpaid_full+(Fnew_Gamma-Fpaid_Gamma); all63 nonGamma contributions retained with original opposite-q raw bath, source, coefficients and exact nodes. No symmetry reduction or new zero modes.',
        arithmetic='Paid full force was accumulated originally in q0..q63 order. Delta replacement uses binary64 addition/subtraction, not a separately re-summed nonGamma bank. Identical Gamma inputs recover identical full-force bytes.',
        reuse_limits='Only unchanged original prescribed nonGamma bath, material coefficients, cubic source and [0,.5,1] nodes. Any change invalidates this cache. Phi000 is a stored original source vertex, not an interpolated physical kernel.')
    save(out/'result.json',spec)
    cache=PaidMeanFeedbackCache(out/'result.json',result['coefficient_manifest'])
    try:
        actual=cache.forcing_for(result['input_dataset']['path'],source)
        consumer=dict(original_material_input=result['input_dataset'],archive=binding(archive),result=binding(out/'result.json'),full_force_same_input_byte_identical=actual['cubic_covariance_force'].tobytes()==full.tobytes(),full_force_same_input_max_error=float(np.max(np.abs(actual['cubic_covariance_force']-full))),nonGamma_contractions=actual['nonGamma_contractions'],Gamma_contractions=actual['Gamma_contractions'],Phi000_norm=float(np.linalg.norm(cache.phi000)),covered_original_internal_q=list(range(64)))
        if not consumer['full_force_same_input_byte_identical']:raise AssertionError(consumer)
        save(out/'saved-consumer.json',consumer)
    finally:cache.close()
    return dict(result=binding(out/'result.json'),archive=binding(archive),consumer=binding(out/'saved-consumer.json'),findings=consumer)
