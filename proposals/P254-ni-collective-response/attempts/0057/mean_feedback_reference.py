"""Zero-reintegration consumer for actual common mean/curvature feedback bytes.

Gamma S and U accept exact saved nodes only; no transition interpolation.
Other-q S is the original explicit raw expm bath, not an updated all-q solution.
"""
from mean_feedback_functional import HERE,binding,authenticate
import json
from pathlib import Path
import numpy as np
from general_moving_reference import GaussianTransitions
from material_electron_ion_inputs import FrozenReference

class MeanFeedbackReference(GaussianTransitions):
    label='actual-common-mean-curvature-feedback-saved-U-GammaS-original-otherq-bath'
    def __init__(self,c,result_path,coefficient_binding):
        spec=json.loads(Path(result_path).read_bytes())
        if spec['schema']!='p254-common-mean-feedback-reference-v1' or spec['status']!='complete' or spec['coefficient_manifest']!=coefficient_binding:
            raise ValueError('complete matching mean-feedback material result required')
        for item in [spec['input_dataset'],spec['baseline_catalog'],*spec['physical_history_bindings'],*spec['source_bindings']]:authenticate(item)
        self._archive=np.load(authenticate(spec['archive']),allow_pickle=False)
        times=self._archive['times'];self.u=self._archive['U'];self.mean=self._archive['mean_gamma'];self.gamma_s=self._archive['S_gamma']
        nx=c.mechanical_sizes[0];nt=len(times)
        if (times.tolist()!=spec['times_inverse_ry'] or self.u.shape!=(nt,c.ns,c.m,c.m) or self.mean.shape!=(nt,nx) or self.gamma_s.shape!=(nt,nx,nx)
                or not all(np.isfinite(v).all() for v in (times,self.u,self.mean,self.gamma_s))):
            raise ValueError('full finite phase mean/U/S at every saved node required')
        if (not np.array_equal(self.u[0],np.broadcast_to(np.eye(c.m),(c.ns,c.m,c.m)))
                or not np.array_equal(self.gamma_s[0],np.eye(nx)) or not np.array_equal(self.mean[0],np.zeros(nx))
                or not np.array_equal(self._archive['rho0'],c.rho0) or not np.array_equal(self._archive['W0_gamma'],c.mechanical.w0[0])):
            raise ValueError('original prepared mean/rho/covariance and literal U0/S0 required')
        self.spec=spec;self.u_integration_calls=0;self.gamma_integration_calls=0
        self.bindings=[binding(result_path),spec['archive'],spec['input_dataset'],spec['baseline_catalog'],*spec['physical_history_bindings'],*spec['source_bindings']]
        super().__init__(c,times)
    def S(self,q,t):
        if int(q)!=0:return super().S(q,t)
        matches=np.flatnonzero(self.times==float(t))
        if len(matches)!=1:raise ValueError('Gamma time absent from saved bank; run producer, never interpolate S')
        return self.gamma_s[int(matches[0])]
    electron=FrozenReference.electron
    def close(self):
        self.cache.clear();self.cache_size=0;self._transitions=None;self._archive.close()
