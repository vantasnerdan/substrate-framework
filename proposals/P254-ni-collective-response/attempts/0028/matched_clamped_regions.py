"""Independent copies of one LOCAL ionic+field ClampedArray system.

Each region is a whole original supercell. There is no inter-region ionic
coupling AND no shared field. N fields at volume Vcell replace one field at
volume N*Vcell. This changes field multiplicity, vacuum and zero-point energy;
it is an explicit alternate CONTROL, not the original global-field experiment.
The tuples named q below are independent real-space blocks, not Fourier modes.
"""
import numpy as np


class IndependentRegions:
    def __init__(self, local, shape):
        if local.n != 1:
            raise ValueError('one local ionic+field block required')
        self.local = local
        self.shape = tuple(shape)
        self.n = int(np.prod(shape))
        self.indices = np.indices(shape).reshape(3,-1).T
        self.d, self.f = local.d, local.f
        self.masses, self.omega = local.masses, local.omega

    def cell(self, cell):
        raw = np.asarray(cell)
        if (raw.shape != (3,) or not np.issubdtype(raw.dtype,np.integer)
                or np.any(raw < 0) or np.any(raw >= self.shape)):
            raise ValueError('original region cell index outside grid')
        return raw

    def kicks(self, cell, coordinate):
        target = np.ravel_multi_index(tuple(self.cell(cell)), self.shape)
        single, _ = self.local.kicks((0,0,0), coordinate)
        local = tuple(single[0].copy() if k == target else np.zeros_like(single[0])
                      for k in range(self.n))
        common = tuple(single[0]/np.sqrt(self.n) for _ in range(self.n))
        return local, common

    def evolve(self, time_inverse_ry):
        return self.local.evolve(time_inverse_ry)

    def mean(self, evolution, initial=None):
        values = [self.local.mean(evolution, None if initial is None else (initial[k],))
                  for k in range(self.n)]
        return dict(ions=np.asarray([v['ions'][0,0,0] for v in values]).reshape(*self.shape,2*self.d),
            field=np.asarray([v['field'] for v in values]),
            q=tuple(v['q'][0] for v in values))

    def ordered_pair(self, left, right, cell, other, coordinates):
        same = np.array_equal(self.cell(cell), self.cell(other))
        local = self.local.ordered_pair(left,right,(0,0,0),(0,0,0),coordinates)
        return local if same else np.zeros_like(local)

    def field_ordered_pair(self, left, right):
        """One representative local field; cross-region connected entries zero."""
        return self.local.field_ordered_pair(left,right)

    def ion_field_ordered_pair(self, left, right, coordinates):
        return self.local.ion_field_ordered_pair(left,right,coordinates)

    def diagnostics(self, evolution, mean_q):
        zero = np.zeros_like(mean_q[0])
        result = self.local.diagnostics(evolution,(zero,))
        h, linear = self.local.h[0], self.local.linear
        result['mean_energy_ry'] = sum(value@h@value/2+linear@value for value in mean_q)
        result['covariance_energy_ry'] *= self.n
        result['energy_ry'] = result['mean_energy_ry']+result['covariance_energy_ry']
        result['mean_q_reality_relative_defect'] = max(
            float(np.linalg.norm(value-value.conj())/max(1.,np.linalg.norm(value)))
            for value in mean_q)
        result['normalization'] = 'sum over independent local ionic+field regions, including N field vacua'
        return result

    def coefficient_diagnostics(self):
        return dict(local=self.local.coefficient_diagnostics(), independent_region_count=self.n,
                    independent_field_count=self.n, shared_global_field=False,
                    inter_region_ionic_couplings=False, repair_applied=False)
