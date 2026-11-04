"""Initial-noise origins of an already pulled-back moving Gamma observable.

This is a decomposition of ONE original product preparation through ONE saved
response, not independent noise added to that response or a dissipation budget.
The source0 operator uses the actual persisted general eigenvector frame. All
ionic/field and retained/omitted cross orders survive without a repair or weight.
"""
import numpy as np

from export_frozen_bath_comparator import car_moments, projected_factors
from reduced_correlation_transpose import wick_edge_pair


def _covariance_lookup(moments):
    return {(entry['left'], entry['right']): entry['value']
            for entry in moments['covariance']}


def split_operator_moments(operators, bank, *, factor_identity, basis_identity,
                           state_authority, units):
    """Reuse the existing finite-bath original-frame/D CAR consumer.

    P=Qa Qa† acts in the declared finite bath, not the original full-PW space.
    No projected state or active-population replacement is made. Physical
    adjoints are formed AFTER the observable's general similarity transform.
    """
    q = bank['projected_active_coordinates']
    matrices = {}
    for name, operator in operators.items():
        retained = q @ (q.conj().T @ operator @ q) @ q.conj().T
        matrices['full:'+name] = operator
        matrices['retained:'+name] = retained
        matrices['omitted:'+name] = operator-retained
    adjoints = {name: matrix.conj().T for name, matrix in matrices.items()}
    factors = projected_factors(matrices, adjoints, bank['projected_frame'],
                                basis_identity, state_authority, units)
    return car_moments(bank['projected_frame'], bank['gram'],
        bank['density_coefficients'], factors, factor_identity, basis_identity,
        state_authority)


def decompose_initial_noise(duals, mechanical_w0, ionic_dimension,
        active_densities, bath_bank, coordinate_bank, *, source_slot,
        factor_identity, basis_identity, state_authority):
    """Return all represented source covariances in saved output-axis order.

    Layout: original canonical mechanical, ascending active slots except the
    promoted slot, full bath density. Density matrices use C-order and the
    bilinear dual, so its physical one-body observable is V Z.T V^-1.
    """
    duals = np.asarray(duals)
    w0 = np.asarray(mechanical_w0)
    rho = np.asarray(active_densities)
    frame = bath_bank['projected_frame']
    if duals.ndim == 1:
        duals = duals[:, None]
    if (duals.ndim != 2 or w0.ndim != 2 or w0.shape[0] != w0.shape[1]
            or rho.ndim != 3 or rho.shape[1] != rho.shape[2]):
        raise ValueError('saved dual, canonical covariance and active density layout required')
    nx, ns, m, b = len(w0), len(rho), rho.shape[1], len(frame)
    if (not 0 <= source_slot < ns or not 0 <= 2*ionic_dimension <= nx
            or duals.shape[0] != nx+(ns-1)*m*m+b*b):
        raise ValueError('exact saved moving-hybrid source/mechanical layout required')
    for value in (duals, w0, rho):
        if not np.isfinite(value).all():
            raise ValueError('finite unmodified preparation and duals required')
    p = duals.shape[1]
    x = duals[:nx]
    split = 2*ionic_dimension
    regions = {'ionic': slice(0, split), 'field': slice(split, nx)}
    mechanical = {a+'_'+bname: x[ia].T @ w0[ia, ib] @ x[ib]
                  for a, ia in regions.items() for bname, ib in regions.items()}
    other_slots = [slot for slot in range(ns) if slot != source_slot]
    end = nx+(ns-1)*m*m
    other = duals[nx:end].reshape(ns-1, m, m, p, order='C')
    electronic = np.empty((ns, p, p), complex)
    for index, slot in enumerate(other_slots):
        electronic[slot] = wick_edge_pair(other[index:index+1],
            other[index:index+1], rho[slot:slot+1], rho[slot:slot+1])
    v = coordinate_bank['right_eigenvectors']
    inverse = coordinate_bank['inverse_right_eigenvectors']
    if v.shape != (b, b) or inverse.shape != (b, b):
        raise ValueError('actual persisted source0 coordinate bank required')
    coefficients = np.moveaxis(duals[end:].reshape(b, b, p, order='C'), -1, 0)
    physical = v @ coefficients.swapaxes(-1, -2) @ inverse
    moments = split_operator_moments(
        {str(index): operator for index, operator in enumerate(physical)}, bath_bank,
        factor_identity=factor_identity, basis_identity=basis_identity,
        state_authority=state_authority, units='bohr')
    covariance = _covariance_lookup(moments)
    def block(left, right):
        return np.asarray([[covariance[left+':'+str(i), right+':'+str(j)]
                            for j in range(p)] for i in range(p)])
    electronic[source_slot] = block('full', 'full')
    source0 = {a+'_'+bname: block(a, bname)
               for a in ('retained', 'omitted') for bname in ('retained', 'omitted')}
    mechanical_total = sum(mechanical.values())
    electronic_total = electronic.sum(axis=0)
    return dict(mechanical=mechanical, electronic_by_slot=electronic,
        electronic_total=electronic_total, source0=source0,
        source0_full=electronic[source_slot],
        source0_reconstruction_defect=electronic[source_slot]-sum(source0.values()),
        mechanical_total=mechanical_total,
        total=mechanical_total+electronic_total,
        source0_means=moments['means'], source0_admissibility=moments['admissibility'],
        scope='initial product-noise origins through the SAME saved coupled response; '
              'finite bath source0 plus other active sectors; not full-PW uncertainty, '
              'independent added noise, physical damping, heat or a lifetime')
