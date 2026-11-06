"""Initial-state source tensors on the ORIGINAL matched slot ledger.

For epsilon=-1 both raw Gram and density coefficient D are conjugated by
qe_real_bloch_grid. This is within-spin orbital conjugation, not a symmetry
of positive-time trajectories. If F_s(u,a) is the fixed-state scalar, then
F_t(u,a) = conjugate_if_needed F_s(R.T u[atom_inverse], epsilon R.T a)
at the transformed geometry. Differentiating gives the rules below. No
Hermitian/real projection, geometry repair, refilling or multiplicity occurs.

Rounded saved geometry is unchanged. Except density and number, transported
values are declared geometric predictions, NOT exact original-target values.
Direct-target comparisons and operator diagnostics live in the certificate
module; neither supplies a tolerance license for unmeasured coefficients.
"""
import numpy as np

from .material_matched_reference import KEYS, slot_identity
from .material_reference_transport import slot_transport

GEOMETRY_POLICY = 'declared-geometric-initial-coefficient-transport'
NORMALIZATION = 'unweighted-one-explicit-spin-slot'
ODD_FIELD_KEYS = ('field_gradient', 'current_hartree', 'mixed_ionic_field')


def identity_representative(grid, spin, source_ik):
    """Original physical representative, allowing its admitted reciprocal fold.

    With R=I, f=0 and epsilon=1, m_rep=m_source-g and k_rep=k_source+g.
    Thus every physical momentum and raw coefficient is unchanged. The
    periodic frame's common exp(-2pi i g.x) cancels from all coefficients.
    A nonzero integer BZ fold is NOT a different physical representative.
    """
    candidates = [s for s in grid.slots if (s.ispin, s.source_ik) == (spin, source_ik)
        and not s.orbital_conjugation
        and np.array_equal(s.operation.rotation, np.eye(3))
        and not np.any(s.operation.fractional_translation)]
    if len(candidates) != 1:
        raise ValueError('one original identity representative required')
    return candidates[0]


def require_arrays(arrays, nat, fft_size):
    shapes = dict(density_bohr3=(fft_size,), positive_gradient_ry_bohr=(3*nat,),
        local_onsite_ry_bohr2=(nat, 3, 3), nonlocal_onsite_ry_bohr2=(nat, 3, 3),
        field_gradient=(3,), current_hartree=(3,), contact=(3, 3),
        mixed_ionic_field=(nat, 3, 3), number=())
    for name, shape in shapes.items():
        value = np.asarray(arrays[name])
        if value.shape != shape or not np.isfinite(value).all():
            raise ValueError('full finite original coefficient shape required: '+name)
    if 'number_pair' not in arrays or not np.isfinite(arrays['number_pair']).all():
        raise ValueError('original independent number_pair required')


def tensor_transport(arrays, transform, *, orbital_conjugation):
    """Polar ionic/field indices; improper rotations have NO determinant factor.

    One field derivative is epsilon-odd, two are even. Complex values are
    conjugated for epsilon=-1 before rotation, including density and number.
    Mixed blocks retain ionic-major/field-minor order on both matrix legs.
    """
    if type(orbital_conjugation) is not bool:
        raise ValueError('explicit within-spin orbital conjugation boolean required')
    require_arrays(arrays, len(transform.atom_pullback), len(transform.fft_pullback))
    values = {k: np.asarray(arrays[k]).conj() if orbital_conjugation else np.asarray(arrays[k])
              for k in (*KEYS, 'number_pair')}
    epsilon = -1 if orbital_conjugation else 1
    r = transform.cartesian_rotation
    if r.shape != (3, 3) or not np.allclose(r.T @ r, np.eye(3), atol=2e-12, rtol=0):
        raise ValueError('orthogonal physical Cartesian rotation required')
    result = dict(density_bohr3=transform.density(values['density_bohr3']),
        positive_gradient_ry_bohr=transform.gradient(values['positive_gradient_ry_bohr']),
        number=values['number'], number_pair=values['number_pair'])
    for name in ('local_onsite_ry_bohr2', 'nonlocal_onsite_ry_bohr2'):
        result[name] = transform.onsite(values[name])
    for name in ('field_gradient', 'current_hartree'):
        result[name] = epsilon * (r @ values[name])
    result['contact'] = r @ values['contact'] @ r.T
    result['mixed_ionic_field'] = epsilon * transform.onsite(values['mixed_ionic_field'])
    return result


def transport_slot(grid, slot, representative, arrays, *, geometry_policy):
    """All original slots/weights remain explicit; this call applies NO weight."""
    grid.require_slot(representative)
    if slot_identity(representative) != slot_identity(identity_representative(
            grid, slot.ispin, slot.source_ik)):
        raise ValueError('source coefficients must be the original identity representative')
    if geometry_policy != GEOMETRY_POLICY:
        raise ValueError('explicit declared geometric approximation required')
    transform = slot_transport(grid, slot)
    result = tensor_transport(arrays, transform, orbital_conjugation=slot.orbital_conjugation)
    metadata = dict(slot=slot_identity(slot), geometry_policy=GEOMETRY_POLICY,
        normalization=NORMALIZATION, coefficient_origin='transformed-original-representative',
        original_representative=slot_identity(representative),
        geometry=transform.diagnostics(), positive_time_symmetry=False,
        exact_kinematic_scope=['density FFT pullback with orbital conjugation', 'number'],
        unbounded_coefficient_defect_scope=[k for k in KEYS if k not in ('density_bohr3', 'number')],
        approximation='rounded atom positions and physical UPF data unchanged; no coefficient-error budget assigned')
    return result, metadata


def coefficient_residuals(predicted, direct):
    """Actual same-state/target differences, not whole-orbit upper bounds."""
    result = {}
    for name in (*KEYS, 'number_pair'):
        a, b = np.asarray(predicted[name]), np.asarray(direct[name])
        if a.shape != b.shape or not all(np.isfinite(v).all() for v in (a, b)):
            raise ValueError('finite identical direct/transported dimensions required')
        delta = a-b
        absolute = float(np.linalg.norm(delta.reshape(-1)))
        scale = float(np.linalg.norm(b.reshape(-1)))
        result[name] = dict(max_absolute=float(np.max(np.abs(delta))),
            frobenius=absolute, direct_frobenius=scale,
            relative_frobenius=absolute/scale if scale else None,
            bound_scope='measured specified initial state and target only; not an operator bound')
    return result
