"""Actual physical-factor diagnostics for initial matched tensor reuse.

The nonlocal zero-field operator bound does NOT bound its ionic/field jets.
Every unmeasured coefficient defect remains explicitly unbounded. Reported
inequalities are floating evaluations, not interval-enclosed certificates.
"""
import numpy as np

from .material_angular_intertwiners import original_channel_intertwiner, intertwiner_diagnostics
from .material_physical_operators import build_physical_factors
from .material_reference_transport import slot_transport, raw_state_transport
from .qe_real_bloch_grid import RealBlochBasis, real_basis


def authenticated_basis(control, slot, source, admission, expected_identity):
    basis = real_basis(control.grid, slot, source, admission)
    if not isinstance(basis, RealBlochBasis) or basis.identity != expected_identity:
        raise ValueError('original RealBlochBasis ledger identity differs')
    return basis


def initial_state_diagnostics(state, slot):
    target = raw_state_transport(state, antiunitary=slot.orbital_conjugation)
    gram = state.gram.conj() if slot.orbital_conjugation else state.gram
    d = state.density_coefficients.conj() if slot.orbital_conjugation else state.density_coefficients
    if not np.array_equal(target.gram, gram) or not np.array_equal(target.density_coefficients, d):
        raise ValueError('orbital conjugation must transform BOTH Gram and D')
    number = np.einsum('ij,ji->', state.gram, state.density_coefficients)
    return dict(source_identity=state.source_identity,
        gram_conjugated=slot.orbital_conjugation, density_coefficients_conjugated=slot.orbital_conjugation,
        spin_unchanged=slot.ispin, positive_time_symmetry=False,
        original_number_real=float(number.real), original_number_imag=float(number.imag))


def nonlocal_factor_bound(source_beta, target_beta, source_d, target_d, q):
    """Bt=CQ+E; norm defect <= ||E||(2||C||+||E||)||Dt||+||C||²||QDtQ*-Ds||.

    Frobenius norms majorize projector spectral norms. Ds may already be
    conjugated by the caller. Active UPF D is never changed or made real.
    """
    residual = target_beta-source_beta @ q
    c, e = float(np.linalg.norm(source_beta)), float(np.linalg.norm(residual))
    defect = q @ target_d @ q.conj().T-source_d
    dnorm, qnorm = float(np.linalg.norm(target_d, 2)), float(np.linalg.norm(defect, 2))
    return dict(projector_residual_frobenius=e,
        d_conjugacy_spectral=qnorm,
        nonlocal_operator_bound_ry=e*(2*c+e)*dnorm+c*c*qnorm)


def physical_diagnostics(control, inputs, slot, source, basis, source_factors):
    """Full actual PW row banks with declared Q, atom images and Bloch phases.

    For each target atom, source irt(a), n=atom_image_offsets[:,a]:
    C = exp[-epsilon*2pi*i*((G+ks).f + ks.n)] conjugate_if_needed(Bs).
    This includes nonsymmorphic translation and atomic image phases. Q(R)
    handles proper/improper angular channels; epsilon is already in K Bs.
    """
    transform = slot_transport(control.grid, slot)
    target = (source_factors if source_factors.sources[-1]['identity'] == basis.identity
              else build_physical_factors(inputs, basis))
    epsilon = slot.epsilon
    phase = np.exp(-epsilon*2j*np.pi*((source.miller.T+slot.source_k_fractional)
                                    @ slot.operation.fractional_translation))
    records, total = [], 0.0
    sources = {s.name: s for s in source_factors.species}
    for species in target.species:
        original = sources[species.name]
        q = original_channel_intertwiner(inputs.tables[species.name].upf,
                                          slot.operation.rotation)
        ds = original.active_d.conj() if slot.orbital_conjugation else original.active_d
        nh = ds.shape[0]
        positions = {atom-1: i for i, atom in enumerate(original.atom_indices)}
        angular = intertwiner_diagnostics(q, ds, species.active_d)
        for j, atom in enumerate(species.atom_indices):
            a = atom-1
            i = positions[int(transform.atom_pullback[a])]
            bs = original.beta[:, i*nh:(i+1)*nh]
            bs = bs.conj() if slot.orbital_conjugation else bs
            image_phase = np.exp(-epsilon*2j*np.pi*np.dot(
                slot.source_k_fractional, slot.atom_image_offsets[:, a]))
            c = (phase*image_phase)[:, None]*bs
            diagnostic = nonlocal_factor_bound(c, species.beta[:, j*nh:(j+1)*nh],
                                                ds, species.active_d, q)
            total += diagnostic['nonlocal_operator_bound_ry']
            records.append(dict(species=species.name, original_atom=int(transform.atom_pullback[a])+1,
                                target_atom=atom, **diagnostic))
        # One angular diagnostic per species; per-atom records retain all defects.
        records.append(dict(species=species.name, angular=angular))
    return dict(basis_identity=basis.identity, physical_factor_provenance=dict(target.sources[-1]),
        geometry=transform.diagnostics(), phase_modulus_defect=float(np.max(abs(abs(phase)-1))),
        cartesian_vs_declared_rotation_max=float(np.max(abs(
            transform.cartesian_rotation-slot.operation.rotation))),
        nonlocal_zero_field_operator_bound_ry=total, atoms_and_channels=records,
        qualification='floating analytic nonlocal zero-field bound; no interval enclosure',
        excluded_from_bound=['local ionic gradient and curvature', 'nonlocal ionic gradient and curvature',
                             'field gradient/current/contact', 'mixed ionic-field coefficient'],
        no_tolerance_license=True, positive_time_symmetry=False)
