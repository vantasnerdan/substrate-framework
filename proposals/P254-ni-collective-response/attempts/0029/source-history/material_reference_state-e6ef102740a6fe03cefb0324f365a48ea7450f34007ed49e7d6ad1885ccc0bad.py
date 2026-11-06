"""Full-column Bloch-source preparation and common-reference moments.

The preparation is the existing explicitly chosen symmetric-frame quasifree
family. It is not a thermal state of the reconstructed H0. Source-orbit reuse
requires identical populations and the admitted unitary row transport; it
never replaces spin/full-grid measures by an IBZ multiplicity implicitly.
"""
from dataclasses import dataclass
from types import MappingProxyType

import numpy as np
from scipy.linalg.blas import zgemm

from .physical_factor_preparation import (
    FactorFunctional, FactorPreparation, coefficients_from_orthonormal_populations,
)


@dataclass(frozen=True)
class SourceState:
    gram: np.ndarray
    populations: np.ndarray
    density_coefficients: np.ndarray
    admissibility: object
    source_identity: str


def prepare_source(coefficients, populations, *, source_identity, basis_identity, authority):
    """Preserve all raw rows/columns; compute the source Gram only once."""
    psi = np.asarray(coefficients)
    if (psi.ndim != 2 or not 0 < psi.shape[1] <= psi.shape[0]
            or not np.isfinite(psi).all() or not source_identity or not basis_identity):
        raise ValueError("finite independent raw frame and exact source/basis identities required")
    gram = zgemm(1., psi, psi, trans_a=2)
    d = coefficients_from_orthonormal_populations(gram, populations, mapping_authority=authority)
    gram.flags.writeable = False
    # Identity's M and T are both S, even in a nonorthogonal raw frame.
    metadata = MappingProxyType({"number": MappingProxyType(dict(
        operator_identity="full-PW-identity", adjoint_identity="full-PW-identity",
        units="electrons", pw_basis_identity=basis_identity, source_authority=authority))})
    functional = FactorFunctional(gram, MappingProxyType({"number": gram}),
        MappingProxyType({("number", "number"): gram}), metadata,
        source_identity, basis_identity, psi.shape[0])
    preparation = FactorPreparation(functional, d, authority)
    populations = np.array(populations, dtype=float, copy=True)
    populations.flags.writeable = False
    return SourceState(gram, populations, preparation.density_coefficients,
                       MappingProxyType(preparation.admissibility), source_identity)


def restore_source(arrays, *, source_identity):
    """Restore authenticated completed state arrays, without refilling.

    The producer first authenticates the original receipt, source/frame
    metadata, and exact current populations. Old diagnostics are retained.
    """
    if set(arrays) != {"gram", "populations", "density_coefficients", "density_eigenvalues"}:
        raise ValueError("complete authenticated source-state arrays required")
    n = len(arrays["populations"])
    if (not source_identity or arrays["gram"].shape != (n, n)
            or arrays["density_coefficients"].shape != (n, n)
            or arrays["density_eigenvalues"].shape != (n,)):
        raise ValueError("restored source-state dimensions or identity differ")
    for value in arrays.values():
        value.flags.writeable = False
    return SourceState(arrays["gram"], arrays["populations"], arrays["density_coefficients"],
                       MappingProxyType({"density_eigenvalues": arrays["density_eigenvalues"]}),
                       source_identity)


def source_population_indices(grid, ledger, spin, source_ik):
    """Return all explicit slots of one source; reject a different orbit state."""
    indices = np.flatnonzero((ledger["spin"] == spin) & (ledger["source_ik"] == source_ik))
    slots = {(s.ispin, s.mesh_index): s for s in grid.slots}
    if not len(indices):
        raise ValueError("requested source has no populated full-grid slot")
    population = ledger["populations"][indices[0]]
    selected = []
    for index in indices:
        key = (int(ledger["spin"][index]), tuple(int(x) for x in ledger["mesh_index"][index]))
        slot = slots[key]
        if (slot.source_ik != source_ik
                or slot.source_binding["sha256"] != str(ledger["source_sha256"][index])
                or slot.weight != ledger["measures"][index]
                or not np.array_equal(population, ledger["populations"][index])):
            raise ValueError("source/orbit population or measure differs; no cache equivalence")
        selected.append((int(index), slot))
    return population, tuple(selected)


def hamiltonian_state_moments(psi, state, hpsi):
    """Unweighted full-row H0 moments and actual commutator norm, without PW².

    For S=Psi†Psi, M=Psi†H Psi and T=(H Psi)†H Psi,
    ||[H,rho]||_F² = 2 Tr(D S D T) - 2 Tr(D M D M).
    Raw complex roundoff is retained; no clipping or stationarity is imposed.
    """
    d, s = state.density_coefficients, state.gram
    m = zgemm(1., psi, hpsi, trans_a=2)
    t = zgemm(1., hpsi, hpsi, trans_a=2)
    ds, dm = d @ s, d @ m
    norm_squared = 2 * (np.einsum("ij,ji->", ds @ d, t) - np.einsum("ij,ji->", dm, dm))
    projected_drift = -1j * (m @ ds - s @ dm)
    return dict(projected_h0_ry=m, h0_action_gram_ry2=t,
        one_body_energy_ry=np.asarray(np.einsum("ij,ji->", d, m)),
        density_commutator_norm_squared_ry2=np.asarray(norm_squared),
        projected_density_drift_ry=projected_drift)


def sector_density_moments(frame, state, coupling, *, check_direction=None):
    """Unweighted density/positive gradient, independently checkpointable."""
    psi = frame.coefficients
    if (frame.basis.slot.source_binding["sha256"] != state.source_identity
            or psi.shape[1] != len(state.populations)
            or coupling.npw != psi.shape[0]):
        raise ValueError("raw frame, prepared source and checked coupling differ")
    density, gradient = coupling.prepared_moments(psi, state.density_coefficients)
    arrays = dict(density_bohr3=density, positive_gradient_ry_bohr=gradient,
        number_pair=np.asarray([coupling.grid_measure_bohr3*np.sum(density),
                                np.einsum("ij,ji->", state.density_coefficients, state.gram)]))
    if check_direction is not None:
        action = coupling.ionic_action(check_direction, psi)
        projected = zgemm(1., psi, action, trans_a=2)
        arrays["gradient_action_pair"] = np.asarray([
            check_direction @ gradient, np.einsum("ij,ji->", state.density_coefficients, projected)])
    return arrays


def sector_hamiltonian_moments(frame, state, hamiltonian):
    """Full-row H0 moments/drift, separate from the density/force reduction."""
    if (frame.basis.slot.source_binding["sha256"] != state.source_identity
            or hamiltonian.ispin != frame.ispin
            or hamiltonian.factors.sources[-1].get("identity") != frame.basis.identity):
        raise ValueError("prepared source, spin and full-row Hamiltonian basis differ")
    psi = frame.coefficients
    hpsi = hamiltonian.apply(psi)
    arrays = hamiltonian_state_moments(psi, state, hpsi)
    arrays["saved_label_residual_norm_ry"] = np.linalg.norm(
        hpsi-psi*frame.eigenvalues_ry[None, :], axis=0)
    return arrays
