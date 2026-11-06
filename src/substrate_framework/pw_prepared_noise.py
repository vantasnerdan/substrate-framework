"""Prepare actual PW bilinears for an auxiliary independent-H0 quasi-free state.

The finite-BvK convention is A_q = Ncells**(-1/2) sum_k c†_{k+q} O c_k;
unit-cell PW matrix elements therefore carry weight 1/Ncells in two-point
functions, not in Bernoulli probabilities. See Giustino (2017), Sec. III.B.2,
Eqs. (37)-(38), and Appendix A (original manuscript SHA256
5d980627ad272b2a92d1f760f0e6a0883b1b3fcc941436903a361024ce330251).
This normalization citation does not identify our bare positive ionic gradients
with screened phonon-mode vertices or establish coupled thermal equilibrium.
"""
from __future__ import annotations

from itertools import product
from typing import Any

import numpy as np

from .pw_response_data import PWResponseModel
from .pw_transition_space import PWTransitionSpace, finite_array
from .qe_pw_response_frames import _coordinate_equal, _fermi


def _finite_bvk(model: PWResponseModel, reciprocal_basis, mesh):
    """Validate geometry and explicit spin identities, not just weight sums."""
    dimensions = np.asarray(mesh)
    if dimensions.shape != (3,) or dimensions.dtype.kind not in "iu" or np.any(dimensions < 1):
        raise ValueError("mesh must contain three positive integers")
    basis = finite_array(reciprocal_basis, (3, 3), "reciprocal basis")
    if basis.dtype.kind not in "fiu":
        raise ValueError("reciprocal basis must be real")
    try:
        q_fractional = np.linalg.solve(basis, model.q_cart_2pi_alat)
    except np.linalg.LinAlgError as error:
        raise ValueError("reciprocal basis must be invertible") from error
    q_index = q_fractional * dimensions
    if not _coordinate_equal(q_index, np.rint(q_index)):
        raise ValueError("q is incompatible with the finite BvK mesh")
    cells = int(np.prod(dimensions))
    records = model.metadata.get("source_channels", [])
    by_identity = {record["identity"]: record for record in records}
    identities = [channel.identity for channel in model.channels]
    if (len(by_identity) != len(records) or len(set(identities)) != len(identities)
            or set(by_identity) != set(identities)):
        raise ValueError("source_channels must match unique model channel identities exactly")
    expected = set(product(*(range(int(n)) for n in dimensions)))
    seen: list[set[tuple[int, ...]]] = [set(), set()]
    validated = []
    for channel in model.channels:
        record = by_identity[channel.identity]
        if str(record["native_spin"]) != str(channel.spin + 1):
            raise ValueError("source metadata and explicit channel spin disagree")
        k = finite_array(record["k_source_2pi_alat"], (3,), "source k")
        k_fractional = np.linalg.solve(basis, k)
        scaled = k_fractional * dimensions
        if not _coordinate_equal(scaled, np.rint(scaled)):
            raise ValueError("source k is outside the unshifted finite BvK mesh")
        index = tuple((np.rint(scaled).astype(np.int64) % dimensions).tolist())
        if index in seen[channel.spin]:
            raise ValueError("duplicated reciprocal-equivalent source k within explicit spin")
        seen[channel.spin].add(index)
        weight = float(record["quadrature_weight"])
        if not np.isfinite(weight) or not _coordinate_equal(np.array(weight), np.array(1 / cells)):
            raise ValueError("finite BvK source quadrature must be uniform 1/Ncells")
        validated.append({"identity": channel.identity, "spin": channel.spin,
                          "native_ik": record["native_ik"],
                          "k_source_2pi_alat": np.asarray(k).tolist(),
                          "k_fractional": k_fractional.tolist(), "mesh_index": list(index),
                          "quadrature_weight": weight})
    if any(points != expected for points in seen):
        raise ValueError("each explicit spin requires the complete finite BvK source mesh")
    return dimensions, cells, q_fractional, validated


def prepare_pw_noise(model: PWResponseModel, reciprocal_basis, mesh=(2, 2, 2)):
    """Return plus/minus complete-target spectral recipes and preparation metadata.

    Every recipe contains identity and the ordered_transition_correlation keyword
    arguments except times. Source model arrays must remain immutable during use;
    this adapter neither changes their flags nor modifies their contents. Existing
    representably-zero FD source omissions are inherited; no further cutoff is
    introduced. Both target spectra and both actual gradient actions are retained,
    with their target plane-wave rows projected into the target energy eigenbasis.
    """
    space = PWTransitionSpace(model)
    dimensions, cells, q_fractional, records = _finite_bvk(model, reciprocal_basis, mesh)
    mu = float(model.metadata["common_mu_Ry"])
    width = float(model.metadata["actual_degauss_Ry"])
    if not np.isfinite(mu) or not np.isfinite(width) or width <= 0:
        raise ValueError("common_mu must be finite and actual_degauss finite and positive")
    branches: dict[str, list[dict[str, Any]]] = {"plus": [], "minus": []}
    for channel, maps, record in zip(model.channels, space.maps, records):
        # Recompute from energies: dividing weighted subnormal tails loses data.
        raw_source = _fermi((mu - channel.source_energies) / width)
        for branch, coordinates in zip(("plus", "minus"), maps):
            energies = getattr(channel, f"{branch}_energies")
            vectors = getattr(channel, f"{branch}_vectors")
            adjoint = vectors.conj().T
            operators = np.zeros((8, energies.size, raw_source.size), dtype=np.complex128)
            # Channel g stores PW action coefficients, not target spectral rows.
            operators[:6] = adjoint @ getattr(channel, f"g_{branch}")
            # Unit-cell average of u_target* u_source. No omega or spin factor.
            # Only source-supported grids are FFT'd; no full target grid exists.
            operators[6 + channel.spin] = adjoint @ space.to_pw(
                channel.source_grid, coordinates)
            branches[branch].append({
                "identity": channel.identity,
                "source_energies": channel.source_energies,
                "target_energies": energies,
                "source_occupations": raw_source,
                "target_occupations": _fermi((mu - energies) / width),
                "operators": operators,
                "weight": record["quadrature_weight"],
            })
    metadata = {
        "preparation": "auxiliary grand-canonical quasi-free independent reconstructed-H0 FD state",
        "common_mu_Ry": mu, "actual_degauss_Ry": width,
        "ensemble_boundary": "fixed mean N, not exact canonical N; degauss is not claimed physical temperature",
        "normalization": "A_q=Ncells^(-1/2) sum_k c_target^dag O c_source; correlation weight=1/Ncells",
        "Ncells": cells, "mesh": dimensions.tolist(),
        "reciprocal_basis_2pi_alat": np.asarray(reciprocal_basis).tolist(),
        "q_fractional": q_fractional.tolist(),
        "q_cart_2pi_alat": np.asarray(model.q_cart_2pi_alat).tolist(),
        "source_channels": records,
        "geometry_validation": "complete unshifted mesh modulo reciprocal vectors separately for each explicit spin; compatible q",
        "observables": ["g_0", "g_1", "g_2", "g_3", "g_4", "g_5", "number_up", "number_down"],
        "operator_units": ["Ry/bohr"] * 6 + ["dimensionless"] * 2,
        "gradient_sign": "g=+partial h/partial u; mixed gradient-density is not force-density until g is negated",
        "density_probe": "spin-selected periodic-uniform number Fourier probe, not pointwise density covariance; no extra volume/spin factor",
        "energy_unit": "Ry", "time_unit": "tau=t*Ry/hbar",
        "evolution": "isolated same H0, no eta or bath rate",
        "probabilities": "raw stable FD from each actual source/complete-target energy; quadrature weight stored separately",
        "zero_gap_policy": "retain every directed transition including exactly zero gaps and degenerate coherences; no division or frequency binning",
        "source_support": "inherits frame f_FD!=0 float64 source support; no additional occupation threshold; complete targets",
        "branches": "independent actual plus/minus target spectra, PW maps, eigenvectors and positive gradients",
        "missing": ["coupled electron-ion equilibrium", "irreversible bath or Gaussian higher moments",
                    "general spatial covariance beyond supplied q sectors", "matched material reference/operators and material transfer"],
        "redundant_coordinates": "no independent noise assigned to redundant delta-psi amplitudes",
        "model_branch": model.metadata.get("branch"),
        "input_manifests": model.metadata.get("input_manifests", {}),
    }
    return branches, metadata
