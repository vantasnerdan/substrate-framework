"""A31 frozen-initial-H NONSTATIONARY finite-bath comparator.

The sole approximation here is holding the supplied literal projected Hamiltonian
fixed. No state is supplied, diagonalized, refilled, or assumed to commute with H.
This is not stationary/FDT response, damping, physical loss, or moving/RPA closure.
Time and H must be in mutually reciprocal units (hbar=1). No material files are read.

For each supplied interval [l,r], x=(t-(l+r)/2)/((r-l)/2), and coefficients are
ordinary NumPy Chebyshev coefficients (the constant term is NOT halved). The
polynomial is integrated on every interval, including all degree-seven terms.
Raw H is never Hermitized: a general eigendecomposition and its inverse retain
its Hermiticity/eigenpair defects. The physical projected H is near-Hermitian;
the inverse-evolution formula is unitary Heisenberg evolution only for Hermitian
H, not a general non-Hermitian physical history. The supplied-adjoint leg evolves
with raw H† and independent supplied Fadj, never a replacement F†.

Source inspected before implementation: installed SciPy 1.11.4
special/_spherical_bessel.py (negative real arguments have no reflection wrapper),
NumPy 1.26.4 polynomial/legendre.py:legmulx and chebyshev.py:chebval. Explicit
reflection below also works with newer SciPy's reflection-decorated wrapper.
"""

from __future__ import annotations

import numpy as np
from numpy.polynomial.legendre import legmulx
from scipy.linalg import eig, inv
from scipy.special import spherical_jn


SCOPE = "frozen-initial-H nonstationary selected-mode finite-basis comparator"


def _chebyshev_to_legendre(count: int) -> np.ndarray:
    """Columns represent T_k in the P_n basis, without a monomial detour.

    Use T_(k+1)=2*x*T_k-T_(k-1) and the stable Legendre multiplication
    recurrence x*P_n=((n+1)*P_(n+1)+n*P_(n-1))/(2*n+1).
    """
    conversion = np.zeros((count, count), dtype=float)
    conversion[0, 0] = 1.0
    if count > 1:
        conversion[1, 1] = 1.0
    for k in range(1, count - 1):
        column = 2.0 * legmulx(conversion[: k + 1, k])
        conversion[: len(column), k + 1] = column
        conversion[:, k + 1] -= conversion[:, k - 1]
    return conversion


def _fourier_legendre_moments(gaps: np.ndarray, half_width: float,
                               midpoint: float, count: int) -> np.ndarray:
    """Integrals of P_n(x(t))*exp(i*gap*t) dt, including phase/Jacobian ONCE.

    There is no inverse-gap recurrence or tolerance-based gap clipping. The
    removable singularity at EXACT zero is supplied analytically. Reflection
    j_n(-z)=(-1)^n*j_n(z) is valid for complex as well as real z.
    """
    z = gaps * half_width
    reflect = z.real < 0
    reflected_z = np.where(reflect, -z, z)
    zero = z == 0
    phase_jacobian = half_width * np.exp(1j * gaps * midpoint)
    moments = np.empty((count,) + gaps.shape, dtype=np.complex128)
    for n in range(count):
        bessel = spherical_jn(n, reflected_z)
        if n % 2:
            bessel = np.where(reflect, -bessel, bessel)
        bessel = np.where(zero, 1.0 if n == 0 else 0.0, bessel)
        moments[n] = (2.0 * (1j ** n)) * bessel * phase_jacobian
    return moments


def _relative(defect: np.ndarray, reference: np.ndarray) -> float:
    denominator = float(np.linalg.norm(reference))
    numerator = float(np.linalg.norm(defect))
    return numerator / denominator if denominator else numerator


def integrate_filtered_operators(hamiltonian, forces, force_adjoints,
                                 segment_bounds, modal_chebyshev):
    """Integrate sum_m a_mp(t) exp(+iHt) F_m exp(-iHt) for each probe p.

    Inputs: H[b,b], F/Fadj[m,b,b], coefficients[s,d,m,p]. Bounds may be
    edges[s+1] or explicit intervals[s,2]; no interval or polynomial term is
    subsampled. Degree d is ascending and at most seven. The independent
    adjoint leg is integral conj(a_mp(t))*exp(+iH†t)*Fadj_m*exp(-iH†t).

    Return a dictionary with operators[p,b,b], adjoint_operators[p,b,b],
    diagnostics, metadata, and raw spectral data. All state contractions and
    retained/omitted cross terms are the caller's responsibility. In particular,
    do NOT replace an original noncommuting state by spectral populations.

    This spectral evaluation requires a numerically invertible eigenvector
    frame. Singular frames fail explicitly; ill conditioning and eigenpair
    residuals are reported, never repaired by clipping or Hermitization.
    """
    h = np.asarray(hamiltonian, dtype=np.complex128)
    f = np.asarray(forces, dtype=np.complex128)
    fa = np.asarray(force_adjoints, dtype=np.complex128)
    coefficients = np.asarray(modal_chebyshev, dtype=np.complex128)
    bounds = np.asarray(segment_bounds)
    if np.iscomplexobj(bounds):
        raise ValueError("segment bounds must be real")
    bounds = np.asarray(bounds, dtype=float)
    if h.ndim != 2 or h.shape[0] != h.shape[1] or not h.shape[0]:
        raise ValueError("hamiltonian must be a nonempty square matrix")
    size = h.shape[0]
    if f.ndim != 3 or f.shape[1:] != (size, size) or fa.shape != f.shape:
        raise ValueError("forces and supplied adjoints must have shape [mode,basis,basis]")
    if coefficients.ndim != 4 or coefficients.shape[2] != f.shape[0]:
        raise ValueError("coefficients must have shape [segment,degree,mode,probe]")
    segments, count, modes, probes = coefficients.shape
    if not segments or not modes or not probes or not 1 <= count <= 8:
        raise ValueError("require nonempty axes and polynomial degree at most seven")
    if bounds.shape == (segments + 1,):
        intervals = np.column_stack((bounds[:-1], bounds[1:]))
    elif bounds.shape == (segments, 2):
        intervals = bounds
    else:
        raise ValueError("bounds must be [segment+1] edges or [segment,2] intervals")
    if not all(np.all(np.isfinite(a)) for a in (h, f, fa, coefficients, intervals)):
        raise ValueError("all operator, polynomial and interval inputs must be finite")
    if np.any(intervals[:, 1] <= intervals[:, 0]):
        raise ValueError("every supplied interval must have positive width")

    energies, vectors = eig(h, check_finite=False)
    inverse_vectors = inv(vectors, check_finite=False)
    gaps = energies[:, None] - energies[None, :]
    eigen_forces = inverse_vectors @ f @ vectors
    # H† = V^(-†) E* V†; neither Fadj nor H is silently replaced/repaired.
    eigen_adjoints = vectors.conj().T @ fa @ inverse_vectors.conj().T
    conversion = _chebyshev_to_legendre(count)
    eigen_operators = np.zeros((probes, size, size), dtype=np.complex128)
    eigen_adjoint_operators = np.zeros_like(eigen_operators)
    real_gaps = np.all(gaps.imag == 0)
    for segment, (left, right) in enumerate(intervals):
        half_width = (right - left) / 2.0
        midpoint = left + half_width
        legendre = np.einsum("dk,kmp->dmp", conversion, coefficients[segment])
        moments = _fourier_legendre_moments(gaps, half_width, midpoint, count)
        adjoint_moments = moments if real_gaps else _fourier_legendre_moments(
            gaps.conj(), half_width, midpoint, count)
        for mode in range(modes):
            weights = np.einsum("dij,dp->pij", moments, legendre[:, mode, :])
            eigen_operators += weights * eigen_forces[mode]
            adjoint_weights = np.einsum(
                "dij,dp->pij", adjoint_moments, legendre[:, mode, :].conj())
            eigen_adjoint_operators += adjoint_weights * eigen_adjoints[mode]
    operators = vectors @ eigen_operators @ inverse_vectors
    adjoint_operators = inverse_vectors.conj().T @ eigen_adjoint_operators @ vectors.conj().T
    right_residual = h @ vectors - vectors * energies[None, :]
    left_residual = inverse_vectors @ h - energies[:, None] * inverse_vectors
    reconstruction_defect = (vectors * energies[None, :]) @ inverse_vectors - h
    force_defects = fa - f.conj().transpose(0, 2, 1)
    integrated_defects = adjoint_operators - operators.conj().transpose(0, 2, 1)
    identity = np.eye(size)
    diagnostics = {
        "hamiltonian_hermiticity_absolute_fro": float(np.linalg.norm(h - h.conj().T)),
        "hamiltonian_hermiticity_relative_fro": _relative(h - h.conj().T, h),
        "hamiltonian_reconstruction_absolute_fro": float(np.linalg.norm(reconstruction_defect)),
        "hamiltonian_reconstruction_relative_fro": _relative(reconstruction_defect, h),
        "right_eigenpair_absolute_fro": float(np.linalg.norm(right_residual)),
        "right_eigenpair_relative_fro": _relative(right_residual, h @ vectors),
        "left_eigenpair_absolute_fro": float(np.linalg.norm(left_residual)),
        "left_eigenpair_relative_fro": _relative(left_residual, inverse_vectors @ h),
        "eigenvector_condition_2": float(np.linalg.cond(vectors)),
        "eigenvector_inverse_defect_fro": float(np.linalg.norm(inverse_vectors @ vectors - identity)),
        "eigenvector_unitarity_defect_fro": float(np.linalg.norm(vectors.conj().T @ vectors - identity)),
        "maximum_absolute_eigenvalue_imaginary_part": float(np.max(np.abs(energies.imag))),
        "force_supplied_adjoint_absolute_fro_by_mode": np.linalg.norm(force_defects, axis=(1, 2)).tolist(),
        "force_supplied_adjoint_relative_fro_by_mode": [_relative(force_defects[m], f[m]) for m in range(modes)],
        "integrated_supplied_adjoint_absolute_fro_by_probe": np.linalg.norm(integrated_defects, axis=(1, 2)).tolist(),
        "integrated_supplied_adjoint_relative_fro_by_probe": [_relative(integrated_defects[p], operators[p]) for p in range(probes)],
        "exact_zero_energy_gap_count": int(np.count_nonzero(gaps == 0)),
        "negative_real_part_energy_gap_count": int(np.count_nonzero(gaps.real < 0)),
        "minimum_nonzero_absolute_energy_gap": float(np.min(np.abs(gaps[gaps != 0]))) if np.any(gaps != 0) else None,
        "all_outputs_finite": bool(np.all(np.isfinite(operators)) and np.all(np.isfinite(adjoint_operators))),
    }
    return {
        "operators": operators,
        "adjoint_operators": adjoint_operators,
        "diagnostics": diagnostics,
        "metadata": {
            "approximation_scope": SCOPE,
            "hamiltonian_policy": "literal raw H, general eigenvectors and inverse; no Hermitization/clipping",
            "adjoint_evolution": "raw H.conj().T with independent supplied force_adjoints and conjugated coefficients",
            "state_policy": "no state consumed or altered; original noncommuting state remains unrestricted",
            "excluded_claims": ["stationary response", "FDT", "thermal refill", "damping", "full moving/RPA closure", "physical loss", "full campaign"],
            "integration": "Chebyshev-to-Legendre basis recurrence; 2*i**n*spherical_jn(n,gap*half_width)",
            "time_phase_convention": "+i*(energy_row-energy_column)*t; hbar=1",
            "segment_count": segments,
            "polynomial_degree": count - 1,
            "mode_count": modes,
            "probe_count": probes,
            "basis_count": size,
            "intervals": intervals.tolist(),
            "numerical_scope": "floating-point spectral evaluation, not rigorous enclosure or a finite-basis error bound",
        },
        "spectral": {
            "eigenvalues": energies,
            "right_eigenvectors": vectors,
            "inverse_right_eigenvectors": inverse_vectors,
            "energy_gaps": gaps,
            "eigenbasis_operators": eigen_operators,
            "eigenbasis_adjoint_operators": eigen_adjoint_operators,
            "chebyshev_to_legendre": conversion,
        },
    }
