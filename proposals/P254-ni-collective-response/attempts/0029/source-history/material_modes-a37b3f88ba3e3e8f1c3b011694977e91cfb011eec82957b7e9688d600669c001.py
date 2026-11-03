"""Convert unweighted QE force tensors to harmonic operators without hidden repairs.

Source convention: QE 7.5 constants.f90 and PHonon/PH/dyndia.f90.
Cartesian ordering is (atom, cartesian); complex q matrices remain complex.
No mass, unit conversion, projection or instability clipping is implicit.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

# NIST2018 decimal values copied from QE 7.5 Modules/constants.f90:31-39.
H_PLANCK_SI = 6.62607015e-34
ELECTRONMASS_SI = 9.1093837015e-31
HARTREE_SI = 4.3597447222071e-18
RYDBERG_SI = HARTREE_SI / 2
BOHR_RADIUS_SI = 0.529177210903e-10
AMU_SI = 1.66053906660e-27
C_SI = 299792458
AMU_RY = AMU_SI / ELECTRONMASS_SI / 2
AU_SEC = H_PLANCK_SI / (2 * np.pi) / HARTREE_SI
RY_TO_THZ = 1 / (AU_SEC * 1e12) / (4 * np.pi)


def atom_masses_amu(geometry):
    """Resolve per-atom masses, retaining nonuniform native type assignments."""
    types = {entry["index"]: entry["mass_amu"] for entry in geometry["types"]}
    return validated_masses([types[atom["type_index"]] for atom in geometry["atoms"]])


def validated_masses(masses_amu):
    """Return a finite, strictly positive one-dimensional amu array."""
    masses = np.asarray(masses_amu, dtype=float)
    if masses.ndim != 1 or not masses.size or not np.all(np.isfinite(masses)):
        raise ValueError("masses must be a nonempty finite one-dimensional array")
    if np.any(masses <= 0):
        raise ValueError("masses must be positive")
    return masses


def force_tensor(matrix_payload):
    """Decode a backend real/imag payload, preserving [i,a,j,b] axes."""
    real = np.asarray(matrix_payload["real"], dtype=float)
    imag = np.asarray(matrix_payload["imag"], dtype=float)
    if real.shape != imag.shape:
        raise ValueError("real/imag tensor shapes differ")
    return validated_force_tensor(real + 1j * imag)


def validated_force_tensor(force):
    """Validate an unweighted Ry/Bohr^2 tensor; never average or modify it."""
    force = np.asarray(force, dtype=complex)
    if force.ndim != 4 or force.shape[1::2] != (3, 3):
        raise ValueError("force tensor axes must be [atom_i,3,atom_j,3]")
    if force.shape[0] != force.shape[2] or not force.shape[0]:
        raise ValueError("force tensor requires equal nonzero atom dimensions")
    if not np.all(np.isfinite(force)):
        raise ValueError("force tensor must be finite")
    return force


def mass_weighted_operator(force_ry_per_bohr2, masses_amu, *, units):
    """Return M^-1/2 C M^-1/2 in s^-2 ('si') or QE dyndia units ('qe_ry').

    The SI route uses C_SI=C*Ry_J/a0_m^2 and M_SI=m_amu*u_kg.
    The QE route divides by sqrt(m_i*m_j)*AMU_RY, exactly as dyndia.
    These routes use independently rounded source constants.
    """
    force = validated_force_tensor(force_ry_per_bohr2)
    masses = validated_masses(masses_amu)
    if len(masses) != force.shape[0]:
        raise ValueError("one mass is required per atom")
    scale = {"si": RYDBERG_SI / BOHR_RADIUS_SI**2 / AMU_SI,
             "qe_ry": 1 / AMU_RY}[units]
    root = np.repeat(np.sqrt(masses), 3)
    return force.reshape(3 * len(masses), -1) * scale / root[:, None] / root[None, :]


def translation_basis(masses_amu):
    """Orthonormal mass-weighted Gamma translations T_ia,b=sqrt(m_i/sum m)δ_ab."""
    masses = validated_masses(masses_amu)
    return (np.sqrt(masses / masses.sum())[:, None, None] * np.eye(3)).reshape(-1, 3)


def validated_operator(operator):
    """Validate a finite square Cartesian operator without modifying entries."""
    matrix = np.asarray(operator, dtype=complex)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or not matrix.size:
        raise ValueError("operator must be nonempty and square")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("operator must be finite")
    return matrix


def operator_residuals(operator, masses_amu):
    """Report absolute raw residuals in input operator units, including non-Gamma probes.

    Translation annihilation is physically required only at Gamma. A non-Gamma
    translation residual is a diagnostic, not a reason to project a matrix.
    """
    matrix = validated_operator(operator)
    translations = translation_basis(masses_amu)
    if matrix.shape[0] != translations.shape[0]:
        raise ValueError("operator and masses disagree")
    return {"hermitian_frobenius": float(np.linalg.norm(matrix - matrix.conj().T)),
            "translation_right_frobenius": float(np.linalg.norm(matrix @ translations)),
            "translation_left_frobenius": float(np.linalg.norm(translations.conj().T @ matrix)),
            "operator_frobenius": float(np.linalg.norm(matrix))}


def project_hermitian(operator):
    """Explicit Frobenius Hermitian projection, matching QE dyndia:56-63."""
    matrix = validated_operator(operator)
    return (matrix + matrix.conj().T) / 2


def project_gamma_translations(operator, masses_amu, *, q_fractional):
    """Explicit P D P acoustic projection at Gamma, not QE's real-space ASR.

    Caller records raw residuals first. This finite-matrix projection changes the
    model; it neither clips eigenvalues nor projects rotational sum rules.
    Exactly Hermitian inputs retain that invariant after multiplication roundoff;
    asymmetric inputs are not silently Hermitian-projected.
    """
    q = np.asarray(q_fractional, dtype=float)
    if q.shape != (3,) or not np.all(np.isfinite(q)) or np.any(q != np.rint(q)):
        raise ValueError("acoustic projection is restricted to exact Gamma modulo reciprocal lattice")
    matrix = validated_operator(operator)
    translations = translation_basis(masses_amu)
    if matrix.shape[0] != translations.shape[0]:
        raise ValueError("operator and masses disagree")
    projector = np.eye(matrix.shape[0]) - translations @ translations.T
    projected = projector @ matrix @ projector
    if np.array_equal(matrix, matrix.conj().T):
        return project_hermitian(projected)
    return projected


@dataclass(frozen=True)
class SignedModes:
    """Signed eigenvalues, signed ordinary THz frequencies, and weighted eigenvectors."""

    eigenvalues: np.ndarray
    frequency_thz: np.ndarray
    mass_weighted_vectors: np.ndarray
    eigenpair_residual: float
    orthogonality_residual: float


def signed_frequencies(eigenvalues, *, units):
    """Negative eigenvalues map to negative sqrt(abs(lambda)); none are clipped."""
    values = np.asarray(eigenvalues, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("eigenvalues must be finite")
    factor = {"si": 1 / (2 * np.pi * 1e12), "qe_ry": RY_TO_THZ}[units]
    return np.sign(values) * np.sqrt(np.abs(values)) * factor


def signed_modes(operator, *, units):
    """Diagonalize an explicitly Hermitian operator; reject unprojected asymmetry.

    Degenerate eigenvectors have no canonical orientation; compare spectra or
    eigenspace projectors, never individual vectors in a degenerate cluster.
    """
    matrix = validated_operator(operator)
    if not np.array_equal(matrix, matrix.conj().T):
        raise ValueError("operator is not exactly Hermitian; inspect raw residuals and explicitly project")
    values, vectors = np.linalg.eigh(matrix)
    return SignedModes(values, signed_frequencies(values, units=units), vectors,
                       float(np.linalg.norm(matrix @ vectors - vectors * values)),
                       float(np.linalg.norm(vectors.conj().T @ vectors - np.eye(len(values)))))
