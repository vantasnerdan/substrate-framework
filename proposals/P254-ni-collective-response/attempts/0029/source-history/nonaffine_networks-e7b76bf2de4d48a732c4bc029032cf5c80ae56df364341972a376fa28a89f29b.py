"""Central-force network elasticity by exact second-variation operators.

Conventions (all arrays float64; no iterative minimization anywhere):

- periodic cubic lattice, spacing ``a=1``, one site per cell;
- nearest-neighbour bonds plus ``f`` body-diagonal families, springs
  linear with stiffness ``k``, stress-free rest lengths;
- simple shear ``gamma`` deforms the cell map
  ``A = [a(ex + gamma*ey), a*ey, a*ez]``;
- affine modulus: ``mu_affine = k/V sum_b alpha_b**2`` with per-bond
  coupling ``alpha_b = d_hat . (dd d/dgamma) = a * D_x * n_y``;
- relaxed modulus (exact for harmonic springs, includes all non-affine
  relaxation): ``mu_relaxed = mu_affine - G.H+ G / V`` where
  ``G_i = k sum_b alpha_b n_b[i]`` is the strain-displacement coupling
  and ``H = k sum_b n_b n_b^T`` the zero-stress displacement Hessian,
  pseudoinverted on its support (null space projected, ratio reported);
- phonon route: monatomic Bloch matrix ``D_ab(q) = 2k sum n_a n_b
  (1 - cos q.R_f)``; ``c**2(q) = lambda_min(q)/q**2`` is extrapolated in
  ``q**2`` to the continuum limit.

Small-ratio practice: the relaxed correction is evaluated on the
pseudoinverse support with the null-space ratio reported; results are
extrapolated in system size; BLAS threads are pinned by the caller
before numpy import when bit-stability matters.
"""

from __future__ import annotations

from itertools import product
from typing import Any

import numpy as np

AXIS_FAMILIES = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
DIAGONAL_FAMILIES = (
    (+1, +1, +1),
    (+1, +1, -1),
    (+1, -1, +1),
    (-1, +1, +1),
)
NULLSPACE_RATIO = 1e-10


class CentralForceNetwork:
    """Periodic cubic central-force network with declared bond families."""

    def __init__(self, n: int, diagonal_families: int = 0, stiffness: float = 1.0):
        if n < 2:
            raise ValueError("lattice size must exceed 1")
        if not 0 <= diagonal_families <= len(DIAGONAL_FAMILIES):
            raise ValueError("diagonal_families must be between 0 and 4")
        self.n = int(n)
        self.stiffness = float(stiffness)
        self.spacing = 1.0
        self.offsets = AXIS_FAMILIES + DIAGONAL_FAMILIES[:diagonal_families]
        self.rest_per_family = np.array(
            [float(np.linalg.norm(off)) * self.spacing for off in self.offsets],
            dtype=np.float64,
        )
        sites = list(product(range(n), repeat=3))
        index = {site: i for i, site in enumerate(sites)}
        sources, targets, deltas_lattice, families = [], [], [], []
        for site in sites:
            for fi, offset in enumerate(self.offsets):
                raw = [site[d] + offset[d] for d in range(3)]
                sources.append(index[site])
                targets.append(index[tuple(r % self.n for r in raw)])
                # The physical bond vector is the family offset itself;
                # wrapping only selects the target's periodic image.
                deltas_lattice.append(np.array(offset, dtype=np.float64))
                families.append(fi)
        self.sources = np.array(sources, dtype=int)
        self.targets = np.array(targets, dtype=int)
        self.deltas_lattice = np.array(deltas_lattice, dtype=np.float64)
        self.family_index = np.array(families, dtype=int)
        self.rest_lengths = self.rest_per_family[self.family_index]
        self.unit_vectors = self.deltas_lattice / np.linalg.norm(
            self.deltas_lattice, axis=1
        )[:, None]


    @property
    def num_sites(self) -> int:
        return self.positions_frac_shape()[0]

    def positions_frac_shape(self) -> tuple[int, int]:
        return (self.n**3, 3)

    @property
    def volume(self) -> float:
        return float(self.n**3 * self.spacing**3)

    @property
    def mean_coordination(self) -> float:
        return float(2.0 * len(self.sources) / self.n**3)

    def affine_shear_modulus(self) -> float:
        """``mu_affine = k/V sum_b alpha_b**2``, ``alpha = a * D_x * n_y``."""

        alpha = self.spacing * self.deltas_lattice[:, 0] * self.unit_vectors[:, 1]
        return float(self.stiffness * np.sum(alpha**2) / self.volume)
    def strain_coupling(self) -> np.ndarray:
        """``G_i = k sum_b alpha_b n_b[i]`` as a flat 3N vector."""

        alpha = self.spacing * self.deltas_lattice[:, 0] * self.unit_vectors[:, 1]
        coupling = np.zeros((self.n**3, 3), dtype=np.float64)
        weighted = self.stiffness * alpha[:, None] * self.unit_vectors
        np.add.at(coupling, self.sources, weighted)
        np.add.at(coupling, self.targets, weighted)
        return coupling.ravel()

    def hessian(self) -> np.ndarray:
        """Zero-stress displacement Hessian ``k sum_b n n^T`` blocks."""

        size = 3 * self.n**3
        matrix = np.zeros((size, size), dtype=np.float64)
        blocks = self.stiffness * np.einsum(
            "bi,bj->bij", self.unit_vectors, self.unit_vectors
        )
        for source, target, block in zip(self.sources, self.targets, blocks):
            i0, j0 = 3 * int(source), 3 * int(target)
            matrix[i0 : i0 + 3, i0 : i0 + 3] += block
            matrix[j0 : j0 + 3, j0 : j0 + 3] += block
            matrix[i0 : i0 + 3, j0 : j0 + 3] -= block
            matrix[j0 : j0 + 3, i0 : i0 + 3] -= block
        return matrix

    def relaxed_shear_modulus(self) -> dict[str, Any]:
        """Exact relaxed modulus via the second-variation operator.

        Returns the modulus, the non-affine correction, and the Hessian
        conditioning ratio lambda_min_support/lambda_2 on the support.
        """

        matrix = self.hessian()
        eigenvalues, eigenvectors = np.linalg.eigh(matrix)
        scale = float(eigenvalues.max())
        if scale <= 0.0:
            return {
                "modulus": self.affine_shear_modulus(),
                "nonaffine_correction": 0.0,
                "nullspace_ratio": 0.0,
            }
        support = eigenvalues > NULLSPACE_RATIO * scale
        coupling = self.strain_coupling()
        reduced = eigenvectors[:, support].T @ coupling
        correction = float(np.sum(reduced**2 / eigenvalues[support]))
        affine = self.affine_shear_modulus()
        modulus = affine - correction / self.volume
        supported = eigenvalues[support]
        ratio = float(supported.min() / np.sort(supported)[-2])
        return {
            "modulus": modulus,
            "nonaffine_correction": correction / self.volume,
            "nullspace_ratio": ratio,
        }


def bloch_matrix(net: CentralForceNetwork, q_vector: np.ndarray) -> np.ndarray:
    """Exact monatomic Bloch matrix ``2k sum n n^T (1-cos q.R_f)``."""

    matrix = np.zeros((3, 3), dtype=np.float64)
    for offset in net.offsets:
        ref = np.array(offset, dtype=np.float64) * net.spacing
        unit = ref / np.linalg.norm(ref)
        weight = 1.0 - np.cos(float(np.dot(q_vector, ref)))
        matrix += 2.0 * net.stiffness * weight * np.outer(unit, unit)
    return matrix


def acoustic_speed_squared(
    net: CentralForceNetwork, mass: float = 1.0
) -> dict[str, Any]:
    """Transverse acoustic speed squared along x from the Bloch matrix.

    Selects the smallest positive eigenvalue whose eigenvector is
    orthogonal to the propagation direction. Exactly soft transverse
    modes report ``0.0``; when anisotropic hybridization leaves no pure
    transverse eigenvector at finite q, the value is ``nan`` and
    ``transverse_branch_found`` is false.
    """

    q_mag = 2.0 * np.pi / (net.n * net.spacing)
    q_vector = np.array([q_mag, 0.0, 0.0])
    matrix = bloch_matrix(net, q_vector) / mass
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    scale = max(float(eigenvalues.max()), 1.0)
    transverse = [
        (float(value), vector)
        for value, vector in zip(eigenvalues, eigenvectors.T)
        if abs(float(np.dot(vector, np.array([1.0, 0.0, 0.0])))) < 1e-8
    ]
    if not transverse:
        return {
            "c_squared_finite_q": float("nan"),
            "transverse_branch_found": False,
            "q_squared": q_mag**2,
        }
    softest = min(value for value, _ in transverse)
    if softest <= 1e-12 * scale:
        omega_min_sq = 0.0
    else:
        omega_min_sq = softest
    return {
        "c_squared_finite_q": omega_min_sq / q_mag**2,
        "transverse_branch_found": True,
        "q_squared": q_mag**2,
    }




def extrapolated_acoustic_speed_squared(
    diagonal_families: int, sizes: tuple[int, ...] = (3, 6), mass: float = 1.0
) -> dict[str, Any]:
    """Continuum acoustic speed squared via the exact dispersion variable.

    Every branch eigenvalue here is a sum of ``weight_f*(1 - cos q.R_f)``
    terms, so ``c**2`` is an exact affine function of
    ``t(q) = (1 - cos q)/q**2`` along x. Two distinct sizes determine the
    line; the continuum limit is its value at ``t -> 1/2``.
    """

    points = []
    for n in sizes:
        net = CentralForceNetwork(n, diagonal_families=diagonal_families)
        data = acoustic_speed_squared(net, mass=mass)
        q_squared = data["q_squared"]
        t_value = (1.0 - np.cos(np.sqrt(q_squared))) / q_squared
        points.append((float(t_value), data["c_squared_finite_q"]))
    t_values = np.array([p[0] for p in points], dtype=np.float64)
    c_values = np.array([p[1] for p in points], dtype=np.float64)
    slope, intercept = np.polyfit(t_values, c_values, 1)
    continuum = float(intercept + slope * 0.5)
    return {
        "c_squared_continuum": continuum,
        "slope": float(slope),
        "points": points,
    }

def threshold_record(n: int = 4, max_families: int = 4) -> list[dict[str, Any]]:
    """Affine and relaxed moduli plus phonon data across family counts."""

    records: list[dict[str, Any]] = []
    for families in range(max_families + 1):
        net = CentralForceNetwork(n, diagonal_families=families)
        relaxed = net.relaxed_shear_modulus()
        records.append(
            {
                "families": families,
                "z": net.mean_coordination,
                "mu_affine": net.affine_shear_modulus(),
                **relaxed,
                **acoustic_speed_squared(net),
            }
        )
    return records
