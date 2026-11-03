"""Complete native q-star grids and their discrete periodic force lattice.

QE 7.5 do_q2r:234-299 and matdyn frc_blk:1196-1220 fix the convention:
C_ia,jb(q)=sum_R C_ia,jb(R) exp(-2π i q·R), R=lattice(row)-lattice(column).
No tau phase occurs in this transform. The oriented geometric separation is
R+tau_i-tau_j (matdyn:1171-1173, ordinary PH convention fd=false).
The finite transform is an aliased periodic force lattice, not an infinite pair
potential, a Wigner-Seitz off-grid interpolation, or a polar correction.
"""

from __future__ import annotations

import numpy as np

from .material_modes import force_tensor


def validated_grid_shape(grid_shape):
    """Require an explicit nonempty three-dimensional integer grid."""
    shape = np.asarray(grid_shape)
    if shape.shape != (3,) or np.any(shape <= 0) or np.any(shape != shape.astype(int)):
        raise ValueError("grid_shape must contain three positive integers")
    return tuple(int(n) for n in shape)


def complete_native_q_grid(matrix_records, grid_shape, *, grid_atol):
    """Expand already-native q-star members into every unshifted grid slot exactly once.

    No symmetry rotation is invented: every member must carry its own raw matrix.
    q_cart is in 2π/alat and direct lattice ROWS are in alat, so q_frac=at@q_cart.
    grid_atol is an explicit absolute tolerance on q_frac*grid_shape; duplicates
    and incomplete stars fail rather than being averaged, filled or inferred.
    """
    shape = validated_grid_shape(grid_shape)
    if not 0 <= grid_atol < 0.5:
        raise ValueError("grid_atol must be finite, nonnegative and below half a grid step")
    records = list(matrix_records)
    if not records:
        raise ValueError("no matrix records")
    geometry = records[0]["geometry"]
    at = np.asarray(geometry["at_cartesian_by_lattice_vector_alat"], dtype=float)
    if at.shape != (3, 3) or not np.all(np.isfinite(at)) or np.linalg.det(at) == 0:
        raise ValueError("direct lattice must be finite and nonsingular")
    nat = geometry["number_of_atoms"]
    grid = np.empty(shape + (nat, 3, nat, 3), dtype=complex)
    occupied = np.zeros(shape, dtype=bool)
    for record in records:
        if record["geometry"] != geometry:
            raise ValueError("q records do not share identical native geometry and masses")
        for member in record["q_star"]:
            q = np.asarray(member["q_cartesian_2pi_alat"], dtype=float)
            if q.shape != (3,) or not np.all(np.isfinite(q)):
                raise ValueError("invalid q vector")
            coordinate = (at @ q) * shape
            nearest = np.rint(coordinate)
            if np.max(np.abs(coordinate - nearest)) > grid_atol:
                raise ValueError("native q vector is not on the declared unshifted grid")
            index = tuple((nearest.astype(int) % shape).tolist())
            if occupied[index]:
                raise ValueError("duplicate native q-star grid member")
            tensor = force_tensor(member["force_matrix_ry_per_bohr2"])
            if tensor.shape != (nat, 3, nat, 3):
                raise ValueError("matrix atom count disagrees with geometry")
            grid[index] = tensor
            occupied[index] = True
    if not np.all(occupied):
        raise ValueError(f"incomplete native q grid: {int(occupied.sum())}/{occupied.size}")
    return grid


def validated_force_grid(values):
    """Validate [n1,n2,n3,atom_i,cart_i,atom_j,cart_j] finite complex data."""
    grid = np.asarray(values, dtype=complex)
    if grid.ndim != 7 or grid.shape[4] != 3 or grid.shape[6] != 3:
        raise ValueError("expected three grid dimensions followed by [atom_i,3,atom_j,3]")
    if min(grid.shape[:4]) < 1 or grid.shape[3] != grid.shape[5]:
        raise ValueError("nonzero grid and equal atom dimensions required")
    if not np.all(np.isfinite(grid)):
        raise ValueError("force grid must be finite")
    return grid


def q_grid_to_force_lattice(force_q):
    """C(R)=N^-1 sum_q C(q) exp(+2π i q·R); no Hermitian/acoustic repair."""
    return np.fft.ifftn(validated_force_grid(force_q), axes=(0, 1, 2))


def force_lattice_to_q_grid(force_r):
    """Inverse reconstruction on the original complete grid only."""
    return np.fft.fftn(validated_force_grid(force_r), axes=(0, 1, 2))


def lattice_residuals(force_r):
    """Report raw reality, reciprocal-pair and Gamma translation residuals in Ry/Bohr².

    Reciprocity: C_ia,jb(R)=conj(C_jb,ia(-R)). Reality is reported separately;
    no imaginary part or residual is discarded by this interface.
    """
    force = validated_force_grid(force_r)
    opposite = force
    for axis, n in enumerate(force.shape[:3]):
        opposite = np.take(opposite, (-np.arange(n)) % n, axis=axis)
    reciprocal = opposite.transpose(0, 1, 2, 5, 6, 3, 4).conj()
    gamma = force.sum(axis=(0, 1, 2))
    return {"imaginary_frobenius": float(np.linalg.norm(force.imag)),
            "reciprocity_frobenius": float(np.linalg.norm(force - reciprocal)),
            "translation_frobenius": float(np.linalg.norm(gamma.sum(axis=2)))}
