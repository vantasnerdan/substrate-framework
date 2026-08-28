"""P247 attempt 0004, gate C3: extended 4x4 radial functional (torch).

Evaluates the conditional P239 Hamiltonian density on the boosted radial
ansatz M(x) = R(chi(x)) blkdiag(0, S(x)) R(chi(x))^T, with S the committed
biaxial-hedgehog spectral form and chi a fourth radial modal channel (boost
rapidity along the local director).

Structure (each term reduces exactly to the committed functional at chi=0):

  curvature   4 * sum_{l<r} ||[D_l, D_r]_eta||^2      (D_mu = d_mu M)
  inertia     4 * sum_i     ||[W, D_i]_eta||^2        (W = boosted clock
              response blkdiag(0, [N_hat, S]_resp); Omega eliminated, the
              committed J = 1 fixed-J convention E_J = 1/(4 I))
  potential   committed pinned-spectrum potential on the 3x3 block
  sigma       kappa * sum_i q(d_i u, d_i u),  u = R(-chi)^T xi the boosted
              0-eigenline, q(a, b) = -a . (eta b);  vanishes at chi = 0
  boost mass  (Lambda^2/2) sinh^2(chi)                  (candidate B, B1)

Constants: kappa = 1 (C-M5S-001 census normalization; boost metric entry
1/16 verified in c1b-kappa.json), Lambda^2 = 0.26847204181661866
(C-M5S-002).  E_J = E_stat + 1/(4 I), omega = 1/(2 I) (committed).
"""

from __future__ import annotations

import sys
from pathlib import Path

import torch

HERE = Path(__file__).resolve().parent
A1 = HERE.parent / "0001"
sys.path.insert(0, str(A1))
sys.path.insert(0, str(HERE.parent / "0002"))

import debox_common as base  # noqa: E402  (installs P240 paths, pins threads)

from cpu_energy import (  # noqa: E402
    DEVICE,
    DTYPE,
    chebyshev_stack,
    frobenius_squared,
    gauss_grid,
)

ETA_INV = torch.tensor(
    [[-1.0, 0, 0, 0], [0, 1.0, 0, 0], [0, 0, 1.0, 0], [0, 0, 0, 1.0]],
    dtype=DTYPE,
    device=DEVICE,
)
ETA = ETA_INV.clone()
LAMBDA_SQ = 0.26847204181661866  # C-M5S-002
KAPPA = 1.0  # C-M5S-001 census normalization


def boosted_boost_matrix(director: torch.Tensor, chi: torch.Tensor):
    """Return (R, R^-T) for boost of rapidity chi along local director.

    R boosts the 0-axis into the director; the 0-eigenline of the boosted
    block field is u = R^{-T} xi = (cosh chi, -sinh chi * director).
    """
    c, s = torch.cosh(chi), torch.sinh(chi)
    nx, ny, nz = director.unbind(-1)
    zero = torch.zeros_like(c)
    one = torch.ones_like(c)
    top = torch.stack((c, s * nx, s * ny, s * nz), dim=-1)
    row1 = torch.stack(
        (s * nx, one + (c - 1) * nx * nx, (c - 1) * nx * ny, (c - 1) * nx * nz), dim=-1
    )
    row2 = torch.stack(
        (s * ny, (c - 1) * ny * nx, one + (c - 1) * ny * ny, (c - 1) * ny * nz), dim=-1
    )
    row3 = torch.stack(
        (s * nz, (c - 1) * nz * nx, (c - 1) * nz * ny, one + (c - 1) * nz * nz), dim=-1
    )
    R = torch.stack((top, row1, row2, row3), dim=-2)
    # R^{-T} for a boost along n: boost of rapidity -chi along n
    RmT = boost_matrix_only(director, -chi)
    return R, RmT


def boost_matrix_only(director: torch.Tensor, chi: torch.Tensor) -> torch.Tensor:
    c, s = torch.cosh(chi), torch.sinh(chi)
    nx, ny, nz = director.unbind(-1)
    zero = torch.zeros_like(c)
    one = torch.ones_like(c)
    top = torch.stack((c, s * nx, s * ny, s * nz), dim=-1)
    row1 = torch.stack(
        (s * nx, one + (c - 1) * nx * nx, (c - 1) * nx * ny, (c - 1) * nx * nz), dim=-1
    )
    row2 = torch.stack(
        (s * ny, (c - 1) * ny * nx, one + (c - 1) * ny * ny, (c - 1) * ny * nz), dim=-1
    )
    row3 = torch.stack(
        (s * nz, (c - 1) * nz * nx, (c - 1) * nz * ny, one + (c - 1) * nz * nz), dim=-1
    )
    return torch.stack((top, row1, row2, row3), dim=-2)


def embed_spatial(spatial: torch.Tensor) -> torch.Tensor:
    shape = spatial.shape[:-2] + (4, 4)
    big = torch.zeros(shape, dtype=spatial.dtype, device=spatial.device)
    big[..., 1:, 1:] = spatial
    return big


def eta_commutator(a: torch.Tensor, b: torch.Tensor) -> torch.Tensor:
    return a @ ETA_INV @ b - b @ ETA_INV @ a


def hedgehog_spatial(c_spatial: torch.Tensor, radial_order: int, radius_grid, mu_grid, radius: float):
    """Committed hedgehog spectral form; c_spatial shape (3, order, 1)."""
    normalized = radius_grid / radius
    radial_coordinate = 2 * normalized**2 - 1
    radial_basis = chebyshev_stack(radial_coordinate, tuple(range(radial_order)))
    angular_basis = chebyshev_stack(mu_grid, (0,))
    modal = torch.einsum("...i,cij,...j->...c", radial_basis, c_spatial, angular_basis)
    q = normalized**2 + normalized**2 * (1 - normalized**2) * modal[..., 0]
    tangent = (1 - normalized**2) * (
        torch.tensor(1 / 3, dtype=DTYPE, device=DEVICE) + modal[..., 1]
    )
    split_amplitude = normalized**4 * (1 - normalized**2) * modal[..., 2]
    sine = torch.sqrt(torch.clamp(1 - mu_grid**2, min=0.0))
    delta = split_amplitude * sine**2
    zero = torch.zeros_like(sine)
    director = torch.stack((sine, zero, mu_grid), dim=-1)
    polar = torch.stack((mu_grid, zero, -sine), dim=-1)
    azimuthal = torch.stack((zero, torch.ones_like(zero), zero), dim=-1)

    def outer(vector):
        return vector[..., :, None] * vector[..., None, :]

    lambda_n = tangent + q
    spatial = (
        lambda_n[..., None, None] * outer(director)
        + (tangent + delta)[..., None, None] * outer(polar)
        + (tangent - delta)[..., None, None] * outer(azimuthal)
    )
    return spatial, director


def extended_static(
    flat: torch.Tensor,
    *,
    radial_order: int,
    radial_nodes: int,
    angular_nodes: int,
    radius: float,
):
    """Extended static energy (no explicit Omega; inertia returned for the
    fixed-J term).  flat = (4, radial_order) coefficient rows: committed
    q, tangent, split channels then the chi channel."""
    channels = flat.reshape(4, radial_order, 1)
    c_spatial = channels[:3]
    c_chi = channels[3:4]

    radial, radial_weight, mu, angular_weight = gauss_grid(
        radial_nodes, angular_nodes, radius
    )
    radius_grid = radial[:, None].repeat(1, angular_nodes).clone().requires_grad_(True)
    mu_grid = mu[None, :].repeat(radial_nodes, 1).clone().requires_grad_(True)
    normalized = radius_grid / radius

    spatial, director = hedgehog_spatial(c_spatial, radial_order, radius_grid, mu_grid, radius)

    radial_coordinate = 2 * normalized**2 - 1
    radial_basis = chebyshev_stack(radial_coordinate, tuple(range(radial_order)))
    chi = (radial_basis @ c_chi)[..., 0]
    chi = chi * normalized * (1 - normalized)

    R, RmT = boosted_boost_matrix(director, chi)
    M = R @ embed_spatial(spatial) @ R.transpose(-1, -2)
    u = RmT @ torch.tensor([1.0, 0, 0, 0], dtype=DTYPE, device=DEVICE)

    spatial_two = spatial @ spatial
    trace_two = torch.diagonal(spatial_two, dim1=-2, dim2=-1).sum(-1)
    trace_three = torch.diagonal(spatial_two @ spatial, dim1=-2, dim2=-1).sum(-1)
    potential = -0.5 * trace_two - trace_three + trace_two**2 + 0.5
    s_m = 0.5 * LAMBDA_SQ * torch.sinh(chi) ** 2

    def spatial_derivative(coordinate):
        rows = []
        for left in range(4):
            columns = []
            for right in range(4):
                columns.append(
                    torch.autograd.grad(
                        M[..., left, right],
                        coordinate,
                        torch.ones_like(M[..., left, right]),
                        retain_graph=True,
                        create_graph=True,
                    )[0]
                )
            rows.append(torch.stack(columns, dim=-1))
        return torch.stack(rows, dim=-2)

    derivative_r = spatial_derivative(radius_grid)
    derivative_mu = spatial_derivative(mu_grid)
    sine = torch.sqrt(torch.clamp(1 - mu_grid**2, min=0.0))
    derivative_theta = -sine[..., None, None] * derivative_mu / radius_grid[..., None, None]
    rotation_z = torch.tensor(
        [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype=DTYPE, device=DEVICE
    )
    big_rot = torch.zeros(4, 4, dtype=DTYPE, device=DEVICE)
    big_rot[1:, 1:] = rotation_z
    derivative_phi = (big_rot @ M - M @ big_rot) / (radius_grid * sine)[..., None, None]
    derivatives = (derivative_r, derivative_theta, derivative_phi)

    curvature_density = 4 * sum(
        frobenius_squared(eta_commutator(derivatives[l], derivatives[r]))
        for l in range(3)
        for r in range(l + 1, 3)
    )

    # boosted clock response W = R blkdiag(0, [N_hat, S]_resp) R^T
    nx, ny, nz = director.unbind(-1)
    zero = torch.zeros_like(nx)
    generator = torch.stack(
        (
            torch.stack((zero, -nz, ny), dim=-1),
            torch.stack((nz, zero, -nx), dim=-1),
            torch.stack((-ny, nx, zero), dim=-1),
        ),
        dim=-2,
    )
    response = generator @ spatial + spatial @ generator.transpose(-1, -2)
    W = R @ embed_spatial(response) @ R.transpose(-1, -2)
    inertia_density = 4 * sum(
        frobenius_squared(eta_commutator(W, d)) for d in derivatives
    )

    # sigma: kappa sum_i q(d_i u, d_i u), q(a, b) = -a . (eta b)
    def u_derivative(coordinate):
        rows = []
        for comp in range(4):
            g = torch.autograd.grad(
                u[..., comp],
                coordinate,
                torch.ones_like(u[..., comp]),
                retain_graph=True,
                create_graph=True,
                allow_unused=True,
            )[0]
            if g is None:  # u independent of this coordinate (e.g. fixed axis)
                g = torch.zeros_like(u[..., comp])
            rows.append(g)
        return torch.stack(rows, dim=-1)

    du_r = u_derivative(radius_grid)
    du_theta = -sine[..., None] * u_derivative(mu_grid) / radius_grid[..., None]
    sigma_density = KAPPA * (
        -torch.sum(du_r * (du_r @ ETA), dim=-1)
        - torch.sum(du_theta * (du_theta @ ETA), dim=-1)
    )

    weights = 2 * torch.pi * radius_grid**2 * radial_weight[:, None] * angular_weight[None, :]
    curvature = torch.sum(weights * curvature_density)
    potential_energy = torch.sum(weights * potential)
    sigma_energy = torch.sum(weights * sigma_density)
    boost_mass = torch.sum(weights * s_m)
    inertia = torch.sum(weights * inertia_density)
    static = potential_energy + curvature + sigma_energy + boost_mass
    fixed_j = 1.0 / (4.0 * inertia)
    omega = 1.0 / (2.0 * inertia)
    total = static + fixed_j
    return total, {
        "curvature": curvature.detach(),
        "potential": potential_energy.detach(),
        "sigma": sigma_energy.detach(),
        "boost_mass": boost_mass.detach(),
        "inertia": inertia.detach(),
        "fixed_j": fixed_j.detach(),
        "omega": omega.detach(),
    }
