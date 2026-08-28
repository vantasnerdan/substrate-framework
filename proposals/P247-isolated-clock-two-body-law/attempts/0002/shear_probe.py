"""P247 attempt 0002: shear-channel curvature probe (potential AND static).

The G1 construction proves the np and na shear channels are exactly flat in
the POTENTIAL part at the rank-1 projector background n-hat n-hat^T. This
probe evaluates, along those shear directions:

- the potential-part quadratic coefficient a2_pot, and
- the full static quadratic coefficient a2_static (curvature/gradient term
  included; candidate A's premise is that the shear channels carry
  gradient-only stiffness),

at (i) the projector background - analytically zero in the potential part,
so its numeric value is the probe's own floor - and (ii) each de-boxed root
R=12..18, recording where the exactly-flat channels land on the actual branch
backgrounds.

Note: the rank-1 projector is NOT representable inside the modal ansatz (the
tangential floor (1-x^2)/3 is hard-wired there), so the projector evaluation
is direct field construction, not modal coefficients.

Method: central second differences in the shear amplitude sampled at three
step sizes with a linear fit in h^2 (Delta2(h)/h^2 = a2 + (a4/12) h^2 + ...),
eliminating the leading quartic contamination; radial profile f = x^2(1-x^2);
nodes (96, 48).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
PREV = HERE.parent / "0001"
sys.path.insert(0, str(PREV))

import debox_common as base  # noqa: E402  (installs P240 0041/0042 paths, pins threads)

from cpu_energy import (  # noqa: E402
    chebyshev_stack,
    commutator,
    elementwise_derivative,
    frobenius_squared,
    gauss_grid,
)


STEPS = (2.0e-3, 4.0e-3, 8.0e-3)
NODES = (96, 48)
RADII = [12.0, 14.0, 16.0, 18.0]
ROTATION_Z = torch.tensor(
    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype=torch.float64
)


def shear_matrix(shear: str, director, polar, azimuthal):
    if shear == "np":
        return (
            director[..., :, None] * polar[..., None, :]
            + polar[..., :, None] * director[..., None, :]
        )
    if shear == "na":
        return (
            director[..., :, None] * azimuthal[..., None, :]
            + azimuthal[..., :, None] * director[..., None, :]
        )
    raise ValueError(shear)


def frames(radius: float, with_grad: bool):
    radial, radial_weight, mu, angular_weight = gauss_grid(NODES[0], NODES[1], radius)
    if with_grad:
        radius_grid = radial[:, None].repeat(1, NODES[1]).clone().requires_grad_(True)
        mu_grid = mu[None, :].repeat(NODES[0], 1).clone().requires_grad_(True)
    else:
        radius_grid = radial[:, None].repeat(1, NODES[1])
        mu_grid = mu[None, :].repeat(NODES[0], 1)
    normalized = radius_grid / radius
    sine = torch.sqrt(torch.clamp(1 - mu_grid**2, min=0.0))
    zero = torch.zeros_like(sine)
    director = torch.stack((sine, zero, mu_grid), dim=-1)
    polar = torch.stack((mu_grid, zero, -sine), dim=-1)
    azimuthal = torch.stack((zero, torch.ones_like(zero), zero), dim=-1)
    return (
        radius_grid,
        mu_grid,
        normalized,
        radial_weight,
        angular_weight,
        director,
        polar,
        azimuthal,
    )


def outer_pair(v):
    return v[..., :, None] * v[..., None, :]


def energy_density_integrate(
    spatial, static_density, radius_grid, radial_weight, angular_weight
):
    spatial_two = spatial @ spatial
    trace_two = torch.diagonal(spatial_two, dim1=-2, dim2=-1).sum(-1)
    trace_three = torch.diagonal(spatial_two @ spatial, dim1=-2, dim2=-1).sum(-1)
    potential = -0.5 * trace_two - trace_three + trace_two**2 + 0.5
    extra = static_density if static_density is not None else 0.0
    weights = (
        2 * torch.pi * radius_grid**2 * radial_weight[:, None] * angular_weight[None, :]
    )
    return float(torch.sum(weights * (potential + extra)))


def projector_energy(shear: str | None, t: float, part: str) -> float:
    """Energy at the TRUE rank-1 projector plus optional t*f shear."""
    radius_grid, mu_grid, normalized, rw, aw, director, polar, azimuthal = frames(
        18.0, with_grad=(part == "static")
    )
    spatial = outer_pair(director)
    if shear is not None and t != 0.0:
        w = shear_matrix(shear, director, polar, azimuthal)
        profile = (normalized**2 * (1 - normalized**2))[..., None, None]
        spatial = spatial + t * profile * w
    static_density = None
    if part == "static":
        static_density = static_from_spatial(spatial, radius_grid, mu_grid, rw, aw)
    return energy_density_integrate(spatial, static_density, radius_grid, rw, aw)


def root_energy(flat: np.ndarray, radius: float, shear: str | None, t: float, part: str) -> float:
    """Energy at a de-boxed root background plus optional t*f shear."""
    order = flat.size // 3
    coefficients = torch.tensor(flat.reshape(3, order), dtype=torch.float64)
    radius_grid, mu_grid, normalized, rw, aw, director, polar, azimuthal = frames(
        radius, with_grad=(part == "static")
    )
    radial_coordinate = 2 * normalized**2 - 1
    radial_basis = chebyshev_stack(radial_coordinate, tuple(range(order)))
    modal = torch.einsum("...i,ci->...c", radial_basis, coefficients)
    q = normalized**2 + normalized**2 * (1 - normalized**2) * modal[..., 0]
    tangent = (1 - normalized**2) * (
        torch.tensor(1 / 3, dtype=torch.float64) + modal[..., 1]
    )
    split_amplitude = normalized**4 * (1 - normalized**2) * modal[..., 2]
    sine = torch.sqrt(torch.clamp(1 - mu_grid**2, min=0.0))
    delta = split_amplitude * sine**2

    def outer(v):
        return v[..., :, None] * v[..., None, :]

    lambda_n = tangent + q
    spatial = (
        lambda_n[..., None, None] * outer(director)
        + (tangent + delta)[..., None, None] * outer(polar)
        + (tangent - delta)[..., None, None] * outer(azimuthal)
    )
    if shear is not None and t != 0.0:
        w = shear_matrix(shear, director, polar, azimuthal)
        profile = (normalized**2 * (1 - normalized**2))[..., None, None]
        spatial = spatial + t * profile * w
    static_density = None
    if part == "static":
        static_density = static_from_spatial(spatial, radius_grid, mu_grid, rw, aw)
    return energy_density_integrate(spatial, static_density, radius_grid, rw, aw)


def static_from_spatial(spatial, radius_grid, mu_grid, radial_weight, angular_weight):
    derivative_r = elementwise_derivative(spatial, radius_grid)
    derivative_mu = elementwise_derivative(spatial, mu_grid)
    sine = torch.sqrt(torch.clamp(1.0 - mu_grid**2, min=0.0))
    derivative_theta = -sine[..., None, None] * derivative_mu / radius_grid[..., None, None]
    derivative_phi = (
        ROTATION_Z @ spatial + spatial @ ROTATION_Z.T
    ) / (radius_grid * sine)[..., None, None]
    derivatives = (derivative_r, derivative_theta, derivative_phi)
    return 4 * sum(
        frobenius_squared(commutator(derivatives[left], derivatives[right]))
        for left in range(3)
        for right in range(left + 1, 3)
    )


def quadratic_coefficient(evaluator) -> dict:
    """Evaluator maps the amplitude t to the energy part.

    V(P + tW) = V(P) + (a2/2) t^2 + (a4/24) t^4 + ..., so the central second
    difference obeys Delta2(h)/h^2 = a2 + (a4/12) h^2 + O(h^4). Sampling
    three step sizes and fitting a line in h^2 eliminates the leading quartic
    contamination (which the control shows dominates naive two-step
    Richardson); the intercept is the quoted a2.
    """
    base_value = evaluator(0.0)
    ys = []
    for step in STEPS:
        second = evaluator(step) - 2.0 * base_value + evaluator(-step)
        ys.append(second / step**2)
    hs = np.asarray(STEPS, dtype=np.float64) ** 2
    slope, intercept = np.polyfit(hs, np.asarray(ys), 1)
    fitted = slope * hs + intercept
    return {
        "per_step": {f"h={step:g}": y for step, y in zip(STEPS, ys)},
        "quartic_slope": float(slope),
        "a2_intercept": float(intercept),
        "fit_residual": float(np.max(np.abs(np.asarray(ys) - fitted))),
    }


def main() -> None:
    started = time.time()
    with open(PREV / "clean-ladder.json") as handle:
        previous = json.load(handle)
    prior = {row["radius"]: row for row in previous["rungs"] if row.get("accepted")}

    control: dict = {}
    for shear in ("np", "na"):
        row = {}
        for part in ("potential", "static"):
            row[part] = quadratic_coefficient(
                lambda t, s=shear, p=part: projector_energy(s, t, p)
            )
            print(
                f"[control projector {shear} {part}] "
                f"a2={row[part]['a2_intercept']:+.6e} "
                f"(fit residual {row[part]['fit_residual']:.2e})",
                flush=True,
            )
        control[shear] = row
    floor_pot = max(abs(control[s]["potential"]["a2_intercept"]) for s in control)
    floor_stat = max(abs(control[s]["static"]["a2_intercept"]) for s in control)

    records: dict = {}
    for radius in RADII:
        flat = np.asarray(prior[radius]["values"], dtype=np.float64)
        row: dict = {}
        for shear in ("np", "na"):
            entry: dict = {}
            for part, floor in (("potential", floor_pot), ("static", floor_stat)):
                entry[part] = quadratic_coefficient(
                    lambda t, s=shear, f=flat, r=radius, p=part: root_energy(
                        f, r, s, t, p
                    )
                )
                a2 = entry[part]["a2_intercept"]
                print(
                    f"[R={radius} {shear} {part}] a2={a2:+.6e} "
                    f"(fit residual {entry[part]['fit_residual']:.2e}, "
                    f"floor {floor:.2e}, resolved={abs(a2) > 10 * floor})",
                    flush=True,
                )
            row[shear] = entry
        records[str(radius)] = row

    payload = {
        "control_projector": control,
        "probe_floor_potential": floor_pot,
        "probe_floor_static": floor_stat,
        "records": records,
        "nodes": list(NODES),
        "radial_profile": "x^2 (1 - x^2)",
        "minutes_total": round((time.time() - started) / 60.0, 2),
    }
    (HERE / "g7-shear-probe.json").write_text(json.dumps(payload, indent=2))
    print("WROTE g7-shear-probe.json")


if __name__ == "__main__":
    main()
