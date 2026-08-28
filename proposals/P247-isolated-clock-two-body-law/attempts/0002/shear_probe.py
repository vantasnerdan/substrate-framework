"""P247 attempt 0002: potential-curvature probe along the G1 flat shear channels.

The G1 construction proves the np and na shear channels are exactly flat in
the POTENTIAL part at the rank-1 projector background n-hat n-hat^T. This
probe evaluates the potential-part quadratic coefficient along those shear
directions at (i) the projector background - analytically zero, so its
numeric value is the probe's own roundoff and quadrature floor - and (ii)
each de-boxed root R=12..18, recording where the exactly-flat channels land
on the actual branch backgrounds (candidate A's long-range sector health).

Note: the rank-1 projector is NOT representable inside the modal ansatz (the
tangential floor (1-x^2)/3 is hard-wired there), so the control is built by
direct field construction, not through modal coefficients.

Method: central second differences in the shear amplitude with Richardson in
the step (fine 1e-4, coarse 1e-3); radial profile f = x^2 (1 - x^2); nodes
(96, 48); potential part only, matching the exact G1 statement.
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

from cpu_energy import gauss_grid  # noqa: E402

torch.set_grad_enabled(False)

STEPS = (2.0e-3, 4.0e-3, 8.0e-3)
NODES = (96, 48)
RADII = [12.0, 14.0, 16.0, 18.0]


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


def integrate_potential(spatial, radius_grid, radial_weight, angular_weight):
    spatial_two = spatial @ spatial
    trace_two = torch.diagonal(spatial_two, dim1=-2, dim2=-1).sum(-1)
    trace_three = torch.diagonal(spatial_two @ spatial, dim1=-2, dim2=-1).sum(-1)
    potential = -0.5 * trace_two - trace_three + trace_two**2 + 0.5
    weights = (
        2 * torch.pi * radius_grid**2 * radial_weight[:, None] * angular_weight[None, :]
    )
    return float(torch.sum(weights * potential))


def build_frames(radius: float):
    radial, radial_weight, mu, angular_weight = gauss_grid(NODES[0], NODES[1], radius)
    radius_grid = radial[:, None].repeat(1, NODES[1])
    mu_grid = mu[None, :].repeat(NODES[0], 1)
    sine = torch.sqrt(torch.clamp(1 - mu_grid**2, min=0.0))
    zero = torch.zeros_like(sine)
    director = torch.stack((sine, zero, mu_grid), dim=-1)
    polar = torch.stack((mu_grid, zero, -sine), dim=-1)
    azimuthal = torch.stack((zero, torch.ones_like(zero), zero), dim=-1)
    normalized = radius_grid / radius
    return (
        radius_grid,
        normalized,
        radial_weight,
        angular_weight,
        director,
        polar,
        azimuthal,
    )


def projector_potential_energy(shear: str | None, t: float) -> float:
    """Potential part at the TRUE rank-1 projector plus optional t*f shear."""
    frames = build_frames(18.0)
    radius_grid, normalized, radial_weight, angular_weight = frames[:4]
    director, polar, azimuthal = frames[4:]

    def outer(v):
        return v[..., :, None] * v[..., None, :]

    spatial = outer(director)
    if shear is not None and t != 0.0:
        w = shear_matrix(shear, director, polar, azimuthal)
        profile = (normalized**2 * (1 - normalized**2))[..., None, None]
        spatial = spatial + t * profile * w
    return integrate_potential(spatial, radius_grid, radial_weight, angular_weight)


def root_potential_energy(flat: np.ndarray, radius: float, shear: str | None, t: float) -> float:
    order = flat.size // 3
    coefficients = torch.tensor(flat.reshape(3, order), dtype=torch.float64)
    frames = build_frames(radius)
    radius_grid, normalized, radial_weight, angular_weight = frames[:4]
    director, polar, azimuthal = frames[4:]
    radial_coordinate = 2 * normalized**2 - 1
    angle = torch.acos(torch.clamp(radial_coordinate, -1.0, 1.0))
    radial_basis = torch.stack(tuple(torch.cos(d * angle) for d in range(order)), dim=-1)
    modal = torch.einsum("...i,ci->...c", radial_basis, coefficients)
    q = normalized**2 + normalized**2 * (1 - normalized**2) * modal[..., 0]
    tangent = (1 - normalized**2) * (
        torch.tensor(1 / 3, dtype=torch.float64) + modal[..., 1]
    )
    split_amplitude = normalized**4 * (1 - normalized**2) * modal[..., 2]
    sine_sq = director[..., 0] ** 2
    delta = split_amplitude * sine_sq

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
    return integrate_potential(spatial, radius_grid, radial_weight, angular_weight)


def quadratic_coefficient(evaluator) -> dict:
    """Evaluator maps the amplitude t to the potential part.

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
        control[shear] = quadratic_coefficient(
            lambda t, shear=shear: projector_potential_energy(shear, t)
        )
        print(
            f"[control projector {shear}] a2={control[shear]['a2_intercept']:+.6e} "
            f"(fit residual {control[shear]['fit_residual']:.2e})",
            flush=True,
        )
    floor = max(abs(control[s]["a2_intercept"]) for s in control)

    records: dict = {}
    for radius in RADII:
        flat = np.asarray(prior[radius]["values"], dtype=np.float64)
        row: dict = {}
        for shear in ("np", "na"):
            row[shear] = quadratic_coefficient(
                lambda t, shear=shear, flat=flat, radius=radius: root_potential_energy(
                    flat, radius, shear, t
                )
            )
            print(
                f"[R={radius} {shear}] a2={row[shear]['a2_intercept']:+.6e} "
                f"(fit residual {row[shear]['fit_residual']:.2e}, floor {floor:.2e}, "
                f"resolved={abs(row[shear]['a2_intercept']) > 10 * floor})",
                flush=True,
            )
        records[str(radius)] = row

    payload = {
        "control_projector": control,
        "probe_floor": floor,
        "records": records,
        "nodes": list(NODES),
        "radial_profile": "x^2 (1 - x^2)",
        "minutes_total": round((time.time() - started) / 60.0, 2),
    }
    (HERE / "g7-shear-probe.json").write_text(json.dumps(payload, indent=2))
    print("WROTE g7-shear-probe.json")


if __name__ == "__main__":
    main()
