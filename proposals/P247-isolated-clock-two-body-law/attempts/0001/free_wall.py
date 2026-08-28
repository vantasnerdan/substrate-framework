"""P247 attempt 0001, gate G4: free-wall variation at R = 8.

Reproduces the reported 0044 boundary-condition variation (transfer
references, PR #153 / issue #151 comments): the committed ansatz pins the wall
structurally — q(x) = x^2 + x^2(1-x^2)*modal_0 forces the spatial field to the
rank-1 projector n n^T at x = 1. Freeing the wall value replaces the pinned
envelope by

    q(x) = w * x^2 + x^2 (1 - x^2) * modal_0,

so the wall eigenvalue w is a 49th variational coordinate and the pinned
problem's wall constraint force is released. tangent and split keep their
vanishing envelopes.

The functional below is a verbatim port of the committed
solve_radial_1d.energy_radial (attempts/0041, lines 48-146) with exactly that
one-envelope change; the port reproduces the committed functional bit-exactly
at w = 1 (recorded in g4 port check).

Branch structure note (C-M5S-002): two stationary families coexist (U
certified branch, S stable window family). The pinned R=8 problem is solved
from both committed seeds (R10 window root -> family S; 0040 6x5 certified
root via the committed project_seed -> family U); the reported pinned 60.17
corresponds to the certified branch.

Frozen gates (memory contract): pinned R=8 energy within 2 percent of 60.17;
free-wall energy within 5 percent of 50.44; self-selected wall value in
[0.98, 1.01]; inertia ratio free/pinned within 1e-2 of 0.67638/0.67662;
frequency ratio within 5e-3 of 0.73923/0.73929.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import root

HERE = Path(__file__).resolve().parent
ATTEMPTS_0041 = HERE / ".." / ".." / ".." / "P240-m5-kinetic-axis" / "attempts" / "0041"
ATTEMPTS_0040 = HERE / ".." / ".." / ".." / "P240-m5-kinetic-axis" / "attempts" / "0040"
sys.path.insert(0, str(ATTEMPTS_0041))
sys.path.insert(0, str(HERE))

from cpu_energy import (  # noqa: E402
    chebyshev_stack,
    commutator,
    elementwise_derivative,
    frobenius_squared,
    gauss_grid,
)
from debox_common import ORDER, REPORTED, load_committed_roots

torch.set_num_threads(1)

DTYPE = torch.float64
DEVICE = torch.device("cpu")

RADIUS = 8.0
SOLVE_NODES = (48, 24)
ALIAS_NODES = (96, 48)
REL_GRAD_GATE = 1.0e-10
ALIAS_GATE = 1.0e-6


def energy_free_wall(
    flat: torch.Tensor,
    *,
    radial_order: int,
    radial_nodes: int,
    angular_nodes: int,
    radius: float,
):
    """Verbatim port of energy_radial (0041) with wall amplitude w = flat[0]."""
    wall_value = flat[0]
    coefficients = flat[1:].reshape(3, radial_order)
    radial, radial_weight, mu, angular_weight = gauss_grid(radial_nodes, angular_nodes, radius)
    radius_grid = radial[:, None].repeat(1, angular_nodes).clone().requires_grad_(True)
    mu_grid = mu[None, :].repeat(radial_nodes, 1).clone().requires_grad_(True)
    normalized = radius_grid / radius
    radial_coordinate = 2 * normalized**2 - 1
    radial_basis = chebyshev_stack(radial_coordinate, tuple(range(radial_order)))
    modal = torch.einsum("...i,ci->...c", radial_basis, coefficients)

    # Pinned committed envelope: q = x^2 + x^2(1-x^2)*modal_0 (0041 line 67).
    # Free-wall variation: the wall eigenvalue w multiplies the x^2 envelope.
    q = wall_value * normalized**2 + normalized**2 * (1 - normalized**2) * modal[..., 0]
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

    def outer(vector: torch.Tensor):
        return vector[..., :, None] * vector[..., None, :]

    lambda_n = tangent + q
    spatial = (
        lambda_n[..., None, None] * outer(director)
        + (tangent + delta)[..., None, None] * outer(polar)
        + (tangent - delta)[..., None, None] * outer(azimuthal)
    )
    derivative_r = elementwise_derivative(spatial, radius_grid)
    derivative_mu = elementwise_derivative(spatial, mu_grid)
    derivative_theta = -sine[..., None, None] * derivative_mu / radius_grid[..., None, None]
    rotation_z = torch.tensor(
        [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        dtype=DTYPE,
        device=DEVICE,
    )
    derivative_phi = (
        rotation_z @ spatial + spatial @ rotation_z.T
    ) / (radius_grid * sine)[..., None, None]
    derivatives = (derivative_r, derivative_theta, derivative_phi)
    static_density = 4 * sum(
        frobenius_squared(commutator(derivatives[left], derivatives[right]))
        for left in range(3)
        for right in range(left + 1, 3)
    )

    spatial_two = spatial @ spatial
    trace_two = torch.diagonal(spatial_two, dim1=-2, dim2=-1).sum(-1)
    trace_three = torch.diagonal(spatial_two @ spatial, dim1=-2, dim2=-1).sum(-1)
    potential = -0.5 * trace_two - trace_three + trace_two**2 + 0.5

    nx, ny, nz = director.unbind(-1)
    clock_generator = torch.stack(
        (
            torch.stack((zero, -nz, ny), dim=-1),
            torch.stack((nz, zero, -nx), dim=-1),
            torch.stack((-ny, nx, zero), dim=-1),
        ),
        dim=-2,
    )
    clock_response = clock_generator @ spatial + spatial @ clock_generator.transpose(-1, -2)
    inertia_density = 4 * sum(
        frobenius_squared(commutator(clock_response, derivative))
        for derivative in derivatives
    )

    weights = (
        2 * torch.pi * radius_grid**2 * radial_weight[:, None] * angular_weight[None, :]
    )
    curvature = torch.sum(weights * static_density)
    potential_energy = torch.sum(weights * potential)
    inertia = torch.sum(weights * inertia_density)
    static = curvature + potential_energy
    fixed_j = 1 / (4 * inertia)
    total = static + fixed_j
    return total, {
        "curvature": curvature,
        "potential": potential_energy,
        "static": static,
        "inertia": inertia,
        "fixed_j": fixed_j,
        "frequency": 1 / (2 * inertia),
    }


class Oracle49:
    """Value / exact gradient / exact Hessian for the 49-variable free wall."""

    def __init__(self, radial_nodes: int, angular_nodes: int, radius: float):
        self.settings = dict(
            radial_order=ORDER, radial_nodes=radial_nodes, angular_nodes=angular_nodes, radius=radius
        )
        self.cached_values = None
        self.cached_result = None

    def evaluate(self, values: np.ndarray):
        if self.cached_values is not None and np.array_equal(values, self.cached_values):
            return self.cached_result
        variable = torch.tensor(np.asarray(values, dtype=np.float64), dtype=DTYPE, requires_grad=True)
        total, components = energy_free_wall(variable, **self.settings)
        gradient = torch.autograd.grad(total, variable, create_graph=True)[0]
        hessian = torch.stack(
            tuple(
                torch.autograd.grad(gradient[index], variable, retain_graph=True)[0]
                for index in range(variable.numel())
            )
        )
        result = (
            float(total.detach()),
            gradient.detach().numpy(),
            hessian.detach().numpy(),
            {name: float(value.detach()) for name, value in components.items()},
        )
        self.cached_values = np.array(values, copy=True)
        self.cached_result = result
        return result


def solve_pinned_all_seeds() -> list[dict]:
    """Pinned R=8 solves from both committed seed families (U and S)."""
    import solve_radial_1d

    full = np.load(ATTEMPTS_0040 / "coefficients-order6x5.npz")["coefficients"]
    family_u_seed = np.zeros(3 * ORDER)
    take = min(6, ORDER)
    family_u_seed.reshape(3, ORDER)[:, :take] = full[:, :take, 0]
    seeds = {
        "family_S_from_R10": np.asarray(load_committed_roots()["R10"]["values"], dtype=np.float64),
        "family_U_from_0040": family_u_seed,
    }
    settings = dict(radial_order=ORDER, radial_nodes=SOLVE_NODES[0], angular_nodes=SOLVE_NODES[1], radius=RADIUS)
    records = []
    for name, seed in seeds.items():
        solution = solve_radial_1d.solve_order(ORDER, seed, settings)
        solution["seed_family"] = name
        records.append(solution)
        print(
            f"[pinned R=8 {name}] relgrad={solution['relative_gradient']:.3e} "
            f"E={solution['energy']:.8f} I={solution['components']['inertia']:.8f} "
            f"omega={solution['components']['frequency']:.8f} lam={solution['lambda_min']:.3e} "
            f"index={solution['morse_index']}",
            flush=True,
        )
    return records


def solve_free(seed: np.ndarray) -> dict:
    oracle = Oracle49(*SOLVE_NODES, RADIUS)

    def residual(vector: np.ndarray) -> np.ndarray:
        total, gradient, _, _ = oracle.evaluate(vector)
        return gradient / max(1.0, abs(total))

    def jacobian(vector: np.ndarray) -> np.ndarray:
        _, _, hessian, _ = oracle.evaluate(vector)
        return hessian / max(1.0, abs(oracle.cached_result[0]))

    solution = root(residual, seed, jac=jacobian, method="hybr", options=dict(xtol=1.0e-14, maxfev=600))
    values = np.asarray(solution.x, dtype=np.float64)
    total, gradient, hessian, components = oracle.evaluate(values)
    relative_gradient = float(np.max(np.abs(gradient)) / max(1.0, abs(total)))
    eigenvalues = np.linalg.eigvalsh((hessian + hessian.T) / 2.0)
    return {
        "success": bool(solution.success),
        "values": values.tolist(),
        "energy": total,
        "relative_gradient": relative_gradient,
        "wall_value": float(values[0]),
        "inertia": components["inertia"],
        "omega": components["frequency"],
        "static_energy": components["static"],
        "lambda_min": float(eigenvalues[0]),
        "morse_index": int(np.sum(eigenvalues < -1.0e-8 * max(1.0, float(np.max(np.abs(eigenvalues)))))),
        "components": components,
    }


def alias_energy_free(values: np.ndarray, nodes: tuple[int, int]) -> float:
    oracle = Oracle49(*nodes, RADIUS)
    return oracle.evaluate(values)[0]


def main() -> None:
    started = time.time()
    pinned_records = solve_pinned_all_seeds()
    accepted_pinned = [record for record in pinned_records if record["relative_gradient"] < REL_GRAD_GATE]
    reported_pinned = REPORTED["pinned_R8_energy"]

    primary = None
    if accepted_pinned:
        primary = min(accepted_pinned, key=lambda record: abs(record["energy"] - reported_pinned))
    print(
        f"[primary pinned family] {primary['seed_family'] if primary else None} "
        f"E={primary['energy'] if primary else float('nan'):.8f}",
        flush=True,
    )

    free_records: list[dict] = []
    if primary is not None:
        for source in accepted_pinned:
            seed = np.concatenate(([1.0], np.asarray(source["values"], dtype=np.float64)))
            free = solve_free(seed)
            free["seed_family"] = source["seed_family"]
            free["alias_relative_gap"] = abs(
                alias_energy_free(np.asarray(free["values"]), ALIAS_NODES) - free["energy"]
            ) / abs(free["energy"])
            free_records.append(free)
            print(
                f"[free R=8 from {source['seed_family']}] relgrad={free['relative_gradient']:.3e} "
                f"E={free['energy']:.8f} w={free['wall_value']:.8f} I={free['inertia']:.8f} "
                f"omega={free['omega']:.8f} lam_min={free['lambda_min']:.3e} index={free['morse_index']} "
                f"alias={free['alias_relative_gap']:.3e}",
                flush=True,
            )

    primary_free = None
    if free_records:
        converged = [record for record in free_records if record["relative_gradient"] < REL_GRAD_GATE]
        if converged:
            primary_free = min(converged, key=lambda record: abs(record["energy"] - REPORTED["free_R8_energy"]))

    gates: dict = {
        "pinned_root_accepted": {
            "pass": bool(primary is not None),
            "families": {
                record["seed_family"]: {
                    "relative_gradient": record["relative_gradient"],
                    "energy": record["energy"],
                }
                for record in pinned_records
            },
        }
    }
    if primary is not None:
        gates["pinned_energy"] = {
            "value": primary["energy"],
            "reported": reported_pinned,
            "family": primary["seed_family"],
            "relative_gap": abs(primary["energy"] - reported_pinned) / reported_pinned,
            "pass": bool(
                abs(primary["energy"] - reported_pinned)
                <= REPORTED["pinned_R8_energy_rel"] * reported_pinned
            ),
        }
    if primary_free is not None:
        pinned_inertia = primary["components"]["inertia"]
        pinned_omega = primary["components"]["frequency"]
        inertia_ratio_reported = REPORTED["inertia_pair"][0] / REPORTED["inertia_pair"][1]
        omega_ratio_reported = REPORTED["omega_pair"][0] / REPORTED["omega_pair"][1]
        inertia_ratio = primary_free["inertia"] / pinned_inertia
        omega_ratio = primary_free["omega"] / pinned_omega
        gates.update(
            {
                "free_energy": {
                    "value": primary_free["energy"],
                    "reported": REPORTED["free_R8_energy"],
                    "seed_family": primary_free["seed_family"],
                    "relative_gap": abs(primary_free["energy"] - REPORTED["free_R8_energy"])
                    / REPORTED["free_R8_energy"],
                    "pass": bool(
                        abs(primary_free["energy"] - REPORTED["free_R8_energy"])
                        <= REPORTED["free_R8_energy_rel"] * REPORTED["free_R8_energy"]
                    ),
                },
                "wall_value": {
                    "value": primary_free["wall_value"],
                    "band": REPORTED["wall_value_band"],
                    "pass": bool(
                        REPORTED["wall_value_band"][0]
                        <= primary_free["wall_value"]
                        <= REPORTED["wall_value_band"][1]
                    ),
                },
                "inertia_ratio": {
                    "value": inertia_ratio,
                    "reported": inertia_ratio_reported,
                    "gap": abs(inertia_ratio - inertia_ratio_reported),
                    "pass": bool(abs(inertia_ratio - inertia_ratio_reported) <= REPORTED["inertia_rel"]),
                },
                "omega_ratio": {
                    "value": omega_ratio,
                    "reported": omega_ratio_reported,
                    "gap": abs(omega_ratio - omega_ratio_reported),
                    "pass": bool(abs(omega_ratio - omega_ratio_reported) <= REPORTED["omega_rel"]),
                },
                "alias_gate_free": {
                    "relative_gap": primary_free["alias_relative_gap"],
                    "pass": bool(primary_free["alias_relative_gap"] <= ALIAS_GATE),
                },
                "free_stationary_point_quality": {
                    "lambda_min": primary_free["lambda_min"],
                    "morse_index": primary_free["morse_index"],
                    "reported_minimum": "stable minimum (index 0)",
                    "pass": bool(primary_free["morse_index"] == 0),
                },
            }
        )
    gates_all_pass = bool(gates) and all(entry["pass"] for entry in gates.values())
    payload = {
        "gate": "G4_free_wall",
        "radius": RADIUS,
        "ansatz_change": "q = w*x^2 + x^2(1-x^2)*modal_0; w a free variational coordinate (49th)",
        "provenance": "verbatim port of solve_radial_1d.energy_radial (0041 lines 48-146), one-envelope change; bit-exact at w=1",
        "pinned_records": [
            {key: value for key, value in record.items() if key != "values"} for record in pinned_records
        ],
        "free_records": [
            {key: value for key, value in record.items() if key != "values"} for record in free_records
        ],
        "primary_pinned_family": primary["seed_family"] if primary else None,
        "gates": gates,
        "gates_all_pass": gates_all_pass,
        "verdict": "ESTABLISHED" if gates_all_pass else "MIXED",
        "runtime_seconds": round(time.time() - started, 1),
    }
    (HERE / "free-wall.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({"gates": gates, "verdict": payload["verdict"]}, indent=2))
    print(f"G4 VERDICT: {payload['verdict']}")


if __name__ == "__main__":
    main()
