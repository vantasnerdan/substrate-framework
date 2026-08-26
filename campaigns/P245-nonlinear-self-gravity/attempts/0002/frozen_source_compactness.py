"""P245 attempt 0002: repaired nonlinear compactness and homothetic gate.

Units, signs, and geometric factors are frozen in the campaign manifest:
``c=1``, ``m'=4*pi*r**2*rho``, ``f=1-2*G*m/r``.  The source is exactly the
C-M5S-004 static curvature-plus-potential density; the global fixed-J term is
not localized.  Cell masses use exact spherical shell volumes, so no small
quantity is formed by subtracting large energies.
"""

from __future__ import annotations

import os

for _name in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

import json
import sys
from pathlib import Path

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "proposals/P240-m5-kinetic-axis/attempts/0041"))

from cpu_energy import (  # noqa: E402
    chebyshev_stack,
    commutator,
    elementwise_derivative,
    frobenius_squared,
    gauss_grid,
)
from solve_radial_1d import Oracle  # noqa: E402
from substrate_framework.nonlinear_clock_gravity import (  # noqa: E402
    homothetic_compactness_profile,
    integrate_spherical_density_cells,
    minimize_homothetic_max_compactness,
)
from substrate_framework.verification import CheckLedger  # noqa: E402

DTYPE = torch.float64
G_TOTAL = 46.80699908016004
RADIUS = 12.0
ACCEPTED_MASS = 54.70900884959007
RADIAL_CELLS = (400, 800, 1600)
ANGULAR_NODES = 64
SCALE_BOUNDS = (0.12, 1200.0)
MUTATION_FACTOR = 1.0e-6

torch.set_num_threads(1)


def sector_components(
    values: np.ndarray,
    radii: np.ndarray,
    *,
    angular_nodes: int = ANGULAR_NODES,
) -> tuple[np.ndarray, np.ndarray]:
    """Return C-M5S-004 curvature and potential density separately."""

    order = len(values) // 3
    coefficients = torch.tensor(values.reshape(3, order), dtype=DTYPE)
    _, _, mu_x, mu_w = gauss_grid(8, angular_nodes, RADIUS)
    weights = mu_w.clone().detach()
    curvature_out = np.empty(len(radii), dtype=np.float64)
    potential_out = np.empty(len(radii), dtype=np.float64)

    for index, radius_value in enumerate(radii):
        radius_grid = torch.full(
            (angular_nodes,), float(radius_value), dtype=DTYPE, requires_grad=True
        )
        mu_grid = mu_x.clone().detach().requires_grad_(True)
        normalized = radius_grid / RADIUS
        radial_coordinate = 2.0 * normalized**2 - 1.0
        radial_basis = chebyshev_stack(radial_coordinate, tuple(range(order)))
        modal = torch.einsum("...i,ci->...c", radial_basis, coefficients)
        q = normalized**2 * (1.0 + (1.0 - normalized**2) * modal[..., 0])
        tangent = (1.0 - normalized**2) * (
            torch.tensor(1.0 / 3.0, dtype=DTYPE) + modal[..., 1]
        )
        split_amplitude = normalized**4 * (1.0 - normalized**2) * modal[..., 2]
        sine = torch.sqrt(torch.clamp(1.0 - mu_grid**2, min=0.0))
        delta = split_amplitude * sine**2
        zero = torch.zeros_like(sine)
        director = torch.stack((sine, zero, mu_grid), dim=-1)
        polar = torch.stack((mu_grid, zero, -sine), dim=-1)
        azimuthal = torch.stack((zero, torch.ones_like(zero), zero), dim=-1)

        def outer(vector: torch.Tensor) -> torch.Tensor:
            return vector[..., :, None] * vector[..., None, :]

        lambda_n = tangent + q
        spatial = (
            lambda_n[..., None, None] * outer(director)
            + (tangent + delta)[..., None, None] * outer(polar)
            + (tangent - delta)[..., None, None] * outer(azimuthal)
        )
        derivative_r = elementwise_derivative(spatial, radius_grid)
        derivative_mu = elementwise_derivative(spatial, mu_grid)
        derivative_theta = (
            -sine[..., None, None] * derivative_mu / radius_grid[..., None, None]
        )
        rotation_z = torch.tensor(
            [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            dtype=DTYPE,
        )
        derivative_phi = (
            rotation_z @ spatial + spatial @ rotation_z.T
        ) / (radius_grid * sine)[..., None, None]
        derivatives = (derivative_r, derivative_theta, derivative_phi)
        curvature = 4.0 * sum(
            frobenius_squared(commutator(derivatives[left], derivatives[right]))
            for left in range(3)
            for right in range(left + 1, 3)
        )
        spatial_two = spatial @ spatial
        trace_two = torch.diagonal(spatial_two, dim1=-2, dim2=-1).sum(-1)
        trace_three = torch.diagonal(
            spatial_two @ spatial, dim1=-2, dim2=-1
        ).sum(-1)
        potential = -0.5 * trace_two - trace_three + trace_two**2 + 0.5
        curvature_out[index] = float(torch.sum(weights * curvature)) / 2.0
        potential_out[index] = float(torch.sum(weights * potential)) / 2.0

    return curvature_out, potential_out


def evaluate_grid(values: np.ndarray, cells: int) -> dict[str, object]:
    edges = np.linspace(0.0, RADIUS, cells + 1, dtype=np.float64)
    centers = 0.5 * (edges[:-1] + edges[1:])
    curvature, potential_raw = sector_components(values, centers)
    potential_floor = float(np.min(potential_raw))
    potential = np.maximum(potential_raw, 0.0)
    total = curvature + potential

    profile = integrate_spherical_density_cells(edges, total, G_TOTAL)
    curvature_profile = integrate_spherical_density_cells(edges, curvature, G_TOTAL)
    potential_profile = integrate_spherical_density_cells(edges, potential, G_TOTAL)
    crossing = profile.first_trapped_surface_radius()
    crossing_index = None
    if crossing is not None:
        crossing_index = int(np.searchsorted(edges, crossing))

    dimensionless = edges / RADIUS
    homothetic = minimize_homothetic_max_compactness(
        dimensionless,
        curvature_profile.enclosed_mass,
        potential_profile.enclosed_mass,
        reference_radius=RADIUS,
        newton_constant=G_TOTAL,
        scale_bounds=SCALE_BOUNDS,
    )
    compactness_at_reference = homothetic_compactness_profile(
        dimensionless,
        curvature_profile.enclosed_mass,
        potential_profile.enclosed_mass,
        reference_radius=RADIUS,
        scale_radius=RADIUS,
        newton_constant=G_TOTAL,
    )
    mutated = integrate_spherical_density_cells(
        edges, total, G_TOTAL * MUTATION_FACTOR
    )

    return {
        "cells": cells,
        "mass": profile.total_mass,
        "curvature_mass": curvature_profile.total_mass,
        "potential_mass": potential_profile.total_mass,
        "maximum_compactness": profile.maximum_compactness,
        "minimum_f": profile.minimum_radial_metric_function,
        "horizon_radius": crossing,
        "crossing_index": crossing_index,
        "critical_newton_constant": profile.critical_newton_constant,
        "exterior_horizon_radius": profile.exterior_horizon_radius,
        "potential_roundoff_floor": potential_floor,
        "reference_compactness_agreement": float(
            np.max(np.abs(compactness_at_reference - profile.compactness))
        ),
        "homothetic": {
            "scale_radius": homothetic.scale_radius,
            "minimum_maximum_compactness": homothetic.maximum_compactness,
            "critical_newton_constant": homothetic.critical_newton_constant,
            "lower_scale_compactness": homothetic.lower_scale_compactness,
            "upper_scale_compactness": homothetic.upper_scale_compactness,
            "success": homothetic.optimizer_success,
            "evaluations": homothetic.function_evaluations,
        },
        "weak_coupling_mutation": {
            "factor": MUTATION_FACTOR,
            "minimum_f": mutated.minimum_radial_metric_function,
            "crossing": mutated.first_trapped_surface_radius(),
        },
    }


def relative_difference(left: float, right: float) -> float:
    return abs(left - right) / max(abs(right), np.finfo(np.float64).tiny)


def main() -> int:
    ledger = CheckLedger("P245-attempt-0002-nonlinear-compactness")
    roots = json.loads(
        (
            REPO
            / "proposals/P240-m5-kinetic-axis/attempts/0042/largeR-roots.json"
        ).read_text()
    )
    record = roots["R12"]
    values = np.asarray(record["values"], dtype=np.float64)
    oracle = Oracle(
        dict(radial_nodes=32, angular_nodes=16, radius=RADIUS, radial_order=16)
    )
    total_energy, _gradient, _hessian, components = oracle.evaluate(values)
    static_energy = float(components["static"])
    ledger.check(
        "root_energy_transfer",
        abs(total_energy - float(record["energy"])) < 1.0e-9,
        f"recomputed={total_energy:.12f} recorded={float(record['energy']):.12f}",
    )

    grids = [evaluate_grid(values, cells) for cells in RADIAL_CELLS]
    previous, finest = grids[-2], grids[-1]
    mass_relative = relative_difference(float(finest["mass"]), ACCEPTED_MASS)
    mass_refinement = relative_difference(
        float(finest["mass"]), float(previous["mass"])
    )
    horizon_refinement = relative_difference(
        float(finest["horizon_radius"]), float(previous["horizon_radius"])
    )
    compactness_refinement = relative_difference(
        float(finest["maximum_compactness"]),
        float(previous["maximum_compactness"]),
    )
    f_error = abs(float(finest["minimum_f"]) - float(previous["minimum_f"]))
    f_margin = abs(float(finest["minimum_f"])) / max(
        f_error, np.finfo(np.float64).eps * abs(float(finest["minimum_f"]))
    )
    finest_hom = finest["homothetic"]
    previous_hom = previous["homothetic"]
    hom_refinement = relative_difference(
        float(finest_hom["minimum_maximum_compactness"]),
        float(previous_hom["minimum_maximum_compactness"]),
    )
    hom_error = abs(
        float(finest_hom["minimum_maximum_compactness"])
        - float(previous_hom["minimum_maximum_compactness"])
    )
    hom_margin = (
        float(finest_hom["minimum_maximum_compactness"]) - 1.0
    ) / max(
        hom_error,
        np.finfo(np.float64).eps
        * float(finest_hom["minimum_maximum_compactness"]),
    )
    crossing_index = int(finest["crossing_index"])
    boundary_separation = min(crossing_index, RADIAL_CELLS[-1] - crossing_index)

    ledger.check(
        "source_mass_identity",
        mass_relative < 2.0e-4 and abs(float(finest["mass"]) - static_energy) / static_energy < 2.0e-4,
        f"mass={float(finest['mass']):.12f} accepted={ACCEPTED_MASS:.12f} "
        f"static={static_energy:.12f} rel={mass_relative:.3e}",
    )
    ledger.check(
        "density_nonnegative_up_to_roundoff",
        min(float(grid["potential_roundoff_floor"]) for grid in grids) > -1.0e-10,
        f"minimum raw potential={min(float(grid['potential_roundoff_floor']) for grid in grids):.3e}",
    )
    ledger.check(
        "radial_refinement",
        mass_refinement < 2.0e-4
        and horizon_refinement < 2.0e-3
        and compactness_refinement < 2.0e-4,
        f"mass={mass_refinement:.3e} horizon={horizon_refinement:.3e} "
        f"compactness={compactness_refinement:.3e}",
    )
    ledger.check(
        "trapped_surface_crossing",
        finest["horizon_radius"] is not None
        and boundary_separation >= 10
        and f_margin >= 100.0,
        f"r_h={float(finest['horizon_radius']):.9f} min_f={float(finest['minimum_f']):.6e} "
        f"margin={f_margin:.3e} boundary_cells={boundary_separation}",
    )
    ledger.check(
        "schwarzschild_exterior_contains_source",
        float(finest["exterior_horizon_radius"]) > RADIUS
        and float(finest["maximum_compactness"]) > 1.0,
        f"2GM={float(finest['exterior_horizon_radius']):.9f} R={RADIUS:.1f} "
        f"maxC={float(finest['maximum_compactness']):.6e}",
    )
    ledger.check(
        "homothetic_scale_escape_refuted",
        bool(finest_hom["success"])
        and hom_refinement < 2.0e-3
        and hom_margin >= 100.0
        and float(finest_hom["minimum_maximum_compactness"]) > 1.0
        and float(finest_hom["lower_scale_compactness"])
        > float(finest_hom["minimum_maximum_compactness"])
        and float(finest_hom["upper_scale_compactness"])
        > float(finest_hom["minimum_maximum_compactness"]),
        f"R*={float(finest_hom['scale_radius']):.9f} "
        f"min maxC={float(finest_hom['minimum_maximum_compactness']):.6e} "
        f"refine={hom_refinement:.3e} margin={hom_margin:.3e}",
    )
    reference_budget = (
        16.0
        * np.finfo(np.float64).eps
        * max(1.0, float(finest["maximum_compactness"]))
    )
    ledger.check(
        "reference_scale_identity",
        float(finest["reference_compactness_agreement"]) <= reference_budget,
        f"max absolute difference={float(finest['reference_compactness_agreement']):.3e} "
        f"budget={reference_budget:.3e}",
    )
    mutation = finest["weak_coupling_mutation"]
    ledger.check(
        "coupling_mutation_is_horizonless",
        mutation["crossing"] is None and float(mutation["minimum_f"]) > 0.9,
        f"factor={MUTATION_FACTOR:.1e} min_f={float(mutation['minimum_f']):.9f}",
    )
    try:
        integrate_spherical_density_cells(
            np.array([0.0, 1.0]), np.array([1.0]), -G_TOTAL
        )
        negative_rejected = False
    except ValueError:
        negative_rejected = True
    ledger.check("wrong_gravity_sign_rejected", negative_rejected)
    vacuum = integrate_spherical_density_cells(
        np.array([0.0, 1.0]), np.array([0.0]), G_TOTAL
    )
    ledger.check(
        "vacuum_limit",
        vacuum.total_mass == 0.0
        and vacuum.first_trapped_surface_radius() is None
        and np.array_equal(vacuum.radial_metric_function, np.ones(2)),
    )

    payload = {
        "campaign": "P245",
        "attempt": "0002",
        "verdict": "FROZEN_AND_HOMOTHETIC_TRAPPED_SURFACE_ESTABLISHED",
        "scope": (
            "The accepted frozen source and its exact curvature-plus-potential "
            "homothetic family cannot be globally static and horizonless at "
            "the B=0 coupling. Arbitrary covariant M5 profile relaxation is "
            "not ruled out by this attempt."
        ),
        "source": {
            "root": "P240 R12 order-16 family-S",
            "g_total": G_TOTAL,
            "radius": RADIUS,
            "accepted_mass": ACCEPTED_MASS,
            "oracle_total_energy": total_energy,
            "oracle_static_energy": static_energy,
        },
        "grids": grids,
        "error_budget": {
            "mass_relative_to_accepted": mass_relative,
            "mass_refinement_relative": mass_refinement,
            "horizon_refinement_relative": horizon_refinement,
            "maximum_compactness_refinement_relative": compactness_refinement,
            "minimum_f_absolute_refinement": f_error,
            "minimum_f_signal_to_error": f_margin,
            "homothetic_minimum_refinement_relative": hom_refinement,
            "homothetic_minimum_absolute_refinement": hom_error,
            "homothetic_signal_to_error": hom_margin,
        },
    }
    output = HERE / "compactness-verdict.json"
    output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    print(f"[DONE] {output.name} written")
    ledger.finish()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
