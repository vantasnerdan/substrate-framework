"""P245 attempt 0003: independent NumPy/DOP853 compactness route.

This route imports no P240 evaluator.  It differentiates the frozen Chebyshev
profiles and rotating angular frame analytically, reconstructs the density on
nonuniform nodes, integrates the Misner--Sharp mass ODE with DOP853, and finds
the first ``f=1-2Gm/r=0`` crossing with Brent's method.
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
from pathlib import Path

import numpy as np
from numpy.polynomial.chebyshev import chebder, chebval
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline, PchipInterpolator
from scipy.optimize import brentq, minimize_scalar

from substrate_framework.verification import CheckLedger

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
G_TOTAL = 46.80699908016004
RADIUS = 12.0
SOURCE_NODES = 801
ANGULAR_NODES = 64
DENSE_SAMPLES = 4001
SCALE_BOUNDS = (0.12, 1200.0)


def outer(vector: np.ndarray) -> np.ndarray:
    return np.einsum("...i,...j->...ij", vector, vector)


def outer_derivative(vector: np.ndarray, derivative: np.ndarray) -> np.ndarray:
    return np.einsum("...i,...j->...ij", derivative, vector) + np.einsum(
        "...i,...j->...ij", vector, derivative
    )


def analytic_sector_components(
    coefficients_flat: np.ndarray,
    radii: np.ndarray,
    angular_nodes: int = ANGULAR_NODES,
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate the P240 source with analytic radial/angular derivatives."""

    coefficients = np.asarray(coefficients_flat, dtype=np.float64).reshape(3, -1)
    mu, weights = leggauss(angular_nodes)
    sine = np.sqrt(1.0 - mu**2)
    zero = np.zeros_like(mu)
    director = np.stack((sine, zero, mu), axis=-1)
    polar = np.stack((mu, zero, -sine), axis=-1)
    azimuthal = np.stack((zero, np.ones_like(mu), zero), axis=-1)
    director_mu = np.stack((-mu / sine, zero, np.ones_like(mu)), axis=-1)
    polar_mu = np.stack((np.ones_like(mu), zero, mu / sine), axis=-1)
    nn = outer(director)
    pp = outer(polar)
    aa = outer(azimuthal)
    nn_mu = outer_derivative(director, director_mu)
    pp_mu = outer_derivative(polar, polar_mu)
    rotation_z = np.array(
        [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        dtype=np.float64,
    )

    curvature_result = np.empty(len(radii), dtype=np.float64)
    potential_result = np.empty(len(radii), dtype=np.float64)
    derivative_coefficients = np.stack(
        [chebder(row) for row in coefficients], axis=0
    )

    for index, radius in enumerate(np.asarray(radii, dtype=np.float64)):
        x = radius / RADIUS
        z = 2.0 * x**2 - 1.0
        dz_dr = 4.0 * x / RADIUS
        modal = np.array([chebval(z, row) for row in coefficients])
        modal_r = np.array(
            [chebval(z, row) for row in derivative_coefficients]
        ) * dz_dr

        q = x**2 * (1.0 + (1.0 - x**2) * modal[0])
        q_r = (
            2.0 * x / RADIUS * (1.0 + (1.0 - x**2) * modal[0])
            + x**2
            * (
                -2.0 * x / RADIUS * modal[0]
                + (1.0 - x**2) * modal_r[0]
            )
        )
        tangent = (1.0 - x**2) * (1.0 / 3.0 + modal[1])
        tangent_r = (
            -2.0 * x / RADIUS * (1.0 / 3.0 + modal[1])
            + (1.0 - x**2) * modal_r[1]
        )
        split = x**4 * (1.0 - x**2) * modal[2]
        split_r = (
            (4.0 * x**3 * (1.0 - x**2) - 2.0 * x**5)
            * modal[2]
            / RADIUS
            + x**4 * (1.0 - x**2) * modal_r[2]
        )
        delta = split * sine**2
        delta_r = split_r * sine**2
        delta_mu = -2.0 * split * mu
        lambda_n = tangent + q
        lambda_n_r = tangent_r + q_r

        spatial = (
            lambda_n * nn
            + (tangent + delta)[:, None, None] * pp
            + (tangent - delta)[:, None, None] * aa
        )
        derivative_r = (
            lambda_n_r * nn
            + (tangent_r + delta_r)[:, None, None] * pp
            + (tangent_r - delta_r)[:, None, None] * aa
        )
        derivative_mu = (
            lambda_n * nn_mu
            + (tangent + delta)[:, None, None] * pp_mu
            + delta_mu[:, None, None] * (pp - aa)
        )
        derivative_theta = -sine[:, None, None] * derivative_mu / radius
        derivative_phi = (
            rotation_z @ spatial + spatial @ rotation_z.T
        ) / (radius * sine)[:, None, None]
        derivatives = (derivative_r, derivative_theta, derivative_phi)
        curvature = np.zeros(angular_nodes, dtype=np.float64)
        for left in range(3):
            for right in range(left + 1, 3):
                commutator = (
                    derivatives[left] @ derivatives[right]
                    - derivatives[right] @ derivatives[left]
                )
                curvature += 4.0 * np.einsum(
                    "...ij,...ij->...", commutator, commutator
                )
        spatial_two = spatial @ spatial
        trace_two = np.trace(spatial_two, axis1=-2, axis2=-1)
        trace_three = np.trace(spatial_two @ spatial, axis1=-2, axis2=-1)
        potential = -0.5 * trace_two - trace_three + trace_two**2 + 0.5
        curvature_result[index] = float(np.dot(weights, curvature) / 2.0)
        potential_result[index] = float(np.dot(weights, potential) / 2.0)

    return curvature_result, potential_result


def maximum_compactness(
    dense_solution: object, sample_radii: np.ndarray
) -> tuple[float, float]:
    masses = np.sum(dense_solution.sol(sample_radii), axis=0)
    compactness = 2.0 * G_TOTAL * masses / sample_radii
    index = int(np.argmax(compactness))
    left = sample_radii[max(index - 1, 0)]
    right = sample_radii[min(index + 1, sample_radii.size - 1)]

    def negative(radius: float) -> float:
        mass = float(np.sum(dense_solution.sol(radius)))
        return -2.0 * G_TOTAL * mass / radius

    refined = minimize_scalar(
        negative,
        bounds=(float(left), float(right)),
        method="bounded",
        options={"xatol": 1.0e-12},
    )
    return float(-refined.fun), float(refined.x)


def main() -> int:
    ledger = CheckLedger("P245-attempt-0003-independent-compactness")
    root_record = json.loads(
        (
            REPO
            / "proposals/P240-m5-kinetic-axis/attempts/0042/largeR-roots.json"
        ).read_text()
    )["R12"]
    coefficients = np.asarray(root_record["values"], dtype=np.float64)
    primary = json.loads(
        (
            REPO
            / "proposals/P245-nonlinear-self-gravity/attempts/0002/compactness-verdict.json"
        ).read_text()
    )
    primary_finest = primary["grids"][-1]

    angles = np.linspace(np.pi, 0.0, SOURCE_NODES)
    source_radii = 0.5 * RADIUS * (1.0 + np.cos(angles))
    evaluator_radii = source_radii.copy()
    evaluator_radii[0] = 1.0e-8
    curvature, potential_raw = analytic_sector_components(
        coefficients, evaluator_radii
    )
    potential = np.maximum(potential_raw, 0.0)
    curvature_interpolant = PchipInterpolator(source_radii, curvature)
    potential_interpolant = PchipInterpolator(source_radii, potential)

    def mass_rhs(radius: float, _state: np.ndarray) -> np.ndarray:
        factor = 4.0 * np.pi * radius**2
        return factor * np.array(
            [curvature_interpolant(radius), potential_interpolant(radius)],
            dtype=np.float64,
        )

    mass_solution = solve_ivp(
        mass_rhs,
        (0.0, RADIUS),
        np.zeros(2, dtype=np.float64),
        method="DOP853",
        rtol=1.0e-11,
        atol=1.0e-12,
        dense_output=True,
    )
    ledger.check(
        "DOP853_completed",
        mass_solution.success and mass_solution.sol is not None,
        f"status={mass_solution.status} nfev={mass_solution.nfev} message={mass_solution.message}",
    )

    sample_radii = np.linspace(RADIUS / (DENSE_SAMPLES - 1), RADIUS, DENSE_SAMPLES - 1)
    component_masses = mass_solution.sol(sample_radii)
    total_masses = np.sum(component_masses, axis=0)
    metric = 1.0 - 2.0 * G_TOTAL * total_masses / sample_radii
    crossing_positions = np.flatnonzero((metric[:-1] > 0.0) & (metric[1:] < 0.0))
    ledger.check(
        "crossing_bracket_exists",
        crossing_positions.size > 0,
        f"minimum sampled f={float(np.min(metric)):.6e}",
    )
    crossing_index = int(crossing_positions[0])

    def metric_function(radius: float) -> float:
        mass = float(np.sum(mass_solution.sol(radius)))
        return 1.0 - 2.0 * G_TOTAL * mass / radius

    horizon = brentq(
        metric_function,
        float(sample_radii[crossing_index]),
        float(sample_radii[crossing_index + 1]),
        xtol=1.0e-13,
        rtol=1.0e-13,
    )
    maximum, maximum_radius = maximum_compactness(mass_solution, sample_radii)
    terminal_components = mass_solution.y[:, -1]
    total_mass = float(np.sum(terminal_components))

    spline_curvature = CubicSpline(sample_radii, component_masses[0])
    spline_potential = CubicSpline(sample_radii, component_masses[1])
    residual_radii = sample_radii[4:-4:17]
    numerical_derivative = np.vstack(
        [spline_curvature(residual_radii, 1), spline_potential(residual_radii, 1)]
    )
    expected_derivative = np.column_stack(
        [mass_rhs(float(radius), np.zeros(2)) for radius in residual_radii]
    )
    residual = float(
        np.max(
            np.abs(numerical_derivative - expected_derivative)
            / (1.0 + np.abs(expected_derivative))
        )
    )

    x = np.concatenate(([0.0], sample_radii / RADIUS))
    curvature_mass = np.concatenate(([0.0], component_masses[0]))
    potential_mass = np.concatenate(([0.0], component_masses[1]))

    def homothetic_objective(log_scale: float) -> float:
        scale = float(np.exp(log_scale))
        scaled_mass = (
            curvature_mass * (RADIUS / scale)
            + potential_mass * (scale / RADIUS) ** 3
        )
        compactness = np.zeros_like(x)
        compactness[1:] = 2.0 * G_TOTAL * scaled_mass[1:] / (scale * x[1:])
        return float(np.max(compactness))

    homothetic = minimize_scalar(
        homothetic_objective,
        bounds=tuple(np.log(SCALE_BOUNDS)),
        method="bounded",
        options={"xatol": 1.0e-10, "maxiter": 1000},
    )
    homothetic_scale = float(np.exp(homothetic.x))
    homothetic_minimum = float(homothetic.fun)

    primary_homothetic = primary_finest["homothetic"]

    def relative(value: float, reference: float) -> float:
        return abs(value - reference) / max(abs(reference), np.finfo(float).tiny)

    differences = {
        "total_mass": relative(total_mass, float(primary_finest["mass"])),
        "curvature_mass": relative(
            float(terminal_components[0]), float(primary_finest["curvature_mass"])
        ),
        "potential_mass": relative(
            float(terminal_components[1]), float(primary_finest["potential_mass"])
        ),
        "horizon_radius": relative(horizon, float(primary_finest["horizon_radius"])),
        "maximum_compactness": relative(
            maximum, float(primary_finest["maximum_compactness"])
        ),
        "homothetic_minimum": relative(
            homothetic_minimum,
            float(primary_homothetic["minimum_maximum_compactness"]),
        ),
        "homothetic_scale": relative(
            homothetic_scale, float(primary_homothetic["scale_radius"])
        ),
    }
    trapped_cross_route_margin = (min(maximum, homothetic_minimum) - 1.0) / max(
        abs(maximum - float(primary_finest["maximum_compactness"])),
        abs(
            homothetic_minimum
            - float(primary_homothetic["minimum_maximum_compactness"])
        ),
        np.finfo(float).eps * maximum,
    )

    ledger.check(
        "source_and_mass_cross_route_agreement",
        differences["total_mass"] < 2.0e-4
        and differences["curvature_mass"] < 2.0e-4
        and differences["potential_mass"] < 2.0e-4,
        f"total={differences['total_mass']:.3e} curvature={differences['curvature_mass']:.3e} "
        f"potential={differences['potential_mass']:.3e}",
    )
    ledger.check(
        "mass_ODE_residual",
        residual < 2.0e-4,
        f"max normalized spline-DOP853 residual={residual:.3e}",
    )
    ledger.check(
        "horizon_cross_route_agreement",
        differences["horizon_radius"] < 5.0e-3
        and differences["maximum_compactness"] < 5.0e-4,
        f"horizon={differences['horizon_radius']:.3e} maxC={differences['maximum_compactness']:.3e}",
    )
    ledger.check(
        "homothetic_cross_route_agreement",
        homothetic.success
        and differences["homothetic_minimum"] < 5.0e-4
        and differences["homothetic_scale"] < 5.0e-4,
        f"minimum={differences['homothetic_minimum']:.3e} scale={differences['homothetic_scale']:.3e}",
    )
    ledger.check(
        "trapped_verdict_cross_route_margin",
        horizon > 0.0
        and maximum > 1.0
        and homothetic_minimum > 1.0
        and trapped_cross_route_margin > 100.0,
        f"horizon={horizon:.9f} maxC={maximum:.6e} homC={homothetic_minimum:.6e} "
        f"margin={trapped_cross_route_margin:.3e}",
    )

    payload = {
        "campaign": "P245",
        "attempt": "0003",
        "verdict": "INDEPENDENT_TRAPPED_SURFACE_CORROBORATION",
        "method": {
            "source": "analytic NumPy Chebyshev/angular derivatives",
            "radial_density_nodes": SOURCE_NODES,
            "interpolant": "PCHIP",
            "mass_integrator": "DOP853",
            "rtol": 1.0e-11,
            "atol": 1.0e-12,
            "horizon_finder": "Brent",
        },
        "result": {
            "mass": total_mass,
            "curvature_mass": float(terminal_components[0]),
            "potential_mass": float(terminal_components[1]),
            "horizon_radius": horizon,
            "maximum_compactness": maximum,
            "maximum_compactness_radius": maximum_radius,
            "minimum_f": 1.0 - maximum,
            "exterior_horizon_radius": 2.0 * G_TOTAL * total_mass,
            "homothetic_minimizing_radius": homothetic_scale,
            "homothetic_minimum_maximum_compactness": homothetic_minimum,
            "homothetic_critical_newton_constant": G_TOTAL / homothetic_minimum,
            "mass_ode_residual": residual,
        },
        "relative_differences_from_attempt_0002": differences,
        "trapped_cross_route_signal_to_difference": trapped_cross_route_margin,
        "scope": (
            "Independent corroboration of the frozen-source and homothetic "
            "necessary trapped-surface result; arbitrary M5 profile relaxation "
            "remains outside this evidence."
        ),
    }
    output = HERE / "independent-verdict.json"
    output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    print(f"[DONE] {output.name} written")
    ledger.finish()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
