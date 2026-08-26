"""P246 attempt 0003: exact zero-inertia identity and axisymmetric stress."""

from __future__ import annotations

import contextlib
import io
import json
import traceback
from pathlib import Path

import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss

from substrate_framework.m5_self_gravitating_clock import (
    m5_clock_axisymmetric_kinematics,
    m5_clock_profiles_from_chebyshev,
)
from substrate_framework.verification import CheckLedger

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
RADIUS = 12.0


def radial_angular_integral(
    values: np.ndarray,
    radius: np.ndarray,
    radial_weights: np.ndarray,
    angular_weights: np.ndarray,
) -> float:
    radial_integral = np.einsum("r,ra->a", radial_weights * radius**2, values)
    return float(2.0 * np.pi * np.dot(angular_weights, radial_integral))


def execute() -> None:
    ledger = CheckLedger("P246-attempt-0003-axisymmetric-source")

    theta, q_symbol, tangent_symbol = sp.symbols("theta q tangent", real=True)
    director = sp.Matrix([sp.sin(theta), 0, sp.cos(theta)])
    nx, ny, nz = director
    generator = sp.Matrix([[0, -nz, ny], [nz, 0, -nx], [-ny, nx, 0]])
    uniaxial = tangent_symbol * sp.eye(3) + q_symbol * (director * director.T)
    exact_response = sp.simplify(generator * uniaxial + uniaxial * generator.T)
    ledger.check(
        "split_free_clock_response_exactly_zero",
        exact_response == sp.zeros(3),
        f"response={exact_response}",
    )

    root = json.loads(
        (
            REPO / "proposals/P240-m5-kinetic-axis/attempts/0042/largeR-roots.json"
        ).read_text()
    )["R12"]
    coefficients = np.asarray(root["values"], dtype=np.float64).reshape(3, -1)
    nodes, raw_weights = leggauss(128)
    radius = 0.5 * RADIUS * (nodes + 1.0)
    radial_weights = 0.5 * RADIUS * raw_weights
    profiles = m5_clock_profiles_from_chebyshev(coefficients, radius, RADIUS)
    angular = m5_clock_axisymmetric_kinematics(
        radius, profiles, angular_quadrature_count=64
    )
    inertia = radial_angular_integral(
        angular.flat_inertia_density,
        radius,
        radial_weights,
        angular.angular_weights,
    )
    omega = 1.0 / (2.0 * inertia)

    e_r = omega**2 * angular.clock_r
    e_theta = omega**2 * angular.clock_theta
    e_phi = omega**2 * angular.clock_phi
    b_rtheta = angular.curvature_rtheta
    b_rphi = angular.curvature_rphi
    b_thetaphi = angular.curvature_thetaphi
    potential = angular.potential
    rho = e_r + e_theta + e_phi + b_rtheta + b_rphi + b_thetaphi + potential
    p_theta = e_r - e_theta + e_phi + b_rtheta - b_rphi + b_thetaphi - potential
    p_phi = e_r + e_theta - e_phi - b_rtheta + b_rphi + b_thetaphi - potential
    mass_scale = radial_angular_integral(
        rho, radius, radial_weights, angular.angular_weights
    )
    pressure_anisotropy = (
        radial_angular_integral(
            np.abs(p_theta - p_phi),
            radius,
            radial_weights,
            angular.angular_weights,
        )
        / mass_scale
    )
    energy_per_mu = 2.0 * np.pi * np.einsum("r,ra->a", radial_weights * radius**2, rho)
    angular_variation = float(
        (np.max(energy_per_mu) - np.min(energy_per_mu)) / np.mean(energy_per_mu)
    )

    zero_split = coefficients.copy()
    zero_split[2] = 0.0
    zero_profiles = m5_clock_profiles_from_chebyshev(zero_split, radius, RADIUS)
    zero_angular = m5_clock_axisymmetric_kinematics(
        radius, zero_profiles, angular_quadrature_count=64
    )
    zero_inertia_max = float(np.max(np.abs(zero_angular.flat_inertia_density)))

    ledger.check(
        "split_free_numeric_inertia_vanishes",
        zero_inertia_max < 1.0e-24,
        f"max inertia density={zero_inertia_max:.6e}",
    )
    ledger.check(
        "accepted_clock_inertia_is_nonzero",
        inertia > 0.6,
        f"I={inertia:.12g} omega={omega:.12g}",
    )
    ledger.check(
        "accepted_energy_depends_on_polar_angle",
        angular_variation > 1.0e-6,
        f"angular variation fraction={angular_variation:.6e}",
    )
    ledger.check(
        "accepted_tangential_pressures_are_anisotropic",
        pressure_anisotropy > 1.0e-5,
        f"pressure anisotropy fraction={pressure_anisotropy:.6e}",
    )

    payload = {
        "background": "P240 R12 order-16 family-S root",
        "exact_split_free_clock_response": str(exact_response),
        "split_free_inertia_density_max": zero_inertia_max,
        "accepted_inertia": inertia,
        "accepted_frequency": omega,
        "angular_energy_variation_fraction": angular_variation,
        "pressure_anisotropy_fraction": pressure_anisotropy,
        "classification": (
            "nonzero clock charge requires the tangent-splitting channel in "
            "this accepted ansatz; that channel produces resolved axisymmetric "
            "stress, so the co-solved gravitational continuation must be "
            "axisymmetric rather than spherical"
        ),
    }
    (HERE / "axisymmetric-source.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    ledger.finish()


def main() -> int:
    capture = io.StringIO()
    failed = False
    with contextlib.redirect_stdout(capture), contextlib.redirect_stderr(capture):
        try:
            execute()
        except Exception:
            failed = True
            traceback.print_exc()
    output = capture.getvalue()
    (HERE / "stdout.txt").write_text(output)
    print(output, end="")
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
