"""P246 attempt 0001: reconstruct the clock and test spherical stress closure."""

from __future__ import annotations

import contextlib
import io
import json
import os
from pathlib import Path

for _name in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

import numpy as np
from numpy.polynomial.legendre import leggauss

from substrate_framework.m5_self_gravitating_clock import (
    m5_clock_kinematics,
    m5_clock_profiles_from_chebyshev,
    m5_clock_stress,
)
from substrate_framework.verification import CheckLedger

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
RADIUS = 12.0
EXPECTED_ENERGY = 55.10418259137505
EXPECTED_INERTIA = 0.6325964726638463
EXPECTED_FREQUENCY = 0.7903932785058282


def relative_error(actual: float, expected: float) -> float:
    return abs(actual - expected) / abs(expected)


def integrate(values: np.ndarray, radius: np.ndarray, weights: np.ndarray) -> float:
    return float(np.dot(weights, 4.0 * np.pi * radius**2 * values))


def execute() -> dict[str, object]:
    ledger = CheckLedger("P246-attempt-0001-stress-symmetry")
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
    kinematics = m5_clock_kinematics(radius, profiles, angular_quadrature_count=64)

    curvature = integrate(
        kinematics.curvature_rtheta
        + kinematics.curvature_rphi
        + kinematics.curvature_thetaphi,
        radius,
        radial_weights,
    )
    potential = integrate(kinematics.potential, radius, radial_weights)
    inertia = integrate(kinematics.flat_inertia_density, radius, radial_weights)
    frequency = 1.0 / (2.0 * inertia)
    total_energy = curvature + potential + 1.0 / (4.0 * inertia)
    stress = m5_clock_stress(
        kinematics, np.ones_like(radius), np.ones_like(radius), frequency
    )
    pressure_anisotropy = integrate(
        np.abs(stress.polar_pressure - stress.azimuthal_pressure),
        radius,
        radial_weights,
    )
    mass_scale = integrate(stress.energy_density, radius, radial_weights)
    anisotropy_fraction = pressure_anisotropy / mass_scale

    zero_split = coefficients.copy()
    zero_split[2] = 0.0
    spherical_profiles = m5_clock_profiles_from_chebyshev(zero_split, radius, RADIUS)
    spherical_kinematics = m5_clock_kinematics(
        radius, spherical_profiles, angular_quadrature_count=64
    )
    spherical_inertia = integrate(
        spherical_kinematics.flat_inertia_density, radius, radial_weights
    )
    spherical_frequency = 1.0 / (2.0 * spherical_inertia)
    spherical_stress = m5_clock_stress(
        spherical_kinematics,
        np.ones_like(radius),
        np.ones_like(radius),
        spherical_frequency,
    )
    spherical_anisotropy = integrate(
        np.abs(spherical_stress.polar_pressure - spherical_stress.azimuthal_pressure),
        radius,
        radial_weights,
    ) / integrate(spherical_stress.energy_density, radius, radial_weights)

    ledger.check(
        "inertia_reproduced",
        relative_error(inertia, EXPECTED_INERTIA) < 2.0e-7,
        f"inertia={inertia:.12g}",
    )
    ledger.check(
        "frequency_reproduced",
        relative_error(frequency, EXPECTED_FREQUENCY) < 2.0e-7,
        f"frequency={frequency:.12g}",
    )
    ledger.check(
        "total_energy_reproduced",
        relative_error(total_energy, EXPECTED_ENERGY) < 2.0e-7,
        f"energy={total_energy:.12g}",
    )
    ledger.check(
        "accepted_split_stress_is_anisotropic",
        anisotropy_fraction > 1.0e-5,
        f"integrated pressure anisotropy fraction={anisotropy_fraction:.6e}",
    )
    ledger.check(
        "split_deletion_restores_spherical_stress",
        spherical_anisotropy < 1.0e-11,
        f"zero-split anisotropy fraction={spherical_anisotropy:.6e}",
    )

    payload: dict[str, object] = {
        "background": "P240 R12 order-16 family-S root",
        "quadrature": {"radial": 128, "angular": 64},
        "curvature_energy": curvature,
        "potential_energy": potential,
        "inertia": inertia,
        "frequency": frequency,
        "total_energy": total_energy,
        "pressure_anisotropy_fraction": anisotropy_fraction,
        "zero_split_pressure_anisotropy_fraction": spherical_anisotropy,
        "classification": (
            "the accepted split-active clock has axisymmetric, not spherical, "
            "stress; a co-solved continuation must use axisymmetric gravity or "
            "construct a distinct split-free stationary branch"
        ),
    }
    (HERE / "stress-symmetry.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    ledger.finish()
    return payload


def main() -> int:
    capture = io.StringIO()
    with contextlib.redirect_stdout(capture):
        execute()
    output = capture.getvalue()
    (HERE / "stdout.txt").write_text(output)
    print(output, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
