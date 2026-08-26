"""P246 attempt 0002: refined reconstruction and stress-symmetry verdict."""

from __future__ import annotations

import contextlib
import io
import json
import os
import traceback
from pathlib import Path

for _name in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

import numpy as np  # noqa: E402
from numpy.polynomial.legendre import leggauss  # noqa: E402

from substrate_framework.m5_self_gravitating_clock import (  # noqa: E402
    m5_clock_kinematics,
    m5_clock_profiles_from_chebyshev,
    m5_clock_stress,
)
from substrate_framework.verification import CheckLedger  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
RADIUS = 12.0
EXPECTED_ENERGY = 55.10418259137505
EXPECTED_INERTIA = 0.6325964726638463
EXPECTED_FREQUENCY = 0.7903932785058282
LADDER = ((64, 32), (128, 64), (256, 96))


def relative_error(actual: float, expected: float) -> float:
    return abs(actual - expected) / abs(expected)


def evaluate(
    coefficients: np.ndarray, radial_count: int, angular_count: int
) -> dict[str, float]:
    nodes, raw_weights = leggauss(radial_count)
    radius = 0.5 * RADIUS * (nodes + 1.0)
    weights = 0.5 * RADIUS * raw_weights
    profiles = m5_clock_profiles_from_chebyshev(coefficients, radius, RADIUS)
    kinematics = m5_clock_kinematics(
        radius, profiles, angular_quadrature_count=angular_count
    )

    def integrate(values: np.ndarray) -> float:
        return float(np.dot(weights, 4.0 * np.pi * radius**2 * values))

    curvature = integrate(
        kinematics.curvature_rtheta
        + kinematics.curvature_rphi
        + kinematics.curvature_thetaphi
    )
    potential = integrate(kinematics.potential)
    inertia = integrate(kinematics.flat_inertia_density)
    frequency = 1.0 / (2.0 * inertia)
    stress = m5_clock_stress(
        kinematics, np.ones_like(radius), np.ones_like(radius), frequency
    )
    mass_scale = integrate(stress.energy_density)
    anisotropy = (
        integrate(np.abs(stress.polar_pressure - stress.azimuthal_pressure))
        / mass_scale
    )
    return {
        "radial_count": float(radial_count),
        "angular_count": float(angular_count),
        "curvature_energy": curvature,
        "potential_energy": potential,
        "inertia": inertia,
        "frequency": frequency,
        "total_energy": curvature + potential + 1.0 / (4.0 * inertia),
        "pressure_anisotropy_fraction": anisotropy,
    }


def execute() -> None:
    ledger = CheckLedger("P246-attempt-0002-stress-symmetry")
    root = json.loads(
        (
            REPO / "proposals/P240-m5-kinetic-axis/attempts/0042/largeR-roots.json"
        ).read_text()
    )["R12"]
    coefficients = np.asarray(root["values"], dtype=np.float64).reshape(3, -1)
    records = [evaluate(coefficients, *counts) for counts in LADDER]
    finest = records[-1]
    previous = records[-2]

    zero_split = coefficients.copy()
    zero_split[2] = 0.0
    mutation = evaluate(zero_split, *LADDER[-1])

    finest_drift = relative_error(finest["inertia"], previous["inertia"])
    ledger.check(
        "new_evaluator_refines",
        finest_drift < 2.0e-8,
        f"128x64 to 256x96 inertia drift={finest_drift:.6e}",
    )
    ledger.check(
        "native_inertia_reproduced",
        relative_error(finest["inertia"], EXPECTED_INERTIA) < 3.0e-6,
        f"inertia={finest['inertia']:.12g}",
    )
    ledger.check(
        "native_frequency_reproduced",
        relative_error(finest["frequency"], EXPECTED_FREQUENCY) < 3.0e-6,
        f"frequency={finest['frequency']:.12g}",
    )
    ledger.check(
        "native_total_energy_reproduced",
        relative_error(finest["total_energy"], EXPECTED_ENERGY) < 3.0e-6,
        f"energy={finest['total_energy']:.12g}",
    )
    ledger.check(
        "accepted_split_stress_is_anisotropic",
        finest["pressure_anisotropy_fraction"] > 1.0e-5,
        (f"pressure anisotropy fraction={finest['pressure_anisotropy_fraction']:.6e}"),
    )
    ledger.check(
        "split_deletion_restores_spherical_stress",
        mutation["pressure_anisotropy_fraction"] < 1.0e-11,
        (
            "zero-split pressure anisotropy fraction="
            f"{mutation['pressure_anisotropy_fraction']:.6e}"
        ),
    )

    payload = {
        "background": "P240 R12 order-16 family-S root",
        "ladder": records,
        "finest_inertia_relative_drift": finest_drift,
        "native_comparators": {
            "energy": EXPECTED_ENERGY,
            "inertia": EXPECTED_INERTIA,
            "frequency": EXPECTED_FREQUENCY,
        },
        "zero_split_mutation": mutation,
        "classification": (
            "the accepted split-active clock has axisymmetric, not spherical, "
            "stress; a co-solved continuation must use axisymmetric gravity or "
            "construct a distinct split-free stationary branch"
        ),
    }
    (HERE / "stress-symmetry-refined.json").write_text(
        json.dumps(payload, indent=2) + "\n"
    )
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
