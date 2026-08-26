"""P246 attempt 0009: time-reversal-even counterrotating stress ledger."""

from __future__ import annotations

import contextlib
import io
import json
import traceback
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from substrate_framework.m5_self_gravitating_clock import (
    m5_clock_axisymmetric_kinematics,
    m5_clock_axisymmetric_stress,
    m5_clock_profiles_from_chebyshev,
)
from substrate_framework.verification import CheckLedger

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
RADIUS = 12.0


def execute() -> None:
    ledger = CheckLedger("P246-attempt-0009-counterrotation")
    root = json.loads(
        (
            REPO / "proposals/P240-m5-kinetic-axis/attempts/0042/largeR-roots.json"
        ).read_text()
    )["R12"]
    coefficients = np.asarray(root["values"], dtype=np.float64).reshape(3, -1)
    nodes, raw_weights = leggauss(96)
    radius = 0.5 * RADIUS * (nodes + 1.0)
    radial_weights = 0.5 * RADIUS * raw_weights
    profiles = m5_clock_profiles_from_chebyshev(coefficients, radius, RADIUS)
    kinematics = m5_clock_axisymmetric_kinematics(
        radius, profiles, angular_quadrature_count=64
    )
    measure = (
        2.0
        * np.pi
        * radial_weights[:, None]
        * radius[:, None] ** 2
        * kinematics.angular_weights[None, :]
    )
    inertia = float(np.sum(measure * kinematics.flat_inertia_density))
    omega = 1.0 / (2.0 * inertia)
    positive = m5_clock_axisymmetric_stress(
        kinematics, np.ones_like(radius), np.ones_like(radius), omega
    ).tensor
    negative = m5_clock_axisymmetric_stress(
        kinematics, np.ones_like(radius), np.ones_like(radius), -omega
    ).tensor
    counter = 0.5 * (positive + negative)
    scale = float(np.max(np.abs(positive)))
    momentum_error = float(np.max(np.abs(counter[..., 0, 1:4])) / scale)
    diagonal_indices = np.arange(4)
    diagonal_error = float(
        np.max(
            np.abs(
                counter[..., diagonal_indices, diagonal_indices]
                - positive[..., diagonal_indices, diagonal_indices]
            )
        )
        / scale
    )
    spatial_error = float(
        np.max(np.abs(counter[..., 1:4, 1:4] - positive[..., 1:4, 1:4])) / scale
    )
    positive_energy = float(np.sum(measure * positive[..., 0, 0]))
    counter_energy = float(np.sum(measure * counter[..., 0, 0]))
    energy_error = abs(counter_energy - positive_energy) / positive_energy

    ledger.check("momentum_cancels", momentum_error < 1.0e-12)
    ledger.check("diagonal_stress_is_preserved", diagonal_error < 1.0e-12)
    ledger.check("spatial_shear_is_preserved", spatial_error < 1.0e-12)
    ledger.check("energy_is_preserved", energy_error < 1.0e-12)
    payload = {
        "frequency_magnitude": omega,
        "momentum_cancellation_relative_error": momentum_error,
        "diagonal_preservation_relative_error": diagonal_error,
        "spatial_shear_preservation_relative_error": spatial_error,
        "energy_preservation_relative_error": energy_error,
        "energy": positive_energy,
        "classification": (
            "counterrotation removes frame-dragging momentum exactly but leaves "
            "the energy, pressures, spatial shear and compactness source "
            "unchanged; it is not an escape from the gravity wall"
        ),
    }
    (HERE / "counterrotation.json").write_text(json.dumps(payload, indent=2) + "\n")
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
