"""P246 attempt 0004: select the axisymmetric metric class from full stress."""

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
    ledger = CheckLedger("P246-attempt-0004-complete-stress")
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
    stress = m5_clock_axisymmetric_stress(
        kinematics, np.ones_like(radius), np.ones_like(radius), omega
    )
    reversed_stress = m5_clock_axisymmetric_stress(
        kinematics, np.ones_like(radius), np.ones_like(radius), -omega
    )
    tensor = stress.tensor
    reversed_tensor = reversed_stress.tensor
    energy_scale = float(np.sum(measure * tensor[..., 0, 0]))
    momentum_integrals = np.array(
        [np.sum(measure * np.abs(tensor[..., 0, axis])) for axis in range(1, 4)],
        dtype=np.float64,
    )
    momentum_fractions = momentum_integrals / energy_scale
    spatial_offdiagonal = np.array(
        [
            np.sum(measure * np.abs(tensor[..., left, right])) / energy_scale
            for left in range(1, 4)
            for right in range(left + 1, 4)
        ],
        dtype=np.float64,
    )
    symmetry_error = float(
        np.max(np.abs(tensor - np.swapaxes(tensor, -1, -2))) / np.max(np.abs(tensor))
    )
    trace = np.einsum("ab,...ab->...", np.diag([-1.0, 1.0, 1.0, 1.0]), tensor)
    trace_error = float(
        np.max(np.abs(trace + 4.0 * kinematics.potential))
        / max(1.0, np.max(np.abs(tensor)))
    )
    diagonal_indices = np.arange(4)
    diagonal_error = float(
        np.max(
            np.abs(
                tensor[..., diagonal_indices, diagonal_indices]
                - reversed_tensor[..., diagonal_indices, diagonal_indices]
            )
        )
        / np.max(np.abs(tensor))
    )
    momentum_reversal_error = float(
        np.max(np.abs(tensor[..., 0, 1:4] + reversed_tensor[..., 0, 1:4]))
        / max(1.0, np.max(np.abs(tensor[..., 0, 1:4])))
    )

    ledger.check("stress_tensor_is_symmetric", symmetry_error < 1.0e-12)
    ledger.check("quartic_trace_identity", trace_error < 1.0e-10)
    ledger.check("omega_reversal_preserves_diagonal", diagonal_error < 1.0e-12)
    ledger.check("omega_reversal_flips_momentum", momentum_reversal_error < 1.0e-10)
    maximum_momentum_fraction = float(np.max(momentum_fractions))
    metric_class = (
        "static_axisymmetric"
        if maximum_momentum_fraction <= 1.0e-10
        else "stationary_axisymmetric"
    )
    payload = {
        "inertia": inertia,
        "frequency": omega,
        "energy_scale": energy_scale,
        "momentum_fractions_r_theta_phi": momentum_fractions.tolist(),
        "spatial_offdiagonal_fractions_rtheta_rphi_thetaphi": (
            spatial_offdiagonal.tolist()
        ),
        "symmetry_relative_error": symmetry_error,
        "trace_relative_error": trace_error,
        "omega_reversal_diagonal_relative_error": diagonal_error,
        "omega_reversal_momentum_relative_error": momentum_reversal_error,
        "selected_metric_class": metric_class,
    }
    (HERE / "complete-stress.json").write_text(json.dumps(payload, indent=2) + "\n")
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
