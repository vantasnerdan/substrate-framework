"""P246 attempt 0008: corrected regular angular-profile relaxation."""

from __future__ import annotations

import contextlib
import io
import json
import traceback
from pathlib import Path

import numpy as np

from substrate_framework.m5_self_gravitating_clock import (
    m5_clock_axisymmetric_kinematics,
    m5_clock_profiles_from_chebyshev,
    relax_m5_clock_split_angular_profile,
)
from substrate_framework.nonlinear_clock_gravity import (
    integrate_spherical_density_cells,
)
from substrate_framework.verification import CheckLedger

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
RADIUS = 12.0
G_TOTAL = 46.80699908016004


def execute() -> None:
    ledger = CheckLedger("P246-attempt-0008-angular-relaxation")
    root = json.loads(
        (
            REPO / "proposals/P240-m5-kinetic-axis/attempts/0042/largeR-roots.json"
        ).read_text()
    )["R12"]
    radial_coefficients = np.asarray(root["values"], dtype=np.float64).reshape(3, -1)
    relaxation = relax_m5_clock_split_angular_profile(
        radial_coefficients,
        domain_radius=RADIUS,
        mode_count=4,
        radial_quadrature_count=96,
        angular_quadrature_count=72,
        refinement_angular_count=96,
    )
    energy_reduction = (
        relaxation.base_total_energy - relaxation.total_energy
    ) / relaxation.base_total_energy
    refinement_drift = (
        abs(relaxation.refined_total_energy - relaxation.total_energy)
        / relaxation.refined_total_energy
    )
    ledger.check("optimizer_completed", relaxation.optimizer_success)
    ledger.check(
        "centered_gradient_is_scale_small",
        relaxation.gradient_scale_relative < 1.0e-6,
        f"relative gradient={relaxation.gradient_scale_relative:.6e}",
    )
    ledger.check(
        "angular_energy_refines",
        refinement_drift < 2.0e-7,
        f"relative drift={refinement_drift:.6e}",
    )
    ledger.check(
        "angular_relaxation_is_nontrivial",
        energy_reduction > 1.0e-6,
        f"relative reduction={energy_reduction:.6e}",
    )

    cells = 1600
    edges = np.linspace(0.0, RADIUS, cells + 1)
    midpoints = 0.5 * (edges[:-1] + edges[1:])
    profiles = m5_clock_profiles_from_chebyshev(radial_coefficients, midpoints, RADIUS)
    kinematics = m5_clock_axisymmetric_kinematics(
        midpoints,
        profiles,
        angular_quadrature_count=64,
        split_angular_coefficients=relaxation.coefficients,
    )
    static_density = (
        np.einsum(
            "ra,a->r",
            kinematics.flat_static_density,
            kinematics.angular_weights,
        )
        / 2.0
    )
    compactness = integrate_spherical_density_cells(edges, static_density, G_TOTAL)
    payload = {
        "coefficients": relaxation.coefficients.tolist(),
        "base_total_energy": relaxation.base_total_energy,
        "relaxed_total_energy": relaxation.total_energy,
        "refined_total_energy": relaxation.refined_total_energy,
        "curvature_energy": relaxation.curvature_energy,
        "potential_energy": relaxation.potential_energy,
        "inertia": relaxation.inertia,
        "frequency": relaxation.frequency,
        "relative_energy_reduction": energy_reduction,
        "centered_gradient_inf": relaxation.gradient_inf_norm,
        "centered_gradient_scale_relative": (relaxation.gradient_scale_relative),
        "refinement_relative_drift": refinement_drift,
        "static_mass": compactness.total_mass,
        "maximum_compactness": compactness.maximum_compactness,
        "critical_G": compactness.critical_newton_constant,
        "classification": (
            "horizonless" if compactness.maximum_compactness < 1.0 else "trapped"
        ),
    }
    (HERE / "angular-relaxation.json").write_text(json.dumps(payload, indent=2) + "\n")
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
