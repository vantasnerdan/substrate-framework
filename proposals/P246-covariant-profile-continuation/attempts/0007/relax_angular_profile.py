"""P246 attempt 0007: regular even angular relaxation of the clock split."""

from __future__ import annotations

import contextlib
import io
import json
import traceback
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize

from substrate_framework.m5_self_gravitating_clock import (
    m5_clock_axisymmetric_kinematics,
    m5_clock_axisymmetric_stress,
    m5_clock_profiles_from_chebyshev,
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
    ledger = CheckLedger("P246-attempt-0007-angular-relaxation")
    root = json.loads(
        (
            REPO / "proposals/P240-m5-kinetic-axis/attempts/0042/largeR-roots.json"
        ).read_text()
    )["R12"]
    radial_coefficients = np.asarray(root["values"], dtype=np.float64).reshape(3, -1)
    radial_nodes, raw_weights = leggauss(96)
    radius = 0.5 * RADIUS * (radial_nodes + 1.0)
    radial_weights = 0.5 * RADIUS * raw_weights
    profiles = m5_clock_profiles_from_chebyshev(radial_coefficients, radius, RADIUS)

    def components(
        angular_coefficients: np.ndarray, angular_count: int = 72
    ) -> dict[str, float]:
        kinematics = m5_clock_axisymmetric_kinematics(
            radius,
            profiles,
            angular_quadrature_count=angular_count,
            split_angular_coefficients=angular_coefficients,
        )
        angular_weights = kinematics.angular_weights

        def integrate(values: np.ndarray) -> float:
            return float(
                2.0
                * np.pi
                * np.einsum(
                    "r,a,ra->",
                    radial_weights * radius**2,
                    angular_weights,
                    values,
                )
            )

        curvature = integrate(
            kinematics.curvature_rtheta
            + kinematics.curvature_rphi
            + kinematics.curvature_thetaphi
        )
        potential = integrate(kinematics.potential)
        inertia = integrate(kinematics.flat_inertia_density)
        return {
            "curvature": curvature,
            "potential": potential,
            "inertia": inertia,
            "frequency": 1.0 / (2.0 * inertia),
            "total": curvature + potential + 1.0 / (4.0 * inertia),
        }

    base = components(np.zeros(0))
    records: list[dict[str, object]] = []
    seed = np.zeros(0, dtype=np.float64)
    final_result = None
    for mode_count in (2, 3, 4):
        initial = np.pad(seed, (0, mode_count - seed.size))

        def objective(values: np.ndarray) -> float:
            return components(values)["total"]

        starts = [initial, np.zeros(mode_count)]
        starts.append(np.linspace(0.15, -0.1, mode_count))
        best = None
        for start in starts:
            result = minimize(
                objective,
                start,
                method="L-BFGS-B",
                bounds=[(-3.0, 3.0)] * mode_count,
                options={"ftol": 1.0e-14, "gtol": 5.0e-9, "maxiter": 500},
            )
            if best is None or result.fun < best.fun:
                best = result
        assert best is not None
        seed = np.asarray(best.x, dtype=np.float64)
        final_result = best
        record = components(seed)
        record.update(
            {
                "mode_count": mode_count,
                "coefficients": seed.tolist(),
                "optimizer_success": bool(best.success),
                "optimizer_message": str(best.message),
                "gradient_inf": float(np.max(np.abs(best.jac))),
                "iterations": int(best.nit),
            }
        )
        records.append(record)

    assert final_result is not None
    final_coefficients = np.asarray(final_result.x, dtype=np.float64)
    final = components(final_coefficients)
    refined = components(final_coefficients, angular_count=96)
    energy_reduction = (base["total"] - final["total"]) / base["total"]
    refinement_drift = abs(refined["total"] - final["total"]) / refined["total"]
    ledger.check("optimizer_completed", bool(final_result.success))
    ledger.check(
        "angular_stationarity",
        float(np.max(np.abs(final_result.jac))) < 2.0e-6,
        f"gradient_inf={np.max(np.abs(final_result.jac)):.6e}",
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
    compact_profiles = m5_clock_profiles_from_chebyshev(
        radial_coefficients, midpoints, RADIUS
    )
    compact_kinematics = m5_clock_axisymmetric_kinematics(
        midpoints,
        compact_profiles,
        angular_quadrature_count=64,
        split_angular_coefficients=final_coefficients,
    )
    averaged_static = (
        np.einsum(
            "ra,a->r",
            compact_kinematics.flat_static_density,
            compact_kinematics.angular_weights,
        )
        / 2.0
    )
    compactness = integrate_spherical_density_cells(edges, averaged_static, G_TOTAL)

    final_kinematics = m5_clock_axisymmetric_kinematics(
        radius,
        profiles,
        angular_quadrature_count=72,
        split_angular_coefficients=final_coefficients,
    )
    final_stress = m5_clock_axisymmetric_stress(
        final_kinematics,
        np.ones_like(radius),
        np.ones_like(radius),
        final["frequency"],
    )
    measure = (
        2.0
        * np.pi
        * radial_weights[:, None]
        * radius[:, None] ** 2
        * final_kinematics.angular_weights[None, :]
    )
    energy_scale = float(np.sum(measure * final_stress.energy_density))
    momentum_fraction = float(
        np.sum(measure * np.abs(final_stress.tensor[..., 0, 3])) / energy_scale
    )
    pressure_anisotropy = float(
        np.sum(
            measure
            * np.abs(final_stress.tensor[..., 2, 2] - final_stress.tensor[..., 3, 3])
        )
        / energy_scale
    )

    payload = {
        "base": base,
        "nested_relaxation": records,
        "final_refined": refined,
        "relative_energy_reduction": energy_reduction,
        "refined_energy_relative_drift": refinement_drift,
        "maximum_compactness": compactness.maximum_compactness,
        "critical_G": compactness.critical_newton_constant,
        "mass": compactness.total_mass,
        "momentum_fraction_phi": momentum_fraction,
        "tangential_pressure_anisotropy_fraction": pressure_anisotropy,
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
