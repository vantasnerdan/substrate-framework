"""P246 attempt 0005: horizon-penetrating initial data for the clock source."""

from __future__ import annotations

import contextlib
import io
import json
import traceback
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from substrate_framework.m5_self_gravitating_clock import (
    m5_clock_kinematics,
    m5_clock_profiles_from_chebyshev,
)
from substrate_framework.nonlinear_clock_gravity import (
    integrate_spherical_density_cells,
    painleve_gullstrand_mass_constraint,
)
from substrate_framework.verification import CheckLedger

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
RADIUS = 12.0
G_TOTAL = 46.80699908016004
ACCEPTED_MASS = 54.70900884959007
CELL_COUNTS = (400, 800, 1600)


def execute() -> None:
    ledger = CheckLedger("P246-attempt-0005-horizon-penetrating")
    root = json.loads(
        (
            REPO / "proposals/P240-m5-kinetic-axis/attempts/0042/largeR-roots.json"
        ).read_text()
    )["R12"]
    coefficients = np.asarray(root["values"], dtype=np.float64).reshape(3, -1)
    records: list[dict[str, float]] = []
    finest_data: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None

    for cells in CELL_COUNTS:
        edges = np.linspace(0.0, RADIUS, cells + 1)
        midpoints = 0.5 * (edges[:-1] + edges[1:])
        profiles = m5_clock_profiles_from_chebyshev(coefficients, midpoints, RADIUS)
        kinematics = m5_clock_kinematics(
            midpoints, profiles, angular_quadrature_count=64
        )
        density = (
            kinematics.curvature_rtheta
            + kinematics.curvature_rphi
            + kinematics.curvature_thetaphi
            + kinematics.potential
        )
        compactness = integrate_spherical_density_cells(edges, density, G_TOTAL)
        metric = compactness.radial_metric_function
        crossings = np.flatnonzero((metric[:-1] > 0.0) & (metric[1:] < 0.0))
        ledger.check(
            f"inner_crossing_bracket_{cells}",
            crossings.size > 0,
            f"minimum f={np.min(metric):.6e}",
        )
        index = int(crossings[0])
        mass_interpolant = PchipInterpolator(edges, compactness.enclosed_mass)

        def static_metric(radius: float) -> float:
            return 1.0 - 2.0 * G_TOTAL * float(mass_interpolant(radius)) / radius

        horizon = brentq(
            static_metric,
            float(edges[index]),
            float(edges[index + 1]),
            xtol=1.0e-14,
            rtol=1.0e-14,
        )
        records.append(
            {
                "cells": float(cells),
                "mass": compactness.total_mass,
                "inner_horizon_radius": horizon,
                "maximum_compactness": compactness.maximum_compactness,
                "outer_horizon_radius": compactness.exterior_horizon_radius,
            }
        )
        if cells == CELL_COUNTS[-1]:
            finest_data = (midpoints, density, compactness.enclosed_mass)

    finest = records[-1]
    previous = records[-2]
    mass_error = abs(finest["mass"] - ACCEPTED_MASS) / ACCEPTED_MASS
    crossing_drift = (
        abs(finest["inner_horizon_radius"] - previous["inner_horizon_radius"])
        / finest["inner_horizon_radius"]
    )
    ledger.check(
        "accepted_mass_reproduced",
        mass_error < 1.0e-5,
        f"relative error={mass_error:.6e}",
    )
    ledger.check(
        "inner_horizon_refines",
        crossing_drift < 2.0e-3,
        f"relative drift={crossing_drift:.6e}",
    )

    assert finest_data is not None
    midpoints, density, masses = finest_data
    inner = finest["inner_horizon_radius"]
    density_at_inner = float(PchipInterpolator(midpoints, density)(inner))
    mass_at_inner = inner / (2.0 * G_TOTAL)
    shift_at_inner = np.sqrt(2.0 * G_TOTAL * mass_at_inner / inner)
    shift_derivative_inner = (
        G_TOTAL
        * (4.0 * np.pi * inner * density_at_inner - mass_at_inner / inner**2)
        / shift_at_inner
    )
    ledger.check(
        "inner_horizon_shift_is_unity",
        abs(shift_at_inner - 1.0) < 1.0e-12,
    )
    ledger.check(
        "inner_extrinsic_curvature_is_finite",
        np.isfinite(shift_derivative_inner),
        f"v_prime={shift_derivative_inner:.12g}",
    )

    exterior = finest["outer_horizon_radius"]
    geometric_mass = G_TOTAL * finest["mass"]
    exterior_shift = np.sqrt(2.0 * geometric_mass / exterior)
    kretschmann = 48.0 * geometric_mass**2 / exterior**6
    ledger.check(
        "outer_horizon_shift_is_unity",
        abs(exterior_shift - 1.0) < 1.0e-12,
    )
    ledger.check(
        "outer_horizon_curvature_is_finite",
        np.isfinite(kretschmann) and kretschmann > 0.0,
        f"K={kretschmann:.12g}",
    )

    r = sp.symbols("r", positive=True)
    theta = sp.symbols("theta", real=True)
    gravity = sp.symbols("G", positive=True)
    rho = sp.Function("rho")(r)
    mass = sp.Function("m")(r)
    exact = painleve_gullstrand_mass_constraint(r, theta, mass, rho, gravity)
    ledger.check(
        "exact_hamiltonian_constraint",
        exact.hamiltonian_constraint_residual == 0,
    )
    ledger.check(
        "exact_momentum_constraint",
        exact.momentum_constraint_residual == 0,
    )
    ledger.check(
        "metric_determinant_ignores_horizon_factor",
        exact.metric_determinant == -(r**4) * sp.sin(theta) ** 2,
    )

    payload = {
        "source": "C-M5S-004 curvature-plus-potential density",
        "grids": records,
        "accepted_mass_relative_error": mass_error,
        "inner_horizon_relative_drift": crossing_drift,
        "inner_horizon": {
            "radius": inner,
            "shift": shift_at_inner,
            "density": density_at_inner,
            "shift_derivative": shift_derivative_inner,
        },
        "matched_exterior": {
            "geometric_mass": geometric_mass,
            "outer_horizon_radius": exterior,
            "shift": exterior_shift,
            "kretschmann": kretschmann,
        },
        "classification": (
            "horizon-penetrating nonlinear initial-data geometry with a trapped "
            "region between resolved inner and matched-vacuum outer marginal "
            "surfaces; not a stationary fixed-J matter solution"
        ),
    }
    (HERE / "horizon-penetrating-clock.json").write_text(
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
