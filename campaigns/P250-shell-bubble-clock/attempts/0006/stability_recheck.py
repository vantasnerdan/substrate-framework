#!/usr/bin/env python
"""Recheck C-M5W-008 against the repaired attempt-0006 bag data."""

from __future__ import annotations

import json
import math
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent


def main() -> int:
    checks: list[bool] = []

    def check(name: str, condition: bool) -> None:
        checks.append(bool(condition))
        print(f"[{'PASS' if condition else 'FAIL'}] {name}")

    sigma, pressure, radius = sp.symbols("sigma p R", positive=True)
    routhian = (
        4 * sp.pi * sigma * radius**2
        - sp.Rational(4, 3) * sp.pi * pressure * radius**3
    )
    stationary_radius = 2 * sigma / pressure
    check(
        "thin-wall stationarity is R=2 sigma/p",
        sp.simplify(sp.diff(routhian, radius).subs(
            radius, stationary_radius
        )) == 0,
    )
    curvature = sp.simplify(
        sp.diff(routhian, radius, 2).subs(radius, stationary_radius)
    )
    check(
        "fixed-frequency radial curvature is exactly -8 pi sigma",
        sp.simplify(curvature + 8 * sp.pi * sigma) == 0
        and sp.Lt(curvature, 0) is sp.true,
    )

    data = json.loads((HERE / "bag_results_repaired.json").read_text())
    rows = data["rungs"]
    check(
        "all seven repaired rungs are converged status-zero solutions",
        len(rows) == 7 and all(row["solver_status"] == 0 for row in rows),
    )
    omegas = [row["omega"] for row in rows]
    charges = [row["Q"] for row in rows]
    radii = [row["R_m"] for row in rows]
    charge_secants = [
        (charges[index + 1] - charges[index])
        / (omegas[index + 1] - omegas[index])
        for index in range(len(rows) - 1)
    ]
    check(
        "all six charge secants are strictly negative",
        all(value < 0 for value in charge_secants),
    )
    check(
        "all six radius secants are strictly negative",
        all(radii[index + 1] < radii[index]
            for index in range(len(rows) - 1)),
    )
    check(
        "finite-difference physical-energy envelope stays below 1.9e-4",
        max(row["relative_error"] for row in data["envelope"]) <= 1.9e-4,
    )
    check(
        "physical charge is reconstructed as omega times inertia",
        all(
            math.isclose(row["Q"], row["omega"] * row["I"], rel_tol=2e-15)
            for row in rows
        ),
    )
    check(
        "physical energy is reconstructed as F_omega plus omega Q",
        all(
            math.isclose(
                row["E"], row["E_omega"] + row["omega"] * row["Q"],
                rel_tol=2e-15,
            )
            for row in rows
        ),
    )

    passed = sum(checks)
    print(
        f"ALL {passed} CHECKS PASS"
        if passed == len(checks)
        else f"{passed}/{len(checks)} CHECKS PASS -- FAILURES PRESENT"
    )
    return 0 if passed == len(checks) else 1


if __name__ == "__main__":
    raise SystemExit(main())
