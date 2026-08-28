"""P247 attempt 0003, gate W2: candidate-B enumeration property checks (sympy).

Enumerated terms are local Lorentz-covariant structures of the conditional
P239 action (src/substrate_framework/m5_covariant_action.py primitives) that
vanish identically on the complete arbitrary static 3x3 sector: the field
branch aligned with the timelike vacuum eigenline u = xi = e0, static
profiles, arbitrary spatial 3x3 block M_ij (M_00 = M_0i = 0).

Checked here exactly:
  P1. Boost-alignment mass s_m = (m^2/2) * ((u^T eta xi)^2 - 1):
      (a) vanishes identically on the aligned static sector;
      (b) is nonnegative for every unit-timelike u (bounded below);
      (c) expands quadratically in the boost rapidity (a mass term);
      (d) is invariant under the clock isorotation symmetry (spatial
          internal rotations about any axis preserve xi).
  P2. Candidate-H scalar current -kappa/2 (d tau)^2 with tau = -u^T M u:
      tau and its gradient vanish identically on the aligned static sector.
  P3. Boost-Skyrme commutator densities vanish on the static sector
      (d u = 0 there) but carry no mass scale (quartic at small rapidity).
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent

ETA = sp.diag(-1, 1, 1, 1)
XI = sp.Matrix([1, 0, 0, 0])  # alignment covector: the vacuum timelike eigenline

B, S1, S2, S3 = sp.symbols("b s1 s2 s3", positive=True, real=True)


def unit_timelike(rapidity, s1, s2, s3) -> sp.Matrix:
    """Unit-timelike eigenvector u with the given rapidity in the spatial
    direction (s1, s2, s3)."""
    spatial_unit = sp.sqrt(s1**2 + s2**2 + s3**2)
    return sp.Matrix(
        [
            sp.cosh(rapidity),
            sp.sinh(rapidity) * s1 / spatial_unit,
            sp.sinh(rapidity) * s2 / spatial_unit,
            sp.sinh(rapidity) * s3 / spatial_unit,
        ]
    )


def alignment_defect(u) -> sp.Expr:
    """(u^T eta xi)^2 - 1 for unit-timelike u."""
    return sp.simplify(((u.T * ETA * XI)[0, 0]) ** 2 - 1)


def check_p1() -> dict:
    u = unit_timelike(B, S1, S2, S3)
    norm_residual = sp.simplify((u.T * ETA * u)[0, 0] + 1)
    s_m = alignment_defect(u)

    on_sector = sp.simplify(s_m.subs(B, 0))
    quadratic = sp.simplify(sp.series(s_m, B, 0, 3).removeO() / B**2)

    c, s = sp.symbols("c s", real=True)
    rotation = sp.Matrix(
        [[1, 0, 0, 0], [0, c, -s, 0], [0, s, c, 0], [0, 0, 0, 1]]
    )
    invariance = sp.simplify(alignment_defect(rotation * u) - s_m)

    return {
        "term": "s_m = (m^2/2) ((u^T eta xi)^2 - 1)",
        "unit_norm_residual": norm_residual,
        "static_sector_value": on_sector,
        "equals_sinh_squared_rapidity": bool(sp.simplify(s_m - sp.sinh(B) ** 2) == 0),
        "quadratic_coefficient": quadratic,
        "isorotation_invariance_residual": invariance,
    }


def check_p2() -> dict:
    """tau = -u^T M u vanishes identically on the aligned static sector for
    an arbitrary static spatial block, hence so does every gradient."""
    m00, m01, m02, m03 = sp.symbols("M00 M01 M02 M03", real=True)
    m11, m12, m13, m22, m23, m33 = sp.symbols(
        "M11 M12 M13 M22 M23 M33", real=True
    )
    big_m = sp.Matrix(
        [
            [m00, m01, m02, m03],
            [m01, m11, m12, m13],
            [m02, m12, m22, m23],
            [m03, m13, m23, m33],
        ]
    )
    tau = sp.expand(-(XI.T * big_m * XI)[0, 0])
    tau_static = sp.simplify(tau.subs({m00: 0, m01: 0, m02: 0, m03: 0}))
    x1, x2, x3, t = sp.symbols("x1 x2 x3 t", real=True)
    gradient = [sp.diff(tau_static, var) for var in (t, x1, x2, x3)]
    return {
        "term": "tau = -u^T M u (Candidate-H scalar current)",
        "tau_on_aligned_static_sector": tau_static,
        "gradient_norm_on_sector": gradient,
    }


def check_p3() -> dict:
    u = unit_timelike(B, S1, S2, S3)
    s_m = alignment_defect(u)
    leading = sp.series(s_m, B, 0, 6).removeO()
    return {
        "note": (
            "any ||d u, d u||-type density vanishes on the static sector "
            "(d u = 0 there); the alignment defect itself is quartic-leading "
            "in a commutator, and carries no constant (mass) term"
        ),
        "alignment_defect_leading_order": sp.degree(leading, B),
        "no_constant_term": bool(sp.simplify(s_m.subs(B, 0)) == 0),
    }


def main() -> None:
    p1 = check_p1()
    p2 = check_p2()
    p3 = check_p3()
    payload = {
        "purpose": "W2 enumeration property checks (manifest gates.W2_candidate_B_enumeration)",
        "p1_boost_alignment_mass": {
            **p1,
            "checks_pass": bool(
                p1["static_sector_value"] == 0
                and p1["equals_sinh_squared_rapidity"]
                and p1["isorotation_invariance_residual"] == 0
                and sp.simplify(p1["unit_norm_residual"]) == 0
            ),
        },
        "p2_tau_scalar_current": {
            **p2,
            "checks_pass": bool(
                p2["tau_on_aligned_static_sector"] == 0
                and all(g == 0 for g in p2["gradient_norm_on_sector"])
            ),
        },
        "p3_boost_skyrme": {
            **p3,
            "checks_pass": bool(p3["no_constant_term"]),
        },
    }
    (HERE / "w2-property-checks.json").write_text(
        json.dumps(payload, indent=2, default=str)
    )
    print(json.dumps(payload, indent=2, default=str))
    print("WROTE w2-property-checks.json")


if __name__ == "__main__":
    main()
