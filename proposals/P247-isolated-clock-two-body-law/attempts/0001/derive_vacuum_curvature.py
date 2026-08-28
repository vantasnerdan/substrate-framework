"""P247 attempt 0001, gate G1: exact second variation of the committed potential
about the rank-1 projector background in the local frame.

Reported 0044 claim (transfer reference, PR #153 comment): the potential second
variation around the rank-1 projector background P = n n^T is exactly
diag(5, 3, 3, 0, 0, 6)/2 in the local frame, with two exactly flat shear
channels (np and n-a) whose stiffness is purely gradient-generated.

The committed potential (proposals/P240-m5-kinetic-axis/attempts/0041/
cpu_energy.py line 124) is

    V(S) = -Tr(S^2)/2 - Tr(S^3) + (Tr S^2)^2 + 1/2.

This script derives the exact symmetric-block Hessian of V at P in the local
frame (n, p, a) = (e1, e2, e3) on the six-dimensional symmetric channel basis

    B = [n n, p p, a a, sym(n p), sym(n a), sym(p a)],

asserts the exact rational diagonal diag(5, 3, 3, 0, 0, 6)/2, identifies the
flat channels, and cross-checks the massive/zero entries against an
independent finite-difference Hessian of the same committed functional.

Runs offline; no simulations at import. PYTHONPATH=src not required (no
canonical-module import); the functional is restated verbatim with provenance.
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent

# Local frame: director n-hat, polar p-hat, azimuthal a-hat (cpu_energy.py:90-92
# at mu such that the frame is orthonormal; the local-frame statement is
# frame-covariant, so e1, e2, e3 suffice).
n = sp.Matrix([1, 0, 0])
p = sp.Matrix([0, 1, 0])
a = sp.Matrix([0, 0, 1])


def outer(u: sp.Matrix, v: sp.Matrix) -> sp.Matrix:
    return u * v.T


def sym(u: sp.Matrix, v: sp.Matrix) -> sp.Matrix:
    return outer(u, v) + outer(v, u)


CHANNELS = [
    ("nn", outer(n, n)),
    ("pp", outer(p, p)),
    ("aa", outer(a, a)),
    ("np", sym(n, p)),
    ("na", sym(n, a)),
    ("pa", sym(p, a)),
]


def potential(matrix: sp.MatrixBase) -> sp.Expr:
    """Verbatim committed functional (cpu_energy.py:124)."""
    trace_two = sp.trace(matrix * matrix)
    trace_three = sp.trace(matrix * matrix * matrix)
    return -sp.Rational(1, 2) * trace_two - trace_three + trace_two**2 + sp.Rational(1, 2)


def exact_curvature() -> dict:
    background = outer(n, n)
    t = sp.Symbol("t", real=True)
    # Exact quadratic form Q2(w) = V(P + t*w) - V(P) at t^2, equal to the
    # census-convention kernel V2(w) = w^T K w on the coefficient vector c.
    # Since Q2 = sum_ij k_ij c_i c_j (k symmetric), k_ij = d^2 Q2/dc_i dc_j / 2.
    coefficients = sp.symbols("c0:6")
    variation = sum((coefficients[i] * CHANNELS[i][1] for i in range(6)), sp.zeros(3, 3))
    expanded = sp.expand(potential(background + t * variation) - potential(background))
    quadratic = sp.expand(expanded.coeff(t, 2))
    kernel = sp.Matrix(
        6, 6,
        lambda i, j: sp.nsimplify(sp.simplify(sp.diff(quadratic, coefficients[i], coefficients[j]) / 2)),
    )
    off_diagonal = {
        f"{CHANNELS[i][0]}-{CHANNELS[j][0]}": kernel[i, j]
        for i in range(6)
        for j in range(i + 1, 6)
        if kernel[i, j] != 0
    }
    return {
        "kernel": kernel,
        "diagonal": [kernel[i, i] for i in range(6)],
        "off_diagonal_nonzero": off_diagonal,
        "value_at_background": sp.simplify(potential(background)),
    }


def expected_kernel() -> sp.Matrix:
    return sp.Rational(1, 2) * sp.diag(5, 3, 3, 0, 0, 6)




def finite_difference_check(tol_massive: float = 1.0e-6, tol_flat: float = 1.0e-7) -> list[dict]:
    """Independent numeric cross-check of the exact kernel entries.

    A central second difference approximates V''(0) = 2 * k, where k is the
    kernel entry in the V2(w) = w^T K w coefficient convention. The Richardson
    combination kills the h^4 term so the exactly flat channels resolve below
    the 1e-7 absolute gate.
    """
    import numpy as np

    background = np.diag([1.0, 0.0, 0.0])

    def v_numeric(matrix: np.ndarray) -> float:
        trace_two = np.trace(matrix @ matrix)
        trace_three = np.trace(matrix @ matrix @ matrix)
        return -0.5 * trace_two - trace_three + trace_two**2 + 0.5

    def second_derivative(direction: np.ndarray, step: float) -> float:
        def symmetric(distance: float) -> float:
            return v_numeric(background + distance * direction) + v_numeric(background - distance * direction)

        middle = v_numeric(background)
        s1 = symmetric(step) - 2.0 * middle
        s2 = symmetric(2.0 * step) - 2.0 * middle
        return (16.0 * s1 - s2) / (12.0 * step**2)

    rows = []
    for index, (name, channel) in enumerate(CHANNELS):
        direction = np.asarray(channel, dtype=np.float64)
        step = 1.0e-3
        curvature = second_derivative(direction, step)
        expected = 2.0 * float(expected_kernel()[index, index])
        tolerance = tol_flat if expected == 0 else tol_massive * abs(expected)
        rows.append(
            {
                "channel": name,
                "exact_kernel_entry": float(expected_kernel()[index, index]),
                "expected_second_derivative": expected,
                "finite_difference": curvature,
                "abs_error": abs(curvature - expected),
                "tolerance": tolerance,
                "pass": bool(abs(curvature - expected) <= tolerance),
            }
        )
    return rows


def main() -> None:
    result = exact_curvature()
    expected = expected_kernel()
    kernel_agrees = sp.simplify(result["kernel"] - expected) == sp.zeros(6, 6)
    flat_channels = [CHANNELS[i][0] for i in range(6) if result["kernel"][i, i] == 0]
    numeric = finite_difference_check()
    all_numeric_pass = all(row["pass"] for row in numeric)
    payload = {
        "gate": "G1_symbolic",
        "potential_source": "proposals/P240-m5-kinetic-axis/attempts/0041/cpu_energy.py:124 (verbatim)",
        "background": "rank-1 projector n n^T, local frame (n, p, a)",
        "channel_order": [name for name, _ in CHANNELS],
        "exact_kernel": [[str(entry) for entry in row] for row in result["kernel"].tolist()],
        "expected_kernel": "diag(5,3,3,0,0,6)/2",
        "kernel_matches_reported": bool(kernel_agrees),
        "off_diagonal_nonzero": {key: str(value) for key, value in result["off_diagonal_nonzero"].items()},
        "flat_channels": flat_channels,
        "flat_channels_are_np_and_na": bool(flat_channels == ["np", "na"]),
        "potential_at_background_exact_zero": str(result["value_at_background"]),
        "massive_channel_speeds": {
            CHANNELS[i][0]: str(sp.sqrt(result["kernel"][i, i])) for i in range(6) if result["kernel"][i, i] > 0
        },
        "finite_difference_check": numeric,
        "finite_difference_all_pass": all_numeric_pass,
    }
    passed = (
        payload["kernel_matches_reported"]
        and payload["flat_channels_are_np_and_na"]
        and result["value_at_background"] == 0
        and all_numeric_pass
    )
    payload["verdict"] = "ESTABLISHED" if passed else "FAILED"
    (HERE / "g1-curvature.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    print(f"G1 VERDICT: {payload['verdict']}")
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
