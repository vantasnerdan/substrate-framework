#!/usr/bin/env python3
"""Exact algebra checks for P253/0024; no numerical spectral claims."""

import sympy as sp


def require_zero(matrix: sp.Matrix, label: str) -> None:
    if any(sp.simplify(entry) != 0 for entry in matrix):
        raise AssertionError(label)


# A common positive physical factor cancels from i_X Omega = dH.
f = sp.symbols("f", nonzero=True)
omega = sp.Matrix([[0, 1], [-1, 0]])
gradient = sp.Matrix(sp.symbols("h1 h2"))
x_normalized = omega.inv() * gradient
x_physical = (f * omega).inv() * (f * gradient)
require_zero(x_physical - x_normalized, "physical normalization")

# The finite-rank counterexample is symplectic and has exact exponential powers.
x = sp.symbols("x", positive=True)
j_form = sp.Matrix([[0, 1], [-1, 0]])
h_block = sp.diag(x, 1 / x)
require_zero(h_block.T * j_form * h_block - j_form, "hyperbolic symplectic block")
for n in range(1, 7):
    require_zero(h_block**n - sp.diag(x**n, x**(-n)), f"hyperbolic power n={n}")

# Verify the exact full-flag recurrence on scalar representatives of each block.
a, b, c, d = sp.symbols("a b c d")
m = sp.Matrix([[a, 0, b], [c, 1, d], [0, 0, 1]])
for n in range(1, 7):
    b_n = sum((a**k * b for k in range(n)), sp.Integer(0))
    c_n = sum((c * a**k for k in range(n)), sp.Integer(0))
    d_n = n * d + sum(
        ((n - 1 - k) * c * a**k * b for k in range(n - 1)),
        sp.Integer(0),
    )
    expected = sp.Matrix([[a**n, 0, b_n], [c_n, 1, d_n], [0, 0, 1]])
    require_zero(m**n - expected, f"flag power n={n}")

# The registered coboundary coupling telescopes exactly.
y = sp.symbols("y")
for n in range(1, 7):
    telescoping_sum = sum((a**k * (1 - a) * y for k in range(n)), sp.Integer(0))
    if sp.simplify(telescoping_sum - (1 - a**n) * y) != 0:
        raise AssertionError(f"coboundary telescope n={n}")

print("physical_common_factor: PASS")
print("symplectic_hyperbolic_block: PASS")
print("full_flag_power_identity_n1_to_n6: PASS")
print("coboundary_telescope_n1_to_n6: PASS")
print("checks: 4 passed, 0 failed")
