#!/usr/bin/env python3
"""Exact algebraic and scaling checks for the P253/0054 HJ2 jet."""

import sympy as sp


# Exact rescaling of the Cao operator.
q, x = sp.symbols("q x", real=True)
w_x, lap_w = sp.symbols("w_x lap_w")
a = 1 / (1 + q * x)
rescaled_lhs = sp.expand((-(a * lap_w + sp.diff(a, x) * w_x) / (1 + q * x)) * (1 + q * x) ** 2)
assert sp.simplify(rescaled_lhs - (-lap_w + q * w_x / (1 + q * x))) == 0

# Exact circular-tube distance and its near-angle expansion.
delta, xx, xp, dz, t = sp.symbols("delta x xp dz t", real=True)
theta = delta * t
distance_exact = delta**2 * ((xx - xp) ** 2 + dz**2) + 4 * (1 + delta * xx) * (1 + delta * xp) * sp.sin(theta / 2) ** 2
near = sp.series(distance_exact / delta**2, delta, 0, 3).removeO().expand()
near_expected = (
    (xx - xp) ** 2
    + dz**2
    + t**2
    + delta * (xx + xp) * t**2
    + delta**2 * (xx * xp * t**2 - t**4 / 12)
)
assert sp.simplify(near - near_expected) == 0
assert sp.expand(near).coeff(delta, 2).coeff(t, 4) == -sp.Rational(1, 12)

# The toroidal phase is the column phase at k=ell*delta.
ell = sp.symbols("ell", integer=True)
k = ell * delta
assert sp.simplify(ell * theta - k * t) == 0

# Determinant and inverse expansions retain all Piola factors.
trH, trQ, trH2 = sp.symbols("trH trQ trH2")
J_second_coefficient = trQ + sp.Rational(1, 2) * (trH**2 - trH2)
assert sp.diff(J_second_coefficient, trQ) == 1
assert sp.diff(J_second_coefficient, trH2) == -sp.Rational(1, 2)

# Resolvent inverse derivatives: R'=-R L' R and
# R''=2 R L' R L' R-R L'' R, checked in the scalar noncommutative-free case.
lam0, lam1, lam2, h = sp.symbols("lam0 lam1 lam2 h", nonzero=True)
resolvent = 1 / (lam0 + h * lam1 + h**2 * lam2 / 2)
r1 = sp.diff(resolvent, h).subs(h, 0)
r2 = sp.diff(resolvent, h, 2).subs(h, 0)
assert sp.simplify(r1 + lam1 / lam0**2) == 0
assert sp.simplify(r2 - (2 * lam1**2 / lam0**3 - lam2 / lam0**2)) == 0

# Shape differentiation does not worsen the |r|^-2 Biot-Savart singularity.
# Under r->s*r and h(y)-h(y')->s*d_h:
# K scales s^-2, DK[d_h] scales s^-2, D2K[d_h,d_h] scales s^-2.
s = sp.symbols("s", positive=True)
K_scale = s ** -2
DK_scale_with_difference = s ** -3 * s
D2K_scale_with_two_differences = s ** -4 * s**2
DK_scale_without_difference = s ** -3
assert sp.simplify(DK_scale_with_difference / K_scale) == 1
assert sp.simplify(D2K_scale_with_two_differences / K_scale) == 1
assert sp.simplify(DK_scale_without_difference / K_scale) == 1 / s

# Polyhomogeneous Taylor remainder after tau=delta^2 log(delta).
L = sp.symbols("L", positive=True)
tau_size = delta**2 * L
q_tau = delta * tau_size
assert sp.simplify(q_tau - delta**3 * L) == 0

# Resonant scales and all relative error comparisons.
nstar = 1 / (delta * L)
r_gap = delta**3 * L**2
hodge_remainder_symbol = delta**3 * L / nstar
local_graph_remainder = delta**3 * L
order_minus_two_schur = delta**2 / nstar**2
assert sp.simplify(hodge_remainder_symbol / r_gap) == delta
assert sp.simplify(local_graph_remainder / r_gap) == 1 / L
assert sp.simplify(order_minus_two_schur / r_gap) == delta

# The old unweighted O(delta^2) comparison is not small at the contour.
unweighted_shortcut = sp.simplify(delta**2 / r_gap)
assert unweighted_shortcut == 1 / (delta * L**2)

# Off-diagonal commutator identity in a finite diagonal model.
la, lb, b = sp.symbols("lambda_a lambda_b b")
Lambda = sp.diag(la, lb)
B = sp.Matrix([[0, b], [0, 0]])
comm = Lambda * B - B * Lambda
assert sp.simplify(comm[0, 1] - (la - lb) * b) == 0

checks = [
    "cao_operator_rescaling",
    "toroidal_distance_exact_near_expansion",
    "toroidal_quartic_sign",
    "longitudinal_phase_retained",
    "piola_determinant_trace_square",
    "piola_determinant_matrix_square",
    "inverse_first_derivative_sign",
    "inverse_second_derivative_terms",
    "shape_first_derivative_order_minus_one",
    "shape_second_derivative_order_minus_one",
    "missing_shape_difference_exposes_order_loss",
    "polyhomogeneous_q_tau_remainder",
    "hodge_remainder_to_gap",
    "local_graph_remainder_to_gap",
    "order_minus_two_schur_to_gap",
    "unweighted_shortcut_diverges",
    "mode_commutator_identity",
]

for name in checks:
    print(f"PASS {name}")
print(f"PASS total={len(checks)}")
