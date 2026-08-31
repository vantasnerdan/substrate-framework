"""Exact verifier for the P249 canonical auxiliary-amplitude seed."""

from __future__ import annotations

import ast
from pathlib import Path

import sympy as sp


CHECKS = 0


def check(condition: object, message: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.true:
        raise AssertionError(message)
    CHECKS += 1


f, x, q, n = sp.symbols("f x q n", positive=True, real=True)
omega = sp.symbols("Omega", positive=True, real=True)
W = 3 * f**2 - 6 * f**4 + 4 * f**6
mass_sq = sp.diff(W, f, 2).subs(f, 0)
nonlinear = sp.expand(W - mass_sq * f**2 / 2)

# W-i and W-ii.
positive_factor = sp.expand(f**2 * (4 * (f**2 - sp.Rational(3, 4)) ** 2 + sp.Rational(3, 4)))
check(sp.expand(W - positive_factor) == 0, "positive factorization")
check(sp.Rational(3, 4) > 0, "strict positive bracket margin")
check(W.subs(f, 0) == 0, "vacuum value")
check(sp.diff(W, f).subs(f, 0) == 0, "vacuum stationarity")
check(mass_sq == 6, "positive asymptotic mass")

# W-iii and W-iiii(b).
f0 = sp.sqrt(3) / 2
f1 = sp.Integer(2)
check(sp.simplify(nonlinear - (-6 * f**4 + 4 * f**6)) == 0, "nonlinear remainder")
check(sp.simplify(nonlinear.subs(f, f0)) == -sp.Rational(27, 16), "negative hylomorphy witness")
check(sp.simplify(sp.diff(nonlinear, f) - 24 * f**3 * (f**2 - 1)) == 0, "remainder derivative")
check(f1 > f0, "growth witness ordering")
check(sp.diff(nonlinear, f).subs(f, f1) > 0, "W-iiii(b) derivative witness")

# Exact frequency window and an interior frequency.
ratio = sp.cancel(2 * W / f**2)
ratio_square = sp.expand(sp.Rational(3, 2) + 8 * (f**2 - sp.Rational(3, 4)) ** 2)
check(sp.expand(ratio - ratio_square) == 0, "frequency-ratio square completion")
omega0_sq = sp.Rational(3, 2)
chosen_omega_sq = sp.Integer(4)
check(omega0_sq < chosen_omega_sq < mass_sq, "Omega=2 is strictly inside the existence interval")
effective = sp.expand(W - chosen_omega_sq * f**2 / 2)
check(sp.expand(effective - f**2 * (1 - 6 * f**2 + 4 * f**4)) == 0, "effective potential")
roots = sp.solve(4 * x**2 - 6 * x + 1, x)
check(roots == [(sp.Integer(3) - sp.sqrt(5)) / 4, (sp.Integer(3) + sp.sqrt(5)) / 4], "negative-well edges")
check(effective.subs(f, f0) < 0, "exact negative effective-potential witness")

# Noether/Legendre normalization and fixed-Q derivatives, represented by N=n.
charge = omega * n
energy_rot = sp.Rational(1, 2) * omega**2 * n
check(sp.diff(energy_rot, omega) == charge, "regular phase Legendre map")
fixed_charge_term = q**2 / (2 * n)
omega_fixed = q / n
check(sp.diff(fixed_charge_term, n) == -omega_fixed**2 / 2, "fixed-Q first derivative")
check(sp.diff(fixed_charge_term, n, 2) == omega_fixed**2 / n, "fixed-Q positive rank-one coefficient")
check(sp.diff(fixed_charge_term, q) == omega_fixed, "dE/dQ envelope derivative")

# Tail and the fixed-Q scale virial.
tail_sq = sp.simplify(mass_sq - chosen_omega_sq)
check(tail_sq == 2, "action-derived positive tail mass")
R, T, V, Erot = sp.symbols("R T V Erot", positive=True, real=True)
scaled = R * T + R**3 * V + Erot / R**3
check(sp.diff(scaled, R).subs(R, 1) == T + 3 * V - 3 * Erot, "complete fixed-Q scale virial")

# Defining-object mutations: removing the sextic destroys boundedness and
# removing Legendre curvature destroys the rank-one term.
mutated_ratio = sp.cancel(2 * (3 * f**2 - 6 * f**4) / f**2)
check(sp.limit(mutated_ratio, f, sp.oo) == -sp.oo, "sextic boundedness mutation detected")
check(sp.diff(q**2 / (2 * n), n, 2) != 0, "rank-one deletion mutation detected")

# The verifier itself must remain assertion-free; every predicate goes through
# the counted check function and a terminal tally.
tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = [node for node in ast.walk(tree) if isinstance(node, ast.Assert)]
check(len(assert_nodes) == 0, "no bare Python asserts")

print(f"ALL {CHECKS} CHECKS PASS [P249-O3-0004]")
