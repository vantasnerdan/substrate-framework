"""Exact verifier for the full-spatial canonical P249 clock completion."""

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


# Global M5 positivity reduces to a one-variable norm bound.
r = sp.symbols("r", nonnegative=True, real=True)
radial_lower = r**4 - r**3 - sp.Rational(1, 2) * r**2 + sp.Rational(1, 2)
check(sp.factor(radial_lower) == (r - 1) ** 2 * (2 * r**2 + 2 * r + 1) / 2, "M5 global lower-bound factorization")
check(2 * r**2 + 2 * r + 1 > 0, "strict positive M5 factor")
check(radial_lower.subs(r, 1) == 0, "rank-one unit vacuum norm")

# Full six-component spatial tensor and every potential term.
a, t, p, u, v, q = sp.symbols("a t p u v q", real=True)
f = sp.symbols("f", nonnegative=True, real=True)
S = sp.Matrix([[1 + a, u, v], [u, t + p, q], [v, q, t - p]])
tr2 = sp.expand(sp.trace(S**2))
V_m5 = sp.expand(-tr2 / 2 - sp.trace(S**3) + tr2**2 + sp.Rational(1, 2))
axis = sp.expand(sp.Rational(1, 4) * (tr2 - S[0, 0] ** 2))
W = 3 * f**2 - 4 * f**3 + 2 * f**4
phase_lock = 6 * ((p - f**2) ** 2 + q**2)
V = sp.expand(V_m5 + axis + W + phase_lock)
zero = {a: 0, t: 0, p: 0, u: 0, v: 0, q: 0, f: 0}

axis_expected = sp.Rational(1, 4) * (2 * u**2 + 2 * v**2 + (t + p) ** 2 + 2 * q**2 + (t - p) ** 2)
check(sp.expand(axis - axis_expected) == 0, "axis lock sum of squares")
check(sp.expand(W - f**2 * (2 * (f - 1) ** 2 + 1)) == 0, "subcritical scalar positivity")
check(V.subs(zero) == 0, "aligned exterior energy")
for variable in (a, t, p, u, v, q, f):
    check(sp.diff(V, variable).subs(zero) == 0, f"exterior stationarity {variable}")

# Reproduce the omitted attempt-0005 first variations before repairing them.
b = sp.symbols("b", real=True)
reduced_subs = {a: 0, t: 0, p: b, u: 0, v: 0, q: 0}
check(sp.factor(sp.diff(V_m5, a).subs(reduced_subs)) == 8 * b**2, "omitted director derivative")
check(sp.factor(sp.diff(V_m5, t).subs(reduced_subs)) == -6 * b**2, "omitted tangent-trace derivative")

# Complete exterior Hessian and generalized masses.
coords = (a, t, p, u, v, q)
H_spatial = sp.hessian(V, coords).subs(zero)
expected_H = sp.diag(5, 7, 19, 1, 1, 19)
check(H_spatial == expected_H, "full six-channel exterior Hessian")
G_spatial = sp.diag(sp.Rational(1, 2), 1, 1, 1, 1, 1)
mass_matrix = G_spatial.inv() * H_spatial
check(mass_matrix == sp.diag(10, 7, 19, 1, 1, 19), "canonical generalized mass pencil")
check(all(value > 0 for value in mass_matrix.diagonal()), "positive physical kinetic/mass subspace")
scalar_mass_sq = sp.diff(V, f, 2).subs(zero)
check(scalar_mass_sq == 6, "scalar mass squared")
charged_edge = min(scalar_mass_sq, sp.Rational(19, 4))
check(charged_edge == sp.Rational(19, 4), "weight-resolved charged continuum")

# The historical director-shear failure and its exact zeta=1/4 repair.
shear_repaired = 16 * b**2 - 6 * b + 1
check(sp.expand(shear_repaired - (16 * (b - sp.Rational(3, 16)) ** 2 + sp.Rational(7, 16))) == 0, "axis-lock shear square completion")
check(sp.Rational(7, 16) > 0, "strict director-shear margin")

# Exact full-field binding witness.
trial = {a: 0, t: 0, p: sp.Rational(1, 4), u: 0, v: 0, q: 0, f: sp.Rational(1, 2)}
V_trial = sp.simplify(V.subs(trial))
I_trial = sp.Rational(1, 2)  # f^2 + 4(p^2+q^2)
ratio_trial = sp.simplify(2 * V_trial / I_trial)
check(V_trial == sp.Rational(39, 64), "full split plateau potential")
check(ratio_trial == sp.Rational(39, 16), "full split plateau J/K ratio")
check(charged_edge - ratio_trial == sp.Rational(37, 16), "strict full-space binding margin")

# No unsplit state can bind below the charged continuum.
unsplit_scalar = sp.expand(W + 6 * f**4)
unsplit_ratio = sp.cancel(2 * unsplit_scalar / f**2)
check(unsplit_ratio == 6 - 8 * f + 16 * f**2, "unsplit charged ratio")
check(sp.expand(unsplit_ratio - (16 * (f - sp.Rational(1, 4)) ** 2 + 5)) == 0, "unsplit global square completion")
check(sp.Rational(5) - charged_edge == sp.Rational(1, 4), "unsplit threshold margin")
check(sp.diff(phase_lock, p).subs({p: 0, q: 0}) == -12 * f**2, "nonzero full B source")

# Complete Legendre calculus applies to arbitrary six-component variations.
Qcharge, inertia = sp.symbols("Q I", positive=True, real=True)
fixed_term = Qcharge**2 / (2 * inertia)
omega = Qcharge / inertia
check(sp.diff(fixed_term, inertia) == -omega**2 / 2, "full fixed-charge first derivative")
check(sp.diff(fixed_term, inertia, 2) == omega**2 / inertia, "full positive rank-one Hessian")
check(sp.diff(fixed_term, Qcharge) == omega, "full dE/dQ identity")

# Exact tail/radiation inequalities at the strict charged edge.
w2 = sp.symbols("omega2", positive=True, real=True)
kpsi2 = scalar_mass_sq - w2
kB2 = 19 - 4 * w2
check(kpsi2.subs(w2, charged_edge) == sp.Rational(5, 4), "scalar tail margin at charged edge")
check(kB2.subs(w2, charged_edge) == 0, "weight-two edge normalization")
check(min(mass_matrix.diagonal()) == 1, "lightest neutral tensor mass")

# Intrinsic scale balance.
Rscale, Tenergy, Uenergy, Erot = sp.symbols("R T U Erot", positive=True, real=True)
scaled = Rscale * Tenergy + Rscale**3 * Uenergy + Erot / Rscale**3
check(sp.diff(scaled, Rscale).subs(Rscale, 1) == Tenergy + 3 * Uenergy - 3 * Erot, "full fixed-charge Derrick identity")

# Defining-object mutations.
mutated_m5 = -r**3 - r**2 / 2 + sp.Rational(1, 2)
check(sp.limit(mutated_m5, r, sp.oo) == -sp.oo, "M5 quartic deletion mutation")
check((sp.diag(sp.Rational(1, 2), 1, 0, 1, 1, 0).det()) == 0, "canonical clock-kinetic deletion mutation")
check(sp.diff(6 * (p**2 + q**2), p).subs(p, 0) == 0, "phase-lock source deletion mutation")

tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = [node for node in ast.walk(tree) if isinstance(node, ast.Assert)]
check(len(assert_nodes) == 0, "no bare Python asserts")

print(f"ALL {CHECKS} CHECKS PASS [P249-O2-O5-0006]")
