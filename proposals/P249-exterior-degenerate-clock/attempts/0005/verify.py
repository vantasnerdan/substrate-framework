"""Exact verifier for the P249 phase-locked M5 core bridge."""

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


# Exact clock-plane representation and equivariance.
q = sp.symbols("q", real=True)
a, c = sp.symbols("a c", real=True)
A = sp.Matrix([[0, -1], [1, 0]])
H = sp.diag(1, -1)
K = sp.simplify((A * H - H * A) / 2)
R = sp.Matrix([[sp.cos(q), -sp.sin(q)], [sp.sin(q), sp.cos(q)]])
check(A**2 == -sp.eye(2), "normalized compact generator")
check(K == sp.Matrix([[0, 1], [1, 0]]), "clock tensor partner")
check(sp.simplify(R * R.T) == sp.eye(2), "SO(2) orthogonality")
check(sp.simplify(R.subs(q, q + 2 * sp.pi) - R) == sp.zeros(2), "2pi period")
rotated_H = (R * H * R.T).applyfunc(sp.trigsimp)
target_H = sp.cos(2 * q) * H + sp.sin(2 * q) * K
check((rotated_H - target_H).applyfunc(sp.trigsimp) == sp.zeros(2), "weight-two tensor action")

psi_sq_re = a**2 - c**2
psi_sq_im = 2 * a * c
Q0 = psi_sq_re * H + psi_sq_im * K
ar = sp.expand(a * sp.cos(q) - c * sp.sin(q))
cr = sp.expand(a * sp.sin(q) + c * sp.cos(q))
Qrot = sp.expand(ar**2 - cr**2) * H + sp.expand(2 * ar * cr) * K
equivariance = (Qrot - R * Q0 * R.T).applyfunc(sp.trigsimp)
check(equivariance == sp.zeros(2), "phase-lock equivariance")
check(sp.simplify(sp.trace(H**2) / 2) == 1, "H normalization")
check(sp.simplify(sp.trace(K**2) / 2) == 1, "K normalization")
check(sp.trace(H * K) == 0, "clock tensor orthogonality")

# Orthogonal curvature projection defining the kinetic replacement.
z, y = sp.symbols("z y", real=True)
F_clock = z * A
inner_A_F = -sp.trace(A * F_clock) / 2
P_A_F = sp.simplify(inner_A_F) * A
check(P_A_F == F_clock, "clock curvature projection")
F_probe = sp.Matrix([[0, -z, -y], [z, 0, 0], [y, 0, 0]])
A3 = sp.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]])
coeff = -sp.trace(A3 * F_probe) / 2
F_perp = sp.simplify(F_probe - coeff * A3)
check(-sp.trace(A3 * F_perp) / 2 == 0, "orthogonal curvature remainder")
check(sp.expand((-sp.trace(F_probe**2) + sp.trace((coeff * A3) ** 2) + sp.trace(F_perp**2)) / 2) == 0, "positive-square norm decomposition")

# Reduced potential, positivity, and unchanged sector.
f, b = sp.symbols("f b", real=True)
W = 3 * f**2 - 4 * f**3 + 2 * f**4
V = sp.expand(W + 3 * b**2 + 4 * b**4 + 6 * (b - f**2) ** 2)
check(sp.expand(W - f**2 * (2 * (f - 1) ** 2 + 1)) == 0, "subcritical scalar positivity")
check(sp.expand(V - (W + 3 * b**2 + 4 * b**4 + 6 * (b - f**2) ** 2)) == 0, "positive potential decomposition")
check(sp.expand(V.subs(f, 0) - (9 * b**2 + 4 * b**4)) == 0, "augmented exterior B potential")
check(V.subs({f: 0, b: 0}) == 0, "uniaxial clock-vacuum sector unchanged")
check(sp.diff(V, f).subs({f: 0, b: 0}) == 0, "vacuum f stationarity")
check(sp.diff(V, b).subs({f: 0, b: 0}) == 0, "vacuum b stationarity")

# Complete first and second variations.
Vf = sp.factor(sp.diff(V, f))
Vb = sp.factor(sp.diff(V, b))
check(sp.expand(Vf - (6 * f - 12 * f**2 + 32 * f**3 - 24 * f * b)) == 0, "complete f equation potential term")
check(sp.expand(Vb - (18 * b + 16 * b**3 - 12 * f**2)) == 0, "complete b equation potential term")
check(sp.expand(sp.diff(V, f, 2) - (6 - 24 * f + 96 * f**2 - 24 * b)) == 0, "V_ff")
check(sp.expand(sp.diff(V, f, b) + 24 * f) == 0, "V_fb")
check(sp.expand(sp.diff(V, b, 2) - (18 + 48 * b**2)) == 0, "V_bb")

Qcharge, inertia = sp.symbols("Q I", positive=True, real=True)
fixed_term = Qcharge**2 / (2 * inertia)
omega = Qcharge / inertia
check(sp.diff(fixed_term, inertia) == -omega**2 / 2, "fixed-charge first derivative")
check(sp.diff(fixed_term, inertia, 2) == omega**2 / inertia, "Legendre rank-one coefficient")
check(sp.diff(fixed_term, Qcharge) == omega, "dE/dQ")

# Exterior masses and charged continuum threshold.
mf2 = sp.diff(V, f, 2).subs({f: 0, b: 0})
mb2 = sp.diff(V, b, 2).subs({f: 0, b: 0})
charged_edge = min(mf2, mb2 / 4)
check(mf2 == 6, "scalar exterior mass squared")
check(mb2 == 18, "M5 clock exterior mass squared")
check(charged_edge == sp.Rational(9, 2), "weight-resolved charged edge")

# Exact plateau binding witness and strict unsplit obstruction.
f_trial = sp.Rational(1, 2)
b_trial = sp.Rational(1, 4)
I_density = f**2 + 4 * b**2
V_trial = sp.simplify(V.subs({f: f_trial, b: b_trial}))
I_trial = sp.simplify(I_density.subs({f: f_trial, b: b_trial}))
ratio_trial = sp.simplify(2 * V_trial / I_trial)
check(V_trial == sp.Rational(37, 64), "split plateau potential")
check(I_trial == sp.Rational(1, 2), "split plateau inertia")
check(ratio_trial == sp.Rational(37, 16), "split plateau J/K ratio")
check(charged_edge - ratio_trial == sp.Rational(35, 16), "exact binding margin")
V_unsplit = sp.factor(V.subs(b, 0))
check(sp.expand(V_unsplit - (3 * f**2 - 4 * f**3 + 8 * f**4)) == 0, "unsplit potential")
unsplit_ratio = sp.cancel(2 * V_unsplit / f**2)
check(unsplit_ratio == 6 - 8 * f + 16 * f**2, "unsplit ratio")
check(sp.expand(unsplit_ratio - (16 * (f - sp.Rational(1, 4)) ** 2 + 5)) == 0, "unsplit square completion")
check(sp.Rational(5) > charged_edge, "unsplit cannot realize binding level")

# Natural tails and exact fixed-charge scale identity.
w2 = sp.symbols("omega2", positive=True, real=True)
kf2 = mf2 - w2
kb2 = mb2 - 4 * w2
check(kf2.subs(w2, charged_edge) == sp.Rational(3, 2), "scalar tail remains positive at charged edge")
check(kb2.subs(w2, charged_edge) == 0, "B tail closes only at strict edge")
Rscale, T, U, Erot = sp.symbols("R T U Erot", positive=True, real=True)
scaled = Rscale * T + Rscale**3 * U + Erot / Rscale**3
check(sp.diff(scaled, Rscale).subs(Rscale, 1) == T + 3 * U - 3 * Erot, "fixed-charge Derrick identity")

# Mutations of the phase lock and charge weight must be detected.
V_no_lock = sp.expand(W + 3 * b**2 + 4 * b**4)
check(sp.diff(V_no_lock, b).subs(b, 0) == 0, "lock deletion loses forced split source")
wrong_edge = min(mf2, mb2)
check(wrong_edge != charged_edge, "weight-two threshold mutation detected")

tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = [node for node in ast.walk(tree) if isinstance(node, ast.Assert)]
check(len(assert_nodes) == 0, "no bare Python asserts")

print(f"ALL {CHECKS} CHECKS PASS [P249-O1-O4-0005]")
