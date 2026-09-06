#!/usr/bin/env python3
"""Exact algebraic checks for P253/0052."""

import sympy as sp


# Correct two-scale resonance inversion.
h, delta = sp.symbols("h delta", positive=True)
A, B, c, d, beta = sp.symbols("A B c d beta", nonzero=True)
exact_leading = A / (delta * (c / h + d))
series = sp.series(exact_leading, h, 0, 3).removeO().expand()
expected = A * h / (c * delta) - A * d * h**2 / (c**2 * delta)
assert sp.simplify(series - expected) == 0
nstar = expected + B / A - beta
assert sp.expand(nstar).coeff(h, 1) == A / (c * delta)
assert sp.expand(nstar).coeff(h, 2) == -A * d / (c**2 * delta)
assert sp.expand(nstar).subs(h, 0) == B / A - beta

# Preserved finite-window scales at n*=Theta(1/(delta*L)).
L = sp.symbols("L", positive=True)
n0 = 1 / (delta * L)
n = sp.symbols("n", positive=True)
coupling = delta**2 / n
spacing = delta / n**2
assert sp.simplify(coupling.subs(n, n0)) == delta**3 * L
assert sp.simplify(spacing.subs(n, n0)) == delta**3 * L**2
assert sp.simplify((coupling / spacing).subs(n, n0)) == 1 / L

# Action coordinate: I=Ia*s^2 gives constant Jacobian relative to s ds dbeta.
Ia, s = sp.symbols("I_a s", positive=True)
I = Ia * s**2
assert sp.diff(I, s) == 2 * Ia * s
assert sp.simplify(sp.diff(I, s) / s) == 2 * Ia

# Contravariant Piola identity for a diagonal affine map. This exposes the
# determinant factor; the general identity follows by the same cofactor law.
aa, bb, cc = sp.symbols("a b c", positive=True)
X, Y, Z = sp.symbols("X Y Z")
u1 = sp.Function("u1")(X, Y, Z)
u2 = sp.Function("u2")(X, Y, Z)
u3 = sp.Function("u3")(X, Y, Z)
Jac = aa * bb * cc
div_physical = (
    sp.diff(aa * u1 / Jac, X) / aa
    + sp.diff(bb * u2 / Jac, Y) / bb
    + sp.diff(cc * u3 / Jac, Z) / cc
)
div_reference = (sp.diff(u1, X) + sp.diff(u2, Y) + sp.diff(u3, Z)) / Jac
assert sp.simplify(div_physical - div_reference) == 0

# Canonical KKS sign/factor for i_X Omega=dH.
nu = sp.symbols("nu", positive=True, real=True)
Omega = sp.Matrix([[0, 1], [-1, 0]])
Hess = sp.diag(nu**2, 1)
qplus = sp.Matrix([1, sp.I * nu])
qbar = sp.conjugate(qplus)
omega_pair = (qplus.T * Omega * qbar)[0]
hess_pair = (qplus.T * Hess * qbar)[0]
assert sp.simplify(omega_pair + 2 * sp.I * nu) == 0
assert sp.simplify(hess_pair - 2 * nu**2) == 0
assert sp.simplify(hess_pair - sp.I * nu * omega_pair) == 0

# Inverse of the local triangular transport/shear block N^2=0.
z = sp.symbols("z", nonzero=True)
N = sp.Matrix([[0, 1], [0, 0]])
Id = sp.eye(2)
candidate_inverse = Id / z - N / z**2
assert (z * Id + N) * candidate_inverse == Id
assert candidate_inverse * (z * Id + N) == Id

# HJ2 modewise remainders are smaller than the resonant spacing.
eps = sp.symbols("epsilon", positive=True)
r_gap = delta**3 * L**2
tail_to_translation = eps * delta**2 / n
tail_internal = eps * delta**2 / n**2
assert sp.simplify((tail_to_translation / r_gap).subs(n, n0)) == eps / L
assert sp.simplify((tail_internal / r_gap).subs(n, n0)) == delta * eps

# Elliptic-integral Green expansion through rho*log(rho).
rho = sp.symbols("rho", positive=True)
ell = sp.symbols("ell", real=True)  # ell = log(4/kappa_prime)
s2 = sp.symbols("s2", positive=True)  # kappa_prime^2
K = ell + s2 * (ell - 1) / 4
E = 1 + s2 * (ell - sp.Rational(1, 2)) / 2
kappa = 1 - s2 / 2
inv_kappa = 1 + s2 / 2
bracket = sp.expand((2 * inv_kappa - kappa) * K - 2 * inv_kappa * E)
bracket = sp.series(bracket, s2, 0, 2).removeO().expand()
assert sp.simplify(bracket.coeff(s2, 0) - (ell - 2)) == 0
assert sp.simplify(bracket.coeff(s2, 1) - sp.Rational(3, 4) * (ell - 1)) == 0
logrho, log8 = sp.symbols("logrho log8", real=True)
ell0 = log8 + logrho / 2
# s2=rho/4+O(rho^2), while ell=ell0+rho/8+O(rho^2).
rho_coefficient_in_source_bracket = sp.expand(
    2 * (sp.Rational(1, 8) + sp.Rational(3, 16) * (ell0 - 1))
)
expected_rho_coefficient = (
    sp.Rational(3, 16) * logrho
    + sp.Rational(3, 8) * log8
    - sp.Rational(1, 8)
)
assert sp.simplify(rho_coefficient_in_source_bracket - expected_rho_coefficient) == 0

checks = [
    "resonance_leading_inverse",
    "resonance_finite_bending_correction",
    "resonance_phase_correction",
    "coupling_scale_preserved",
    "spacing_scale_preserved",
    "window_ratio_preserved",
    "action_angle_constant_jacobian",
    "piola_divergence_factor",
    "kks_canonical_sign",
    "kks_hessian_identity",
    "triangular_transport_left_inverse",
    "triangular_transport_right_inverse",
    "HJ2_translation_remainder_ratio",
    "HJ2_internal_remainder_ratio",
    "green_elliptic_leading_term",
    "green_elliptic_rho_log_term",
    "green_rho_finite_coefficient",
]
for name in checks:
    print(f"PASS {name}")
print(f"PASS total={len(checks)}")
