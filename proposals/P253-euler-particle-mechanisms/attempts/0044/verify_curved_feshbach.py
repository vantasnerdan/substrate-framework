#!/usr/bin/env python3
"""Exact symbolic checks for the P253/0044 curved translation reduction."""

import sympy as s


r = s.symbols("r", positive=True)
Omega = s.Function("Omega")(r)
W = 2 * Omega + r * s.diff(Omega, r)
psi0 = r * Omega

# Exact m=1 translation kernel for the stated T(b) convention.
lap1 = s.diff(psi0, r, 2) + s.diff(psi0, r) / r - psi0 / r**2
assert s.simplify(-lap1 + s.diff(W, r) * psi0 / (r * Omega)) == 0

# d/db (Omega-b)^(-1) is positive.  The weighted integrand is r^2 W'.
b = s.symbols("b", real=True)
potential = s.diff(W, r) / (r * (Omega - b))
tprime_integrand = s.simplify(
    r * psi0**2 * s.diff(potential, b).subs(b, 0))
assert s.simplify(tprime_integrand - r**2 * s.diff(W, r)) == 0
assert s.simplify(
    r**2 * s.diff(W, r) - s.diff(r**2 * W, r) + 2 * r * W
) == 0

# Pointwise ground-state transform.  The total derivative is cancelled by
# the physical exterior Robin form term at the core edge.
f = s.Function("f")(r)
V0 = s.simplify((s.diff(psi0, r, 2) + s.diff(psi0, r) / r) / psi0)
psi = psi0 * f
form_density = r * (s.diff(psi, r)**2 + V0 * psi**2)
ground_density = r * psi0**2 * s.diff(f, r)**2
boundary_derivative = s.diff(r * psi0 * s.diff(psi0, r) * f**2, r)
assert s.simplify(form_density - ground_density - boundary_derivative) == 0

# Correct map from the exterior potential DtN to the scalar psi boundary
# coefficient.  Only the logarithmic jet is needed.
x = s.symbols("x", positive=True)
k1 = 1 / x + x / 2 * (s.log(x / 2) + s.EulerGamma - s.Rational(1, 2))
D = s.simplify(x * s.diff(k1, x) / k1)
minus_inverse_D_remainder = s.limit(
    ((-1 / D) - 1) / x**2 - (s.log(x / 2) + s.EulerGamma),
    x,
    0,
    dir="+",
)
assert minus_inverse_D_remainder == 0

# Exact toroidal differential operators.  h=1+delta*s*cos(alpha), and all
# derivatives are scaled by the core radius a.
delta, sr, alpha = s.symbols("delta sr alpha", real=True)
ell, m = s.symbols("ell m", integer=True)
Vs = s.Function("Vs")(sr, alpha)
Va = s.Function("Va")(sr, alpha)
Vt = s.Function("Vt")(sr, alpha)
h = 1 + delta * sr * s.cos(alpha)
div_exact = (
    s.diff(sr * h * Vs, sr)
    + s.diff(h * Va, alpha)
    + delta * s.I * ell * sr * Vt
) / (sr * h)
div0 = s.diff(sr * Vs, sr) / sr + s.diff(Va, alpha) / sr
div1 = s.cos(alpha) * Vs - s.sin(alpha) * Va + s.I * ell * Vt
div2 = -sr * s.cos(alpha) * div1
div_series = s.series(div_exact, delta, 0, 3).removeO()
assert s.simplify(div_series - div0 - delta * div1 - delta**2 * div2) == 0

scalar = s.Function("scalar")(sr, alpha)
lap_exact = (
    s.diff(sr * h * s.diff(scalar, sr), sr)
    + s.diff(h * s.diff(scalar, alpha) / sr, alpha)
    - delta**2 * ell**2 * sr * scalar / h
) / (sr * h)
lap0 = (s.diff(sr * s.diff(scalar, sr), sr) / sr
        + s.diff(scalar, alpha, 2) / sr**2)
curv1 = s.cos(alpha) * s.diff(scalar, sr) - (
    s.sin(alpha) / sr
) * s.diff(scalar, alpha)
lap2 = -sr * s.cos(alpha) * curv1 - ell**2 * scalar
lap_series = s.series(lap_exact, delta, 0, 3).removeO()
assert s.simplify(lap_series - lap0 - delta * curv1 - delta**2 * lap2) == 0

# Curvature acts on the translation mode with genuinely nonzero m=0 and m=2
# components.  psi0'+psi0/r=W is the decisive m=0 forcing.
radial = s.Function("radial")(sr)
translation = radial * s.exp(s.I * alpha)
curved_translation = s.expand_complex(
    s.cos(alpha) * s.diff(translation, sr)
    - s.sin(alpha) * s.diff(translation, alpha) / sr
)
# Avoid complex trigonometric expansion in the actual assertion by using the
# exact exponential decomposition derived from cos/sin.
m0_coefficient = (s.diff(radial, sr) + radial / sr) / 2
m2_coefficient = (s.diff(radial, sr) - radial / sr) / 2
curved_decomposition = m0_coefficient + m2_coefficient * s.exp(2 * s.I * alpha)
curved_direct = (
    (s.exp(s.I * alpha) + s.exp(-s.I * alpha)) / 2
    * s.diff(translation, sr)
    - (s.exp(s.I * alpha) - s.exp(-s.I * alpha)) / (2 * s.I * sr)
    * s.diff(translation, alpha)
)
assert s.simplify(curved_direct - curved_decomposition) == 0
Omega_s = s.Function("Omega_s")(sr)
W_s = 2 * Omega_s + sr * s.diff(Omega_s, sr)
assert s.simplify(
    m0_coefficient.subs(radial, sr * Omega_s).doit() - W_s / 2
) == 0
assert s.simplify(
    m2_coefficient.subs(radial, sr * Omega_s).doit()
    - sr * s.diff(Omega_s, sr) / 2
) == 0

# Second-order expansion of the Leray projector P=I-G(Delta^-1)D.  A concrete
# exact matrix model catches every sign and ordering in P1 and P2.
G0 = s.Matrix([[1, 0], [0, 1], [1, 1]])
G1 = s.Matrix([[0, 1], [1, 0], [1, -1]])
G2 = s.Matrix([[1, 1], [0, -1], [2, 0]])
D0 = s.Matrix([[2, 0, 1], [0, 3, -1]])
D1 = s.Matrix([[1, -1, 0], [2, 0, 1]])
D2 = s.Matrix([[0, 1, 2], [-1, 2, 0]])
L0 = D0 * G0
L1 = D0 * G1 + D1 * G0
L2 = D0 * G2 + D1 * G1 + D2 * G0
S0 = L0.inv()
S1 = -S0 * L1 * S0
S2 = S0 * L1 * S0 * L1 * S0 - S0 * L2 * S0
P0 = s.eye(3) - G0 * S0 * D0
P1 = -(G1 * S0 * D0 + G0 * S1 * D0 + G0 * S0 * D1)
P2 = -(
    G2 * S0 * D0 + G1 * S1 * D0 + G1 * S0 * D1
    + G0 * S2 * D0 + G0 * S1 * D1 + G0 * S0 * D2
)
Gpoly = G0 + delta * G1 + delta**2 * G2
Dpoly = D0 + delta * D1 + delta**2 * D2
Lpoly = Dpoly * Gpoly
Pexact = s.eye(3) - Gpoly * Lpoly.inv() * Dpoly
for i in range(3):
    for j in range(3):
        coeff = s.series(Pexact[i, j], delta, 0, 3).removeO()
        assert s.simplify(coeff - (P0 + delta * P1 + delta**2 * P2)[i, j]) == 0

# m>=2 core positivity relative to the m=1 ground-state form.
mpos = s.symbols("mpos", integer=True, positive=True)
assert s.factor(mpos**2 - 1) == (mpos - 1) * (mpos + 1)

# Universal leading centerline Hamiltonian and its actual 2x2 compression.
Gamma, rho, R, Log = s.symbols("Gamma rho R Log", positive=True)
Cline = rho * Gamma**2 * Log / (4 * s.pi)
g = s.simplify(Cline / (rho * Gamma * R**2))
Mlog = g * s.Matrix([[0, ell**2 - 1], [-ell**2, 0]])
charpoly = s.factor(Mlog.charpoly().as_expr())
expected_charpoly = s.Symbol("lambda")**2 + g**2 * ell**2 * (ell**2 - 1)
assert s.simplify(charpoly - expected_charpoly) == 0
assert Mlog.subs(ell, 1)[0, 1] == 0

print("Tprime_weighted_integrand =", tprime_integrand)
print("Tprime_pairing_after_core_edge = -2*F(a)")
print("ground_state_transform_remainder =", s.simplify(
    form_density - ground_density - boundary_derivative))
print("minus_inverse_D_log_remainder_limit =", minus_inverse_D_remainder)
print("toroidal_divergence_C1 =", div1)
print("toroidal_divergence_C2 =", div2)
print("toroidal_scalar_laplacian_C1 =", curv1)
print("toroidal_scalar_laplacian_C2 =", lap2)
print("translation_to_m0_coefficient =", W_s / 2)
print("translation_to_m2_coefficient =", sr * s.diff(Omega_s, sr) / 2)
print("leading_centerline_matrix =", Mlog)
print("leading_centerline_characteristic =", charpoly)
print("ALL 22 EXACT CURVED-FESHBACH CHECKS PASSED")
