from __future__ import annotations

import sympy as sp

from substrate_framework.euler_two_label_lock import (
    full_core_cancellation_first_integral,
    full_core_cancellation_residual,
    lane_emden_center_response_slope,
    lane_emden_full_core_response,
    material_contour_mean_tangent,
    material_radial_displacement_for_tag_tangent,
    weighted_contour_zero_mean,
)

s, t = sp.symbols("s t", positive=True)
p = sp.Integer(6)
P = 1 - s**2
Mp = sp.integrate(t * (1 - t**2) ** p, (t, 0, s))
Mp1 = sp.integrate(t * (1 - t**2) ** (p + 1), (t, 0, s))

# Independent Volterra reduction for unit physical scales.
f = P**p * sp.diff(P, s)
outer = sp.integrate(f.subs(s, t), (t, s, 1))
inner = sp.integrate(t**2 * (1 - t**2) ** p * sp.diff(1 - t**2, t), (t, 0, s))
v = sp.simplify(-sp.Rational(1, 2) * (inner / s + s * outer))
v_reduced = sp.simplify(Mp1 / ((p + 1) * s))
assert sp.simplify(v - v_reduced) == 0
print("PASS 1: exact power-law Volterra moment reduction")

phi_prime = -Mp / s
response_direct = sp.simplify(p * P ** (p - 1) * v + P**p * phi_prime)
response_api = lane_emden_full_core_response(P, p, 1, 1, 1, 1, 1, s, Mp, Mp1)
assert sp.simplify(response_direct - response_api) == 0
print("PASS 2: exact full-core response formula")

D = sp.simplify(p * Mp1 / (p + 1) - P * Mp)
assert sp.simplify(response_api - P ** (p - 1) * D / s) == 0

# Derive the physical Lane--Emden sign limits from generic positive data.
p_g = sp.symbols("p_g", positive=True)
P0 = sp.symbols("P0", positive=True)
Mp_edge, Mp1_edge = sp.symbols("Mp_edge Mp1_edge", positive=True)
D_center_coefficient = sp.simplify(
    p_g * P0 ** (p_g + 1) / (2 * (p_g + 1))
    - P0 * P0**p_g / 2
)
D_edge = sp.simplify(p_g * Mp1_edge / (p_g + 1) - 0 * Mp_edge)
assert sp.simplify(
    D_center_coefficient + P0 ** (p_g + 1) / (2 * (p_g + 1))
) == 0
assert D_center_coefficient.is_negative is True
assert D_edge.is_positive is True
print("PASS 3: generic Lane--Emden sign factor has negative center and positive edge")

slope = lane_emden_center_response_slope(1, p, 1, 1, 1, 1, 1)
assert slope == -sp.Rational(1, 14)
assert sp.simplify(sp.limit(response_api / s, s, 0) - slope) == 0
print("PASS 4: exact negative smooth-center coefficient")

assert response_api.subs(s, sp.Rational(1, 10)) < 0
assert response_api.subs(s, sp.Rational(9, 10)) > 0
print("PASS 5: nonphysical polynomial-profile algebra regression only")

Pfun = sp.Function("P")(s)
Gfun = sp.Function("G")
Ps = sp.diff(Pfun, s)
Pss = sp.diff(Pfun, s, 2)
Gp = sp.diff(Gfun(Pfun), Pfun)
Gpp = sp.diff(Gfun(Pfun), Pfun, 2)
classification_residual = full_core_cancellation_residual(Ps, Pss, Gp, Gpp, s)
classification_invariant = full_core_cancellation_first_integral(s, Ps, Gp)
L1 = lambda expr: sp.diff(expr, s, 2) + sp.diff(expr, s) / s - expr / s**2
assert sp.simplify(
    L1(Ps * Gfun(Pfun)) - Gfun(Pfun) * L1(Ps)
    - Ps * classification_residual
) == 0
assert sp.simplify(
    sp.diff(classification_invariant, s)
    - s * Ps**2 * classification_residual
) == 0
print("PASS 6: regular full-core cancellation has exact first integral")

k, C0 = sp.symbols("k C0", positive=True)
Gp_center = -C0 / (k**3 * s**4)
G_center = sp.integrate(Gp_center * (-k * s), s)
v_center = sp.simplify(-(-k * s) * G_center)
assert sp.limit(s * v_center, s, 0) != 0
print("PASS 7: nonzero first-integral constant produces singular dipole field")

a, amplitude, profile_symbol = sp.symbols("a amplitude profile_symbol", nonzero=True)
exponential = amplitude * sp.exp(a * profile_symbol)
assert sp.simplify(sp.diff(exponential, profile_symbol) - a * exponential) == 0
assert exponential.subs(profile_symbol, 0) != 0
print("PASS 8: regular kernel is exponential and cannot vanish at finite edge")

Fp, flux = sp.symbols("Fp flux")
assert material_contour_mean_tangent(Fp, flux) == -Fp * flux
assert material_contour_mean_tangent(Fp, 0) == 0
print("PASS 9: fixed-leaf weighted contour-flux obstruction")

T, chi_s = sp.symbols("T chi_s", nonzero=True)
xi_s = material_radial_displacement_for_tag_tangent(T, chi_s)
assert sp.simplify(-xi_s * chi_s - T) == 0
print("PASS 10: explicit raw zero-mean material-tag tangent")

u, v, w = sp.symbols("u v w")
finite_target = weighted_contour_zero_mean((u, v, w), (1, 2, 3))
assert sp.simplify(sum(weight * value for weight, value in zip((1, 2, 3), finite_target, strict=True))) == 0
print("PASS 11: finite-Cao weighted contour compatibility row")

print("PASS: 11/11 exact checks")
