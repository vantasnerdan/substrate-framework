#!/usr/bin/env python3
"""Exact algebra checks for the P253/0040 thin-ring/action bridge."""

import sympy as s

l = s.symbols("l", integer=True, positive=True)
L, S, Gamma, rho, R, amp = s.symbols(
    "L S Gamma rho R amp", positive=True, finite=True
)
g = Gamma / (4 * s.pi * R**2)

# Fischer--Schopohl (14), with the harmonic sum S_l kept symbolic.
A = g * (l**2 * (L - 2 * S + s.Rational(1, 2))
         + s.Rational(3, 2) * S)
B = g * ((l**2 - 1) * (L - 2 * S + s.Rational(1, 2))
         - s.Rational(3, 2) * (S - 1))
assert s.simplify(s.diff(A, L) - g * l**2) == 0
assert s.simplify(s.diff(B, L) - g * (l**2 - 1)) == 0
assert s.simplify(s.limit(A / L, L, s.oo) - g * l**2) == 0
assert s.simplify(s.limit(B / L, L, s.oo) - g * (l**2 - 1)) == 0

# The oscillator and its fixed-l leading frequency.
nu = s.sqrt(A * B)
leading_nu = g * l * s.sqrt(l**2 - 1) * L
assert s.simplify(s.limit(nu / leading_nu, L, s.oo) - 1) == 0

# Direct Hamilton equations for q=z and p=rho*Gamma*R*r.
z, r, p_can = s.symbols("z r p_can", real=True)
H_can = (rho * Gamma * R * A * z**2 / 2
         + B * p_can**2 / (2 * rho * Gamma * R))
assert s.simplify(s.diff(H_can, p_can).subs(
    p_can, rho * Gamma * R * r) - B * r) == 0
assert s.simplify(-s.diff(H_can, z) + rho * Gamma * R * A * z) == 0

# Angular impulse of a traveling l-wave.  Orthogonality gives
# integral_0^(2*pi) sin(l*phi)^2 dphi=pi for positive integer l.
phi, time = s.symbols("phi time", real=True)
phase = l * phi - nu * time
z_mode = amp * s.cos(phase)
r_mode = nu * amp * s.sin(phase) / B
integrand = -rho * Gamma * R * r_mode * s.diff(z_mode, phi)
J2 = s.simplify(s.integrate(integrand, (phi, 0, 2 * s.pi)))
expected_J2 = s.pi * rho * Gamma * R * l * nu * amp**2 / B
assert s.simplify(J2 - expected_J2) == 0

# The physical push-forward generator is X_J=-[R,omega].  With
# sigma=Omega_KKS(q_c,q_s) and [R,q_c]=l*q_s, this fixes the positive
# moment-map derivative and agrees with the directly integrated filament J_z.
sigma_symbol = s.symbols("sigma_symbol", nonzero=True, real=True)
kks_matrix = s.Matrix([[0, sigma_symbol], [-sigma_symbol, 0]])
xj_coords = s.Matrix([0, -amp * l])
amplitude_coords = s.Matrix([1, 0])
moment_derivative = s.simplify(
    (xj_coords.T * kks_matrix * amplitude_coords)[0])
assert s.simplify(moment_derivative - amp * l * sigma_symbol) == 0
sigma_filament = 2 * s.pi * rho * Gamma * R * nu / B
j2_filament = s.simplify(l * sigma_filament / 2)
assert s.simplify(j2_filament - expected_J2 / amp**2) == 0

# Rotating-frame sign and automatic transversality coefficient.
nu_symbol = s.symbols("nu_symbol", nonzero=True, real=True)
Omega0 = -nu_symbol / l
kernel_scalar = -s.I * nu_symbol - Omega0 * s.I * l
transversality_scalar = -s.I * l
assert s.simplify(kernel_scalar) == 0
assert transversality_scalar != 0

# Vorticity-side rotation moment map: curl(-|x|^2 e_z/2)=e_z cross x.
x, y, zz = s.symbols("x y zz", real=True)
potential = s.Matrix([0, 0, -(x**2 + y**2 + zz**2) / 2])
curl = s.Matrix([
    s.diff(potential[2], y) - s.diff(potential[1], zz),
    s.diff(potential[0], zz) - s.diff(potential[2], x),
    s.diff(potential[1], x) - s.diff(potential[0], y),
])
assert curl == s.Matrix([-y, x, 0])

# Cao's Lane--Emden core gives the exact Gallay--Smets profile test (25)--(26).
srad = s.symbols("srad", positive=True)
U = s.Function("U")(srad)
W_p = -s.diff(U, srad, 2) - s.diff(U, srad) / srad
F_p = -srad * s.diff(U, srad)
D_p = 2 * F_p - srad**2 * W_p
assert s.simplify(s.diff(F_p, srad) - srad * W_p) == 0
assert s.simplify(D_p - (srad**2 * s.diff(U, srad, 2)
                         - srad * s.diff(U, srad))) == 0
J_p = s.factor(2 * F_p * W_p * srad**4 / D_p**2)

# Dimensionless Lane--Emden reduction of the Richardson monotonicity test.
yvar, avar, tvar, power = s.symbols(
    "yvar avar tvar power", positive=True, finite=True
)
dy = (yvar**2 + avar) / srad
da = avar * (2 - (power - 1) * yvar) / srad
dt_from_ratio = s.simplify((da * yvar - avar * dy) / yvar**2)
assert s.simplify(dt_from_ratio.subs(avar, tvar * yvar)
                  -tvar * (2 - tvar - power * yvar) / srad) == 0
J_dimensionless = 2 * srad**2 * yvar * avar / (2 * yvar - avar)**2
dlogJ = s.diff(s.log(J_dimensionless), srad)
dlogJ += s.diff(s.log(J_dimensionless), yvar) * dy
dlogJ += s.diff(s.log(J_dimensionless), avar) * da
expected_dlogJ = ((8 - tvar**2 - 2 * tvar
                   - power * yvar * (tvar + 2))
                  / (srad * (2 - tvar)))
assert s.simplify(dlogJ.subs(avar, tvar * yvar)-expected_dlogJ) == 0

# The H2 comparison closes without solving U explicitly.
h_comp = power * yvar - 2 * (2 - tvar)
radial_h_derivative = power * yvar * (yvar - tvar) + 2 * tvar * (2 - tvar)
at_first_zero = s.simplify(radial_h_derivative.subs(
    yvar, 2 * (2 - tvar) / power))
assert s.simplify(at_first_zero - 4 * (2 - tvar)**2 / power) == 0

# Derive, rather than enter, the center jets from the radial Lane--Emden
# equation.  With x=s^2 and lambda=U(0)^(p-1), the normalized equation is
# 4(x f_xx+f_x)+lambda*f^p=0.  Its recurrence through x^2 fixes f through
# x^3, which is the order needed to obtain t through s^4.
xrad, lam = s.symbols("xrad lam", positive=True, finite=True)
f_jet = (1 - lam * xrad / 4 + power * lam**2 * xrad**2 / 64
         - power * (3 * power - 2) * lam**3 * xrad**3 / 2304)
le_residual = s.series(
    4 * (xrad * s.diff(f_jet, xrad, 2) + s.diff(f_jet, xrad))
    + lam * f_jet**power, xrad, 0, 3
).removeO().expand()
assert s.simplify(le_residual) == 0
y_for_ratio = s.series(-2 * xrad * s.diff(f_jet, xrad) / f_jet,
                       xrad, 0, 4).removeO().expand()
y_jet = s.series(y_for_ratio, xrad, 0, 3).removeO().expand()
a_jet = s.series(lam * xrad * f_jet**(power - 1),
                 xrad, 0, 4).removeO().expand()
t_jet = s.series(a_jet / y_for_ratio, xrad, 0, 3).removeO().expand()
expected_y_jet = (lam * xrad / 2
                  + (2 - power) * lam**2 * xrad**2 / 16)
expected_t_jet = (2 - power * lam * xrad / 4
                  + power * (3 * power - 4) * lam**2 * xrad**2 / 96)
assert s.simplify(y_jet - expected_y_jet) == 0
assert s.simplify(t_jet - expected_t_jet) == 0
h_jet = s.series(power * y_jet - 2 * (2 - t_jet),
                 xrad, 0, 3).removeO().expand()
center_h2 = s.simplify(h_jet.coeff(xrad, 2) / lam**2)
assert s.simplify(center_h2 - power / 24) == 0
h2_boundary = s.simplify((
    (2 - tvar) * (tvar + 4) - power * yvar * (tvar + 2)
).subs(yvar, 2 * (2 - tvar) / power))
assert s.simplify(h2_boundary + tvar * (2 - tvar)) == 0

# Exact exterior m=1 Dirichlet-to-Neumann threshold coefficient.  The
# displayed jet solves the modified-Bessel equation through the retained
# orders and yields x*K_1'(x)/K_1(x)=-1+x^2(log(x/2)+gamma_E)+O(x^4 log^2 x).
xext = s.symbols("xext", positive=True, finite=True)
k1_jet = (1 / xext + xext / 2
          * (s.log(xext / 2) + s.EulerGamma - s.Rational(1, 2)))
bessel_residual = s.expand(
    xext**2 * s.diff(k1_jet, xext, 2) + xext * s.diff(k1_jet, xext)
    - (xext**2 + 1) * k1_jet)
assert s.limit(bessel_residual / (xext**3 * s.log(xext)),
               xext, 0, dir="+") == -s.Rational(1, 2)
dtn_ratio = s.simplify(xext * s.diff(k1_jet, xext) / k1_jet)
dtn_log_remainder_limit = s.limit(
    (dtn_ratio + 1) / xext**2
    - (s.log(xext / 2) + s.EulerGamma), xext, 0, dir="+")
assert dtn_log_remainder_limit == 0

print("A_l =", A)
print("B_l =", B)
print("fixed_l_leading_nu =", leading_nu)
print("traveling_mode_J_z_quadratic =", expected_J2)
print("physical_orbit_moment_derivative =", moment_derivative)
print("filament_KKS_sigma =", sigma_filament)
print("critical_pattern_speed =", Omega0)
print("CR_transversality_multiplier =", transversality_scalar)
print("rotation_moment_map_curl =", curl.T)
print("Cao_profile_Richardson_J =", J_p)
print("Cao_H2_log_derivative =", expected_dlogJ)
print("H2_comparison_center_coefficient =", center_h2)
print("H2_comparison_zero_derivative =", at_first_zero)
print("H2_log_numerator_at_comparison_boundary =", h2_boundary)
print("exterior_m1_DtN_log_remainder_limit =", dtn_log_remainder_limit)
print("ALL 25 EXACT KELVIN-BRIDGE CHECKS PASSED")
