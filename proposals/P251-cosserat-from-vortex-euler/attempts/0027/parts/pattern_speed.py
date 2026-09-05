"""0027 Part 2 (v3) -- m=2 pattern speed, improved Bessel series and the
displaced-base-pressure condition.

Sheet conditions (first order in eta):
  (1) kinematic interior: u_r^i = -i w eta
  (2) kinematic exterior: u_r^o = -i w eta
  (3) displaced pressure jump: p'_i - p'_o = eta (dp0_o/dr - dp0_i/dr)
      = eta rho a (Om_o^2 - Om_i^2)
  (4) material-vorticity sheet strength: v'_theta,i - v'_theta,o
      = 2 (Om_i - Om_o) eta
      (single-tube limit Om_o -> 0 reproduces the recorded [v_t] = 2 Om eta)

Bessel series to the orders that break the k->0 degeneracy:
  J_2(x) = x^2/8 - x^4/96;  dJ_2/dr = lam^2 a/4 - lam^4 a^3/24
  K_2(x) = 2/x^2 - 1/2;     dK_2/dr = -1/a - 4/(lam^4 a^3)
"""
import sympy as sp

I = sp.I
a, k, w, Omi, Omo, rho = sp.symbols("a k omega Omega_i Omega_o rho", positive=True)
Li, Lo = sp.symbols("Li Lo", positive=True)   # lam_i^2, lam_o^2 (both O(k^2))
m = sp.Integer(2)

wt_i = w - m * Omi
wt_o = w - m * Omo
D_i = wt_i**2 - 4 * Omi**2
D_o = wt_o**2 - 4 * Omo**2

J = Li * a**2 / 8 - Li**2 * a**4 / 96
Jp = Li * a / 4 - Li**2 * a**3 / 24
K = 2 / (Lo * a**2) - sp.Rational(1, 2)
Kp = -1 / a - 4 / (Lo**2 * a**3)

A, B, eta = sp.symbols("A B eta")
ur_i = -I * (wt_i * A * Jp + 2 * Omi * m * A * J / a) / (rho * D_i)
vt_i = (2 * Omi * A * Jp + m * wt_i * A * J / a) / (rho * D_i)
ur_o = -I * (wt_o * B * Kp + 2 * Omo * m * B * K / a) / (rho * D_o)
vt_o = (2 * Omo * B * Kp + m * wt_o * B * K / a) / (rho * D_o)

eqs = [ur_i + I * w * eta, ur_o + I * w * eta,
       A * J - B * K - eta * rho * a * (Omo**2 - Omi**2),
       vt_o - vt_i - 2 * (Omi - Omo) * eta]

rows = []
for expr in eqs:
    expr = sp.cancel(sp.expand(expr))
    rows.append([sp.simplify(expr.coeff(A)), sp.simplify(expr.coeff(B)),
                 sp.simplify(expr.coeff(eta))])
Mat = sp.Matrix(rows)

d012 = sp.factor(sp.cancel(Mat.extract([0, 1, 2], [0, 1, 2]).det()))
d023 = sp.factor(sp.cancel(Mat.extract([0, 2, 3], [0, 1, 2]).det()))
d123 = sp.factor(sp.cancel(Mat.extract([1, 2, 3], [0, 1, 2]).det()))
print("minor 012 (in Li, Lo):", d012)
print("minor 023 (in Li, Lo):", d023)
print("minor 123 (in Li, Lo):", d123)

subs = {Li: k**2 * (4 * Omi**2 - wt_i**2) / wt_i**2,
        Lo: k**2 * (4 * Omo**2 - wt_o**2) / wt_o**2}
for name, d in (("012", d012), ("023", d023), ("123", d123)):
    series = sp.expand(sp.numer(sp.together(d.subs(subs))))
    print(f"\nminor {name}: leading k^2 coefficient:")
    print("  ", sp.factor(series / (k**4)))  # strip k^4 = (k^2)^2 from Li*Lo products
