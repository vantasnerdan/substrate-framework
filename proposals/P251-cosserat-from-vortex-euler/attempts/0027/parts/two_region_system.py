"""0027 Part 1 -- two-region Rankine mode problem (core vorticity 2 Om_i,
ambient solid rotation Om_o), momentum-residual-correct forms (0019).

Background: solid rotation Om_i for r < a, solid rotation Om_o for r > a;
vortex sheet at r = a carrying base jump [v_t]0 = (Om_o - Om_i) a.
Modes exp(i(m theta + k z - w t)); each region has its own Doppler
wt_i = w - m Om_i, wt_o = w - m Om_o, and epicyclic denominator
D = wt^2 - 4 Om^2 (the 0019-corrected Coriolis structure).

Sheet conditions at r = a (first order in eta):
  (1) kinematic (both sides): u_r^i = u_r^o = -i w eta
  (2) pressure continuity:    p_i(a) = p_o(a)
  (3) frozen sheet strength:  [v_t'] = 2 (Om_i - Om_o) eta
      (single-tube limit Om_o -> 0 reproduces the recorded [v_t] = 2 Om eta)

Outputs: the exact linear system; determinant; k->0 pattern-speed expansion
for m = 2 via Bessel series; single-tube limit check against
omega = Om_i[(m-1) - (ka)^2/6 + O((ka)^4)].
"""
import sympy as sp

I = sp.I
a, k, w, Omi, Omo, rho = sp.symbols("a k omega Omega_i Omega_o rho", positive=True)
m = sp.Symbol("m", positive=True, integer=True)

for region, Om in (("i", Omi), ("o", Omo)):
    wt = w - m * Om
    D = wt**2 - 4 * Om**2
    lam2 = k**2 * (4 * Om**2 - wt**2) / wt**2
    print(f"region {region}: wt = w - m*{region}, D = {sp.factor(D)}, "
          f"lam^2 = {sp.simplify(lam2)}")

wt_i = w - m * Omi
wt_o = w - m * Omo
D_i = wt_i**2 - 4 * Omi**2
D_o = wt_o**2 - 4 * Omo**2
lam_i = k * sp.sqrt(4 * Omi**2 - wt_i**2) / wt_i   # symbolic sqrt; J species
lam_o = k * sp.sqrt(4 * Omo**2 - wt_o**2) / wt_o

# Bessel symbols and radial derivatives (prime = d/dr)
J = sp.Symbol("J")     # J_m(lam_i a)
Jp = sp.Symbol("Jp")   # dJ_m/d(r) at a
K = sp.Symbol("K")     # K_m(lam_o a)
Kp = sp.Symbol("Kp")
A, B, eta = sp.symbols("A B eta")

# interior: p = A J_m(lam_i r);  exterior: p = B K_m(lam_o r)
ur_i = -I * (wt_i * A * Jp + 2 * Omi * m * A * J / a) / (rho * D_i)
vt_i = (2 * Omi * A * Jp + m * wt_i * A * J / a) / (rho * D_i)
ur_o = -I * (wt_o * B * Kp + 2 * Omo * m * B * K / a) / (rho * D_o)
vt_o = (2 * Omo * B * Kp + m * wt_o * B * K / a) / (rho * D_o)

eqs = [
    sp.Eq(ur_i, -I * w * eta),                          # kinematic interior
    sp.Eq(ur_o, -I * w * eta),                          # kinematic exterior
    sp.Eq(A * J, B * K),                                # pressure continuity
    sp.Eq(vt_o - vt_i, 2 * (Omi - Omo) * eta),          # frozen sheet strength
]
# simpler: build the 4x3 matrix (4 equations, 3 unknowns A, B, eta)
rows = []
for eq in eqs:
    expr = sp.expand(eq.lhs - eq.rhs)
    rows.append([expr.coeff(A), expr.coeff(B), expr.coeff(eta)])
Mat = sp.Matrix(rows).applyfunc(sp.simplify)
print("\nlinear system coefficients (rows = conditions, cols = A, B, eta):")
sp.pprint(Mat)

# compatibility: 3x4 -> dispersion from any 3x3 minor with the kinematic pair
# eliminate eta from (1) and (2): ur_i = ur_o relation independent of eta
eta_elim = sp.simplify((Mat[0, 0] - Mat[1, 0]) * 0 + (Mat[0, 2] - Mat[1, 2]))
print("\neta-coefficient equality condition (ur_i = ur_o => both = -i w eta):",
      sp.simplify(Mat[0, 2] - Mat[1, 2]) != 0)

# dispersion: eliminate A, B, eta -> det of the 3x3 minors; use rows (0,1,3) and (0,2,3)
det1 = sp.simplify(Mat.extract([0, 1, 3], [0, 1, 2]).det())
det2 = sp.simplify(Mat.extract([0, 2, 3], [0, 1, 2]).det())
print("\ndet(rows 0,1,3):", det1)
print("det(rows 0,2,3):", det2)
print("\n(part 1 structural output; dispersion series in part 2)")
