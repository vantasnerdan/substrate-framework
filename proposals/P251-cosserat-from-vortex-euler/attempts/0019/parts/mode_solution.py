"""0019 Part 1 -- normal-mode solution about solid rotation, built in parts.
Goal: arbitrate the velocity field forms and reduce incompressibility to the
pressure equation. Everything checked by residual, nothing by inspection.

Conventions:
  mode ~ exp(i(m theta + k z - w t)),  Doppler wt = w - m*Om
  background: solid rotation Omega (Rankine core, r < a)
"""
import sympy as sp

I = sp.I
r, Om, rho = sp.symbols("r Omega rho", positive=True)
m_, k, w = sp.symbols("m k omega")
wt = w - m_ * Om
p = sp.Function("p")(r)
pp = p                      # pressure eigenfunction
dpp = sp.diff(p, r)

# Linearized Euler about solid rotation (lab frame), pattern ansatz:
#   r: -i wt vr = -dpp/rho - 2 Om vt          (Coriolis r-comp: 2Om(zhat x v)_r = -2Om vt)
#   t: -i wt vt = -i m pp/(r rho) + 2 Om vr   (Coriolis t-comp: 2Om vr)
#   z: -i wt vz = -i k pp/rho
# PART 1a: solve the momentum system for (vr, vt, vz) in terms of (pp, dpp).
vr, vt, vz = sp.symbols("vr vt vz")
sol = sp.solve(
    [
        sp.Eq(-I * wt * vr, -dpp / rho - 2 * Om * vt),
        sp.Eq(-I * wt * vt, -I * m_ * pp / (r * rho) + 2 * Om * vr),
        sp.Eq(-I * wt * vz, -I * k * pp / rho),
    ],
    [vr, vt, vz],
    dict=True,
)[0]
vr_mom = sol[vr]
vt_mom = sol[vt]
vz_mom = sol[vz]
Dtil = sp.expand(wt**2 - 4 * Om**2)                   # epicyclic denominator (to verify)
print("Dtil form matches solve denominator:",
      sp.simplify(sp.denom(sp.together(vr_mom * rho)) - r * Dtil) == 0)

# PART 1b: candidate compact forms (epicyclic v+- structure):
#   v_r = -i (wt p' + 2 Om m p/r) / (rho Dtil)
#   v_t =    (2 Om p' + m wt p/r) / (rho Dtil)        (real, no i)
#   v_z =    k p / (rho wt)
vr_cand = sp.simplify(-I * (wt * dpp + 2 * Om * m_ * pp / r) / (rho * Dtil))
vt_plus = sp.simplify((2 * Om * dpp + wt * m_ * pp / r) / (rho * Dtil))
vt_minus = sp.simplify(-vt_plus)
vz_plus = k * pp / (rho * wt)
vz_minus = -vz_plus

# PART 1c: arbitrate by direct difference against the simplified solve output
vr_s = sp.simplify(vr_mom)
vt_s = sp.simplify(vt_mom)
vz_s = sp.simplify(vz_mom)
print("vr compact == solved:", sp.simplify(vr_cand - vr_s) == 0)
print("vt_plus  == solved :", sp.simplify(vt_plus - vt_s) == 0)
print("vt_minus == solved :", sp.simplify(vt_minus - vt_s) == 0)
print("vz_plus  == solved :", sp.simplify(vz_plus - vz_s) == 0)
print("vz_minus == solved :", sp.simplify(vz_minus - vz_s) == 0)

# PART 1d: momentum residuals of the winning compact triple (each must be 0)
vt_cand = vt_plus if sp.simplify(vt_plus - vt_s) == 0 else vt_minus
vz_cand = vz_plus if sp.simplify(vz_plus - vz_s) == 0 else vz_minus
print("1d r-eq residual:", sp.simplify(-I * wt * vr_cand - (-dpp / rho - 2 * Om * vt_cand)))
print("1d t-eq residual:", sp.simplify(-I * wt * vt_cand - (-I * m_ * pp / (r * rho) + 2 * Om * vr_cand)))
print("1d z-eq residual:", sp.simplify(-I * wt * vz_cand - (-I * k * pp / rho)))

# PART 1e: incompressibility -> pressure equation (Poincare reduction).
# Expect: p'' + p'/r - m^2 p/r^2 + lam^2 p = 0, lam^2 = k^2 (4 Om^2 - wt^2)/wt^2.
div = sp.diff(r * vr_cand, r) / r + I * m_ * vt_cand / r + I * k * vz_cand
div = sp.simplify(sp.together(div))
num, den = sp.fraction(sp.cancel(div * rho * Dtil * wt))
num = sp.expand(num)
c2 = sp.expand(num).coeff(sp.diff(p, r, 2))
c1 = sp.expand(num).coeff(sp.diff(p, r))
c0 = sp.expand(num).coeff(p)
print("1e p'' coeff / (-I r^2 wt^2):", sp.simplify(c2 / (-I * r**2 * wt**2)))
print("1e p'  coeff / (-I r wt^2) :", sp.simplify(c1 / (-I * r * wt**2)))
lam2_from = sp.simplify(sp.cancel(c0 / c2) + m_**2 / r**2)
lam2_target = k**2 * (4 * Om**2 - wt**2) / wt**2
print("1e lam^2 - Poincare target :", sp.simplify(lam2_from - lam2_target))
