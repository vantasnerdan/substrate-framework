"""0019 addendum -- arbitrate verify_cst002's recorded v-forms vs part-1 forms.

Part-1 proven (momentum residuals zero, sympy):
    v_r = -i (wt p' + 2 Om m p/r) / (rho (wt^2 - 4 Om^2))
    v_t =    (2 Om p' + m wt p/r) / (rho (wt^2 - 4 Om^2))
verify_cst002.py records (its check_poincare_reduction):
    v_r = +i (wt p' - 2 Om m p/r) / (rho (4 Om^2 - wt^2))
    v_t =    (2 Om p' - m wt p/r) / (rho (4 Om^2 - wt^2))
Test both against the linearized-momentum residuals.
"""
import sympy as sp

I = sp.I
r, Om, rho = sp.symbols("r Omega rho", positive=True)
m, k, w = sp.symbols("m k omega")
wt = w - m * Om
p = sp.Function("p")(r)
dpp = sp.diff(p, r)

# their forms:
D_theirs = 4 * Om**2 - wt**2
vr_t = I * (wt * dpp - 2 * Om * m * p / r) / (rho * D_theirs)
vt_t = (2 * Om * dpp - m * wt * p / r) / (rho * D_theirs)
vz_t = k * p / (rho * wt)

# part-1 forms:
D_mine = wt**2 - 4 * Om**2
vr_m = -I * (wt * dpp + 2 * Om * m * p / r) / (rho * D_mine)
vt_m = (2 * Om * dpp + m * wt * p / r) / (rho * D_mine)

def residuals(vr_x, vt_x, vz_x):
    r1 = sp.simplify(-I * wt * vr_x - (-dpp / rho - 2 * Om * vt_x))
    r2 = sp.simplify(-I * wt * vt_x - (-I * m * p / (r * rho) + 2 * Om * vr_x))
    r3 = sp.simplify(-I * wt * vz_x - (-I * k * p / rho))
    return r1, r2, r3

print("verify_cst002 forms residuals (r, theta, z):")
for nm, rr in zip(("r", "t", "z"), residuals(vr_t, vt_t, vz_t)):
    print(f"  {nm}: {sp.simplify(rr) if rr != 0 else 0}")
print("part-1 forms residuals:")
for nm, rr in zip(("r", "t", "z"), residuals(vr_m, vt_m, vz_t)):
    print(f"  {nm}: {sp.simplify(rr) if rr != 0 else 0}")
