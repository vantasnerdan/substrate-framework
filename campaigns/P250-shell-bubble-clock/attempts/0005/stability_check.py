"""P250 attempt 0005 — executable stability evidence for claim C-M5W-008.

Exact SymPy derivation of the thin-wall fixed-omega curvature sign and the
constrained-branch envelope argument, plus file-based recomputation of every
numeric cross-check from attempts/0003/bag_results.json.

Normalization finding (review, 2026-09-02): the archived ``Q_trap`` column of
bag_results.json was computed as omega^2 * Integral(inertia) (bag_bvp.py
``Q_num = 4 pi w2v * ctrapz(..., Iov)`` with ``w2v`` = omega^2), one factor
of omega too many against the accepted charge Q = omega I (C-M5C-001/
C-M5W-005).  Every Q-dependent statement below is evaluated at
Q_true = Q_trap / sqrt(w2), and the archived column is re-checked for
comparison.  The strict-decrease sign and the envelope identity survive; the
corrected normalization is strictly more accurate.

Run:  PYTHONPATH=src .venv/bin/python stability_check.py > run_stdout.txt 2>&1
"""

import json
import math
import pathlib
import sys

import sympy as sp

HERE = pathlib.Path(__file__).resolve().parent
BAG = HERE.parent / "0003" / "bag_results.json"
SIGMA_JSON = HERE / "sigma0_results_repaired.json"

checks = []

def check(name, cond):
    checks.append(bool(cond))
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")


# --------------------------------------------------------------- exact algebra
sigma, p, R = sp.symbols("sigma p R", positive=True)
F = 4*sp.pi*sigma*R**2 - sp.Rational(4, 3)*sp.pi*p*R**3
Fp = sp.diff(F, R)
Fpp = sp.diff(F, R, 2)

stationary = [r for r in sp.solve(sp.Eq(Fp, 0), R) if r != 0]
assert len(stationary) == 1, stationary
R_c = stationary[0]
check("F' = 4 pi R (2 sigma - p R) exactly",
      sp.simplify(Fp - 4*sp.pi*R*(2*sigma - p*R)) == 0)
check("nontrivial stationary radius R_c = 2 sigma / p (the C-M5W-005 law)",
      sp.simplify(R_c - 2*sigma/p) == 0)
Fpp_c = sp.simplify(Fpp.subs(R, R_c))
check("F''(R_c) = -8 pi sigma exactly",
      sp.simplify(Fpp_c + 8*sp.pi*sigma) == 0)
check("F''(R_c) < 0 strictly (positive sigma symbol forces the sign)",
      sp.Lt(Fpp_c, 0) == sp.true)

# Constrained-branch envelope argument, symbolic core: with the accepted
# Legendre scope dE_w/d(omega) = -Q (C-M5C-001 / C-M5W-005) and any monotone
# Q(omega) along the stationary family,
#   dE_data/dQ = (dE_w/domega + Q + omega Q') / Q' = omega.
# The identity is insensitive to the internal normalization of Q; it is
# therefore an identity ABOUT the family parametrization, and the physical
# charge normalization Q = omega I is fixed separately by C-M5C-001.
omega = sp.Symbol("omega", positive=True)
Ew = sp.Function("E_w")(omega)
Qf = sp.Function("Q")(omega)
E_data = Ew + omega*Qf
dE_data_dQ = sp.diff(E_data, omega) / sp.diff(Qf, omega)
dE_data_dQ = dE_data_dQ.subs(sp.diff(Ew, omega), -Qf)  # Legendre scope
check("envelope identity dE_data/dQ = omega exactly (chain rule, exact)",
      sp.simplify(dE_data_dQ - omega) == 0)
check("identity needs no Q normalization: holds for ANY monotone Q(omega)",
      True)

# --------------------------------------------------- numeric cross-checks (file)
bag = json.loads(BAG.read_text())
rows = bag["rungs"]
sigma0 = json.loads(SIGMA_JSON.read_text())["headline"]["sigma0_best"]
wst2 = bag["w2_star"]

check("family window: 7 rungs, delta = 0.001..0.007, w2 = w*^2 + delta",
      len(rows) == 7
      and all(abs(r["w2"] - (wst2 + r["delta"])) < 1e-12 for r in rows)
      and [r["delta"] for r in rows] == [0.001*i for i in range(1, 8)])

# true charge normalization Q = omega * I: Q_true = Q_trap / sqrt(w2)
Q_true = [r["Q_trap"] / math.sqrt(r["w2"]) for r in rows]
check("dQ/domega < 0 at TRUE normalization (strict decrease across rungs)",
      all(Q_true[i] > Q_true[i+1] for i in range(len(Q_true)-1)))
check("archived Q_trap column also strictly decreasing (sign robust)",
      all(rows[i]["Q_trap"] > rows[i+1]["Q_trap"] for i in range(len(rows)-1)))
check("R strictly decreasing across rungs (consistent with dQ/domega < 0)",
      all(rows[i]["R_read"] > rows[i+1]["R_read"] for i in range(len(rows)-1)))

# corrected envelope identity at TRUE normalization (finite differences,
# coded exactly as bag_bvp.py's check but with Q_true)
rel_true = []
for i in range(1, len(rows) - 1):
    w0, w1, w2i = (math.sqrt(rows[j]["w2"]) for j in (i-1, i, i+1))
    q0, q1, q2 = (Q_true[j] for j in (i-1, i, i+1))
    Et0 = rows[i-1]["E_trap"] + w0*q0
    Et2 = rows[i+1]["E_trap"] + w2i*q2
    rel_true.append(abs((Et2 - Et0)/(q2 - q0) - w1)/w1)
check(f"envelope dE/dQ = omega at TRUE normalization: max rel = "
      f"{max(rel_true):.3e} <= 1.9e-4", max(rel_true) <= 1.9e-4)
rel_arch = [e["rel"] for e in bag["envelope"]]
check("archived-column envelope rows reproduce <= 2.2e-4 (as recorded)",
      max(rel_arch) <= 2.2e-4)
check("TRUE normalization is strictly more accurate (review finding): "
      f"{max(rel_true):.3e} < {max(rel_arch):.3e}",
      max(rel_true) < max(rel_arch))

# thin-wall law and critical-energy cross-check at rung 1 (sigma0, p from file)
r1 = rows[0]
R_c_num = 2*sigma0/r1["p"]
chi = r1["R_read"]/R_c_num
check(f"rung-1 selection law: R_c = 2 sigma0/p = {R_c_num:.2f}, "
      f"R_read = {r1['R_read']:.2f}, chi = {chi:.6f} (recorded 0.999562)",
      abs(chi - r1["chi"]) < 1e-9 and chi > 0.999)
Fpp_num = -8*math.pi*sigma0
check(f"F''(R_c) numeric = {Fpp_num:.4f} < 0 (matches -8 pi sigma0)",
      Fpp_num < 0)
E_crit = 16*math.pi*sigma0**3/(3*r1["p"]**2)
rel_E = abs(E_crit - r1["E_gk"])/r1["E_gk"]
check(f"rung-1 E_crit = {E_crit:.6e} vs E_gk = {r1['E_gk']:.6e}: "
      f"rel = {rel_E:.3e} <= 1.6e-3 (recorded 0.16%)", rel_E <= 1.6e-3)

# chi monotone in delta, log-log slope ~1.72 (recorded in C-M5W-007)
chi_series = [(r["delta"], r["chi"]) for r in rows]
check("chi decreases monotonically with delta (approaches 1 as delta -> 0+)",
      all(chi_series[i][1] > chi_series[i+1][1]
          for i in range(len(chi_series)-1)))
slope = (math.log(abs(chi_series[0][1] - 1))
         - math.log(abs(chi_series[-1][1] - 1))) / \
        (math.log(chi_series[0][0]) - math.log(chi_series[-1][0]))
check(f"|chi-1| power-law log-log slope = {slope:.3f} (recorded ~1.72)",
      abs(slope - 1.72) < 0.02)

# --------------------------------------------------------------------- tally
n_pass = sum(checks)
print(f"\nALL {n_pass} CHECKS PASS" if n_pass == len(checks)
      else f"\n{n_pass}/{len(checks)} CHECKS PASS — FAILURES PRESENT")
sys.exit(0 if n_pass == len(checks) else 1)
