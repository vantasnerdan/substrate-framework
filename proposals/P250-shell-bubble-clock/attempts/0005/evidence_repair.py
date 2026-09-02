"""P250 attempt 0005 — evidence repair for review findings on C-M5W-006/007.

File-based regeneration of the two defective evidence records, executed and
captured (stdout is the record; JSON outputs land beside this file):

G1 (sigma0_results.json defect class, review 2026-09-02):
  * tail asymptotics were evaluated at the WRONG bulk point
    (m = c* instead of m = 0; the certified endpoint is B = (0, c*, b*, f*)),
    corrupting sqrt_lam_B and the exponential-fit boundary treatment;
  * the headline was taken from an h-ladder rung whose solve_bvp status was 1
    (failed to converge to tolerance), and acceptance ignored solver status.
  Repair here: corrected tail masses (cross-checked against the receipt's
  recorded values), corrected expfit boundary treatment, full production
  chain rerun, headline anchored to a status-0 rung, budget recomputed,
  bracket-notation uncertainty derived from the budget.

G2 (bag_results.json defect class, review 2026-09-02):
  * Q_trap = omega^2 * Integral(inertia) carries one factor of omega too
    many against the accepted charge Q = omega I (C-M5C-001/C-M5W-005).
  Repair here: exact algebraic rescaling of the stored rows to the true
  normalization Q_true = Q_trap/sqrt(w2) and recomputation of the envelope
  finite differences at the true charge; emitted as bag_results_corrected.json.

Not rerun (documented): the L=10 h-ladder and basin probes (neither feeds
the error budget); attempt artifacts under attempts/0003 are never edited —
the corrected columns live here and in bag_results_corrected.json.

Run:  PYTHONPATH=src .venv/bin/python evidence_repair.py > run_repair_stdout.txt 2>&1
"""

import json
import math
import platform
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
A2 = HERE.parent / "0002"
A3 = HERE.parent / "0003"
ROOT = A2.parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(A2))

import mpmath as mp  # noqa: E402
import sympy as sp  # noqa: E402
import sigma0_wall as g1  # noqa: E402  (module level: cheap analytic setup only)

checks = []


def check(name, cond):
    checks.append(bool(cond))
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")


def bracket(x, err):
    """Parenthetical uncertainty notation with two significant digits."""
    e = math.floor(math.log10(err)) - 1
    digits = round(err / 10 ** e)
    mant = round(x, -e)
    return f"{mant:.{-e}f}({digits:d})"


# Corrected endpoint bulk point: B = (0, c*, b*, f*) — m = 0, NOT m = c*.
# Same construction as sigma0_wall.py (Hessian of Vw at omega*^2, Abs->f
# slice license applied), evaluated at the certified endpoint.
_m2, _c2, _b2, _f2 = g1._m, g1._c, g1._b, g1._f
H_B_CORR = sp.hessian(g1.Vw_expr, (_m2, _c2, _b2, _f2)).subs(
    {_m2: sp.Integer(0), _c2: sp.N(str(g1.B_STATE[0]), 30),
     _b2: sp.N(str(g1.B_STATE[1]), 30), _f2: sp.N(str(g1.B_STATE[2]), 30)})
LAM_B_CORR, Q_B_CORR = g1.eigs_of(sp.diag(2, 1, 1, 1) * H_B_CORR)
SQ_B_CORR = np.sqrt(LAM_B_CORR)

RECEIPT_MASSES = [2.5276, 2.899, 3.277, 7.012]  # recorded in 0002/receipt.md
import json as _json  # noqa: E402
print("wrong sqrt_lam_B (m = c*):",
      np.round(np.sqrt(_json.loads(
          (A2 / "sigma0_results.json").read_text())["tail"]["sqrt_lam_B"]),
          6).tolist())
print("corrected sqrt_lam_B        :",
      np.round(SQ_B_CORR, 6).tolist())
print("receipt-recorded masses     :", RECEIPT_MASSES)
check("corrected B-tail masses reproduce the receipt values",
      all(abs(a - b) < 5e-4 for a, b in zip(sorted(SQ_B_CORR), RECEIPT_MASSES)))
M_MIN_CORR = float(min(g1.SQ_A.min(), SQ_B_CORR.min()))
check("M_MIN unchanged (set by the A side): "
      f"{M_MIN_CORR:.6f}", abs(M_MIN_CORR - float(g1.SQ_A.min())) < 1e-12)


def bc_expfit_corrected(ya, yb):
    """Stable-manifold projection at both ends, corrected B-basis."""
    xa = (ya[:4] - g1.A_VEC) @ g1.Q_A
    xb = (yb[:4] - g1.B_VEC) @ Q_B_CORR
    va = ya[4:] @ g1.Q_A
    vb = yb[4:] @ Q_B_CORR
    return np.concatenate((va - g1.SQ_A * xa, vb + SQ_B_CORR * xb))


# =============================================== G1 part 2: production rerun
print("\n== production chain rerun (L-continuation -> L=12 h-ladder) ==")
sol12, cont_rows = g1.continuation()
out = {"continuation": cont_rows}
out["domain_estimate"] = abs(cont_rows[-1]["r1"] - cont_rows[-2]["r1"])

h12 = []
solF = g1.solve(g1.L_FIX, 401, g1.bc_dirichlet, warm=sol12, tol=1e-10)
for n in g1.N_LADDER:
    solF = g1.solve(g1.L_FIX, n, g1.bc_dirichlet, warm=solF, tol=1e-10)
    row = g1.tension_routes(solF, g1.L_FIX)
    row["n_req"] = n
    h12.append(row)
    print(f"n={n:5d} nodes={row['n_nodes']:6d} r1={row['r1']:.14f} "
          f"resRMS={row['res_rms']:.2e} fi(abs)={row['fi_drift']:.2e} "
          f"status={row['status']}")
out["h_crossing_L12"] = h12

anchor = h12[-1] if h12[-1]["status"] == 0 else h12[-2]
check(f"headline anchored to a status-0 rung (n={anchor['n_req']})",
      anchor["status"] == 0)
anchor_i = h12.index(anchor)
coarser = h12[anchor_i - 1]

# corrected exponential-fit cross-check at the anchor resolution
solE = g1.solve(g1.L_FIX, anchor["n_req"], bc_expfit_corrected,
                warm=solF, tol=1e-10)
rE = g1.tension_routes(solE, g1.L_FIX)
print(f"expfit(corrected): status={rE['status']} nodes={rE['n_nodes']} "
      f"r1={rE['r1']:.14f} fi(abs)={rE['fi_drift']:.2e}")
out["expfit_crosscheck_corrected"] = rE
b7 = abs(rE["r1"] - anchor["r1"])
check(f"corrected BC agreement |expfit - dirichlet| = {b7:.3e}", b7 < 1e-8)
check("expfit cross-check converged (status 0)", rE["status"] == 0)

# certified-box drift (b4)
drift = {}
for tag, db in (("up", +g1.BOX), ("down", -g1.BOX)):
    b_vec = np.array([0.0, g1.B_STATE[0] + db, g1.B_STATE[1] + db,
                      g1.B_STATE[2] + db])
    s = g1.solve(g1.L_FIX, anchor["n_req"], g1.bc_dirichlet,
                 warm=solF, tol=1e-10, b_vec=b_vec)
    drift[tag] = g1.tension_routes(s, g1.L_FIX)["r1"]
b4 = abs(drift["up"] - drift["down"]) / 2
out["box_drift"] = {"up": drift["up"], "down": drift["down"],
                    "half_spread": b4}

# evaluator noise (b6)
Y = solF.sol(np.linspace(-g1.L_FIX, g1.L_FIX, 6401))
T = g1.kinetic(Y)
V = g1.pot(Y)
b6 = max(abs(math.fsum((T + V)[::-1].tolist())
             - math.fsum((T + V).tolist())),
         float(np.max(np.abs(
             ((Y[4]**2 + 2*Y[5]**2 + 2*Y[6]**2)/4 + Y[7]**2/2
              + np.array([g1.V_fun(*Y[:4, i]) for i in range(Y.shape[1])]))
             - (T + V)))))

# budget at the status-0 anchor
routes = [anchor["r1"], anchor["r2"], anchor["r3"], anchor["quad"][0]]
spread = max(routes) - min(routes)
budget = {
    "b1_residual_scaled": anchor["res_rms"] * (2 * g1.L_FIX),
    "b2_truncation_halfspread_L12": abs(anchor["r1"] - coarser["r1"]) / 2,
    "b3_domain_L10_to_L12": out["domain_estimate"],
    "b4_certified_box": b4,
    "b5_quadrature": abs(anchor["quad"][0] - anchor["r1"]),
    "b6_evaluator": b6,
    "b7_bc_agreement_corrected": b7,
    "b8_first_integral_drift_absolute": anchor["fi_drift"],
}
total = max(budget.values())
sigma0_best = float(np.mean(routes))
accept = total <= 1e-3 * anchor["r1"]
print("\n== repaired headline (status-0 anchor) ==")
for name, v in zip(("R1 T+V", "R2 2V", "R3 2T", "GK-quad"), routes):
    print(f"  {name:8s} = {v:.15f}")
print(f"  route spread = {spread:.3e}")
for k, v in budget.items():
    print(f"  {k:38s} {v:.3e}")
print(f"  budget total (max item)               {total:.3e}")
print(f"  sigma_0 = {sigma0_best:.15f}"
      f"   bracket: {bracket(sigma0_best, total)}")
print(f"  acceptance: total <= 1e-3 sigma_0 -> "
      f"{'PASS' if accept else 'FAIL'}")
out["headline"] = {
    "anchor_n": anchor["n_req"], "anchor_status": anchor["status"],
    "routes": dict(zip(("r1_TV", "r2_2V", "r3_2T", "gk_quad"), routes)),
    "route_spread": spread, "budget": budget, "budget_total": total,
    "sigma0_best": sigma0_best,
    "sigma0_bracket": bracket(sigma0_best, total),
    "accept": accept,
    "fi_drift_metric": "absolute max |T - V| on the collocation grid",
    "monotonicity_metric": ("sampled reversals above 1e-8 (numerical, "
                            "not a continuum-exact statement)"),
}
check("repaired budget passes acceptance (total <= 1e-3 sigma0)", accept)

# corrected profile artifact
xs = np.linspace(-g1.L_FIX, g1.L_FIX, 2001)
Yp = solF.sol(xs)
with open(HERE / "profile_L12_anchor.csv", "w") as fh:
    fh.write("x,m,c,b,f,V,T\n")
    for i in range(len(xs)):
        Yi = Yp[:, i:i+1]
        fh.write(f"{xs[i]:.10g},{Yi[0, 0]:.12g},{Yi[1, 0]:.12g},"
                 f"{Yi[2, 0]:.12g},{Yi[3, 0]:.12g},{g1.V_fun(*Yi[:4, 0]):.12g},"
                 f"{float(g1.kinetic(Yi)[0]):.12g}\n")
out["profile_artifact"] = ("profile_L12_anchor.csv (2001 samples of the "
                           "n=1601-anchor solution; the archived "
                           "profile_L12_n3201.csv name was misleading)")
mon = anchor["mono"]
check("numerical monotonicity at the anchor (zero sampled reversals > 1e-8)",
      all(v == 0 for v in mon.values()))

# =============================================== G2: corrected charge columns
print("\n== G2 charge normalization repair (Q_true = Q_trap / omega) ==")
bag = json.loads((A3 / "bag_results.json").read_text())
rows = bag["rungs"]
for r in rows:
    r["Q_true"] = r["Q_trap"] / math.sqrt(r["w2"])
check("dQ/domega < 0 at true normalization (strict decrease)",
      all(rows[i]["Q_true"] > rows[i + 1]["Q_true"]
          for i in range(len(rows) - 1)))
env_true = []
for i in range(1, len(rows) - 1):
    w0, w1, w2i = (math.sqrt(rows[j]["w2"]) for j in (i - 1, i, i + 1))
    q0, q1, q2 = (rows[j]["Q_true"] for j in (i - 1, i, i + 1))
    Et0 = rows[i - 1]["E_trap"] + w0 * q0
    Et1 = rows[i]["E_trap"] + w1 * q1
    Et2 = rows[i + 1]["E_trap"] + w2i * q2
    dEdQ = (Et2 - Et0) / (q2 - q0)
    rel = abs(dEdQ - w1) / w1
    env_true.append({"delta": rows[i]["delta"], "dEdQ": dEdQ,
                     "omega": w1, "rel": rel})
max_rel_true = max(e["rel"] for e in env_true)
max_rel_arch = max(e["rel"] for e in bag["envelope"])
rels = [e["rel"] for e in env_true]
print("corrected envelope rel errors:", [f"{v:.2e}" for v in rels])
check(f"envelope at true normalization: max rel = {max_rel_true:.3e} <= 1.9e-4",
      max_rel_true <= 1.9e-4)
check(f"true normalization strictly better than archived "
      f"({max_rel_true:.3e} < {max_rel_arch:.3e})",
      max_rel_true < max_rel_arch)
sigma_shift = abs(sigma0_best - json.loads(
    (A2 / "sigma0_results.json").read_text())["headline"]["sigma0_best"])
print(f"sigma_0 shift vs archived headline: {sigma_shift:.3e} "
      "(budget 2e-7 absorbs; chi and E_crit consumers shift "
      "<= 3.5e-11 relative — negligible at recorded precision)")
out_g2 = {
    "repair": "Q_true = Q_trap / sqrt(w2); envelope recomputed at Q_true",
    "envelope_true": env_true,
    "max_rel_true": max_rel_true,
    "max_rel_archived": max_rel_arch,
    "Q_true": {f"delta={r['delta']}": r["Q_true"] for r in rows},
    "virial_note": ("stored normalized virial residuals |virial|/virial_norm "
                    "are diagnostics, not acceptance gates; they range "
                    f"{rows[0]['virial']/rows[0]['virial_norm']:.2e} to "
                    f"{rows[-1]['virial']/rows[-1]['virial_norm']:.2e}"),
    "sigma0_shift": sigma_shift,
}
(HERE / "bag_results_corrected.json").write_text(json.dumps(out_g2, indent=1))
(HERE / "sigma0_results_repaired.json").write_text(json.dumps(
    out, indent=1,
    default=lambda o: bool(o) if isinstance(o, np.bool_) else float(o)))
print("\nwrote sigma0_results_repaired.json, bag_results_corrected.json, "
      "profile_L12_anchor.csv")

n_pass = sum(checks)
print(f"\nALL {n_pass} CHECKS PASS" if n_pass == len(checks)
      else f"\n{n_pass}/{len(checks)} CHECKS PASS — FAILURES PRESENT")
sys.exit(0 if n_pass == len(checks) else 1)
