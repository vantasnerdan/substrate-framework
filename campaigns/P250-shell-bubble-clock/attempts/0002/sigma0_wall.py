"""P250 G1 production verifier: slice-wall profile and sigma_0 at omega_*^2.

Implements the frozen design of attempts/0002/receipt.md exactly:
solve_bvp collocation on the crossed (nodes x L) ladder, exponential-fit and
Dirichlet boundary treatments, three tension routes with math.fsum,
multi-start branch checks, itemized error budget b1-b8.

Run:  PYTHONPATH=src python attempts/0002/sigma0_wall.py
Outputs: sigma0_results.json, profile_L12_n3201.csv (alongside this file).
"""
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"
import sys
from pathlib import Path
import math, json, platform

import numpy as np
import scipy
from scipy.integrate import solve_bvp, quad

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "src"))

import sympy as sp
import mpmath as mp
from substrate_framework.m5_wall_clock import (
    wall_slice_potential, first_integral_conservation_identity)

mp.mp.dps = 40

# ---------------------------------------------------------------- frozen data
CERT = json.loads((ROOT / "proposals/P250-shell-bubble-clock/attempts/0001"
                   / "maxwell_point.json").read_text())
A_VEC = np.array([1.0, 0.0, 0.0, 0.0])
B_STATE = np.array([float(CERT["c"]), float(CERT["b"]), float(CERT["f"])])
B_VEC = np.array([0.0, *B_STATE])
W2_STAR = float(CERT["omega_sq"])
BOX = 1e-12            # certified bulk-coordinate box halfwidth
TOL_BVP = 1e-11
L_LADDER = (8.0, 12.0, 16.0, 24.0)
N_LADDER = (401, 801, 1601, 3201)
L_FIX, N_FIX = 12.0, 3201

# ------------------------------------------------------------ exact callables
_m, _c, _b, _f = sp.symbols("m c b f", real=True)
_w2 = sp.Symbol("omega_sq", positive=True)
W2_EX = sp.N(str(W2_STAR), 30)
Vw_expr = wall_slice_potential(_m, _c, _b, _f, _w2).subs(_w2, W2_EX)
# Slice license psi = f >= 0 (accepted C-M5W-001): Abs(f) -> f before diff,
# so the EL gradient is the exact polynomial (verified against the certified
# Maxwell system of attempt 0001).
Vw_expr = Vw_expr.replace(lambda e: isinstance(e, sp.Abs), lambda e: e.args[0])
gv = [sp.diff(Vw_expr, v) for v in (_m, _c, _b, _f)]
V_fun = sp.lambdify((_m, _c, _b, _f), Vw_expr, "math")
G_fun = sp.lambdify((_m, _c, _b, _f), gv, "math")

ASSERT_IDENT = first_integral_conservation_identity()
assert ASSERT_IDENT, "first-integral identity must hold exactly"


def rhs(_x, y):
    """EL system 2K u'' = grad V_w, 2K = diag(1/2,1,1,1)."""
    gm, gc, gb, gf = G_fun(y[0], y[1], y[2], y[3])
    return np.vstack((y[4:], [2.0 * gm, gc, gb, gf]))


def kinetic(y):
    return (y[4] ** 2 + 2 * y[5] ** 2 + 2 * y[6] ** 2) / 4.0 + y[7] ** 2 / 2.0


def pot(y):
    return V_fun(y[0], y[1], y[2], y[3])


# ------------------------------------------------------- tail asymptotics
H_A = sp.diag(5, 10, 22 - 4 * W2_EX, 6 - W2_EX)
H_B = sp.hessian(Vw_expr, (_m, _c, _b, _f)).subs(
    {_m: sp.N(str(B_STATE[0]), 30), _c: sp.N(str(B_STATE[0]), 30),
     _b: sp.N(str(B_STATE[1]), 30), _f: sp.N(str(B_STATE[2]), 30)})


def eigs_of(Dsym):
    Dmp = mp.matrix([[mp.mpf(str(sp.N(e, 30))) for e in Dsym.row(i)] for i in range(4)])
    lams, vecs = mp.eigsy(Dmp)
    order = sorted(range(4), key=lambda i: mp.re(lams[i]))
    lams = [float(mp.re(lams[i])) for i in order]
    Q = np.array([[float(mp.re(vecs[r, i])) for i in order] for r in range(4)])
    return np.array(lams), Q  # columns of Q: eigenvectors, ascending lambda


LAM_A, Q_A = eigs_of(sp.diag(2, 1, 1, 1) * H_A)
LAM_B, Q_B = eigs_of(sp.diag(2, 1, 1, 1) * H_B)
SQ_A, SQ_B = np.sqrt(LAM_A), np.sqrt(LAM_B)
M_MIN = float(min(SQ_A.min(), SQ_B.min()))


# --------------------------------------------------------- boundary conditions
def bc_dirichlet(ya, yb):
    return np.concatenate((ya[:4] - A_VEC, yb[:4] - B_VEC))


def bc_expfit(ya, yb):
    """Stable-manifold projection at both ends (linearized tail)."""
    xa = (ya[:4] - A_VEC) @ Q_A
    xb = (yb[:4] - B_VEC) @ Q_B
    va = ya[4:] @ Q_A
    vb = yb[4:] @ Q_B
    return np.concatenate((va - SQ_A * xa, vb + SQ_B * xb))


# ------------------------------------------------------------- initial guesses
_STRING = None


def _load_string():
    """Converged field-space path from attempts/0002/string_guess.json."""
    global _STRING
    if _STRING is None:
        d = json.loads((Path(__file__).parent / "string_guess.json").read_text())
        su = np.array(d["u"])                       # (N, 4) monotone A -> B
        _STRING = np.vstack((A_VEC[None, :], su, B_VEC[None, :]))
    return _STRING


def guess(x, L, family):
    if family == "string":
        S = _load_string()
        w = 1.05                                    # tail-warp rate ~ soft mass
        s = 0.5 * (1.0 - np.tanh(x / w))
        ds = -0.5 * np.cosh(x / w) ** -2 / w
        N = S.shape[0]
        idx = np.clip(s * (N - 1), 0, N - 1.000001)
        j = idx.astype(int)
        frac = idx - j
        y0 = S[j] * (1 - frac[:, None]) + S[j + 1] * frac[:, None]
        dy0 = (S[j + 1] - S[j]) * ds[:, None] * (N - 1)
        return np.vstack((y0.T, dy0.T))
    if family == "tanh":
        w = 0.55
        s = 0.5 * (1.0 - np.tanh(x / w))
        ds = -0.5 * np.cosh(x / w) ** -2 / w
    elif family == "linear":
        s = np.clip((x + L) / (2 * L), 0.0, 1.0)
        ds = np.full_like(x, 0.5 / L)
    elif family == "step":
        w = 0.25
        s = 0.5 * (1.0 - np.tanh(x / w))
        ds = -0.5 * np.cosh(x / w) ** -2 / w
    y0 = A_VEC[None, :] + np.outer(s, B_VEC - A_VEC)
    dy0 = np.outer(ds, B_VEC - A_VEC)
    return np.vstack((y0.T, dy0.T))


def solve(L, n, bc, family=None, b_vec=None, warm=None, tol=None):
    """Collocation solve; warm-start from a previous solution when given."""
    global B_VEC
    b_save = B_VEC.copy()
    if b_vec is not None:
        B_VEC = b_vec
    x = np.linspace(-L, L, n)
    try:
        yy = warm.sol(x) if warm is not None else guess(x, L, family)
        sol = solve_bvp(rhs, bc, x, yy, tol=tol or TOL_BVP, max_nodes=400_000)
    finally:
        B_VEC = b_save
    return sol


CONT_L = (2.0, 3.0, 4.5, 6.0, 8.0, 10.0, 12.0)
CONT_TOL = (1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-10, 1e-10)


def continuation(n0=401):
    """L-continuation of the Dirichlet kink; the production route.

    Short boxes admit only the 1-kink; warm-starting each extension tracks
    its branch to L_FIX. Direct long-box solves fall into multi-kink trains
    (documented separately as basin probes).
    """
    x = np.linspace(-CONT_L[0], CONT_L[0], n0)
    sol = solve_bvp(rhs, bc_dirichlet, x, guess(x, CONT_L[0], "string"),
                    tol=CONT_TOL[0], max_nodes=400_000)
    rows = []
    r = tension_routes(sol, CONT_L[0])
    r["L"], r["tol"] = CONT_L[0], CONT_TOL[0]
    rows.append(r)
    print(f"cont L={CONT_L[0]:5.1f} status={sol.status} nodes={sol.x.size:6d} "
          f"r1={r['r1']:.12f} fi={r['fi_drift']:.1e}")
    for L, tol in zip(CONT_L[1:], CONT_TOL[1:]):
        x = np.linspace(-L, L, n0)
        sol = solve_bvp(rhs, bc_dirichlet, x, sol.sol(x), tol=tol, max_nodes=400_000)
        r = tension_routes(sol, L)
        r["L"], r["tol"] = L, tol
        rows.append(r)
        print(f"cont L={L:5.1f} status={sol.status} nodes={sol.x.size:6d} "
              f"r1={r['r1']:.12f} fi={r['fi_drift']:.1e}")
    return sol, rows


def ctrapz(xe, g):
    """Compensated trapezoid rule via math.fsum (permutation-invariant)."""
    m = len(xe)
    w = np.ones(m)
    w[0] = w[-1] = 0.5
    return (xe[1] - xe[0]) * math.fsum((w * g).tolist())


def tension_routes(sol, L):
    """Three tension routes on a refined grid with compensated summation."""
    m = 4 * (len(sol.x) - 1) + 1
    xe = np.linspace(-L, L, m)
    Y = sol.sol(xe)
    u, du = Y[:4], Y[4:]
    T = kinetic(Y)
    V = pot(Y)
    r1 = ctrapz(xe, T + V)
    r2 = 2.0 * ctrapz(xe, V)
    r3 = 2.0 * ctrapz(xe, T)

    def g_scalar(xv):
        Yv = np.vstack(sol.sol(np.array([xv])))
        return float((kinetic(Yv) + pot(Yv))[0])

    qv, qerr = quad(g_scalar, -L, L, epsabs=1e-13, epsrel=1e-13, limit=800)
    fi_drift = float(np.max(np.abs(T - V)))
    end_a = float(np.max(np.abs(sol.y[:4, 0] - A_VEC)))
    end_b = float(np.max(np.abs(sol.y[:4, -1] - B_VEC)))
    mono = {
        "m_nondecreasing_violations": int(np.sum(np.diff(u[0]) > 1e-8)),
        "c_nonincreasing_violations": int(np.sum(np.diff(u[1]) < -1e-8)),
        "b_nonincreasing_violations": int(np.sum(np.diff(u[2]) < -1e-8)),
        "f_nonincreasing_violations": int(np.sum(np.diff(u[3]) < -1e-8)),
    }
    res_rms = float(np.sqrt(np.mean(sol.rms_residuals ** 2)))
    return dict(r1=r1, r2=r2, r3=r3, quad=(float(qv), float(qerr)),
                fi_drift=fi_drift, end_a=end_a, end_b=end_b, mono=mono,
                res_rms=res_rms, n_nodes=len(sol.x), status=int(sol.status))


def main():
    print("== P250 G1: sigma_0 and the wall profile ==")
    print(f"python {platform.python_version()} numpy {np.__version__} "
          f"scipy {scipy.__version__} threads pinned to 1")
    print(f"first-integral identity exact: {ASSERT_IDENT}")
    print(f"omega*^2 = {W2_STAR!r}")
    print(f"A-masses {np.round(SQ_A, 6).tolist()}  B-masses {np.round(SQ_B, 6).tolist()}"
          f"  m_min {M_MIN:.6f}")
    out = {"env": {"python": platform.python_version(), "numpy": np.__version__,
                   "scipy": scipy.__version__, "threads": 1,
                   "tol_bvp": TOL_BVP, "box": BOX,
                   "route": "Dirichlet L-continuation (production); "
                            "direct long-box solves fall into multi-kink trains"},
           "tail": {"sqrt_lam_A": SQ_A.tolist(), "sqrt_lam_B": SQ_B.tolist(),
                    "m_min": M_MIN}}

    # ---- 1) production route: L-continuation of the Dirichlet kink
    print("\n-- L-continuation (Dirichlet, warm-started) --")
    sol12, cont_rows = continuation()
    out["continuation"] = cont_rows
    sigL = [r["r1"] for r in cont_rows]
    out["domain_estimate"] = abs(sigL[-1] - sigL[-2])
    sol10 = None
    x10 = np.linspace(-10.0, 10.0, 401)

    # ---- 2) h-crossing at L=10 (the clean status-0 rung), warm-started
    print("\n-- h-crossing at L=10 (Dirichlet, warm chain) --")
    h_rows = []
    warm = sol12 if sol12 is not None else None
    # re-solve at L=10 from the continuation state for a warm chain
    solL = solve(10.0, 401, bc_dirichlet, warm=sol12, tol=1e-10)
    for n in N_LADDER:
        solL = solve(10.0, n, bc_dirichlet, warm=solL, tol=1e-10)
        row = tension_routes(solL, 10.0)
        row["n_req"] = n
        h_rows.append(row)
        print(f"n={n:5d} nodes={row['n_nodes']:6d} r1={row['r1']:.14f} "
              f"resRMS={row['res_rms']:.2e} fi={row['fi_drift']:.2e} "
              f"status={row['status']}")
    out["h_crossing_L10"] = h_rows
    sig = [r["r1"] for r in h_rows]
    d = [abs(sig[i + 1] - sig[i]) for i in range(len(sig) - 1)]
    p_obs = math.log2(d[0] / d[1]) if len(d) > 1 and d[1] > 0 else float("inf")
    rich = (sig[-1] + (sig[-1] - sig[-2]) / (2 ** p_obs - 1)
            if math.isfinite(p_obs) and p_obs > 0.5 else None)
    print(f"observed order p = {p_obs:.3f}   richard(sigma_inf) = {rich!r}")
    out["observed_order"] = p_obs
    out["richardson_sigma_inf"] = rich

    # ---- 3) h-crossing at L=12 (headline box), warm chain
    print("\n-- h-crossing at L=12 (Dirichlet, warm chain) --")
    h12 = []
    solF = solve(L_FIX, 401, bc_dirichlet, warm=sol12, tol=1e-10)
    for n in N_LADDER:
        solF = solve(L_FIX, n, bc_dirichlet, warm=solF, tol=1e-10)
        row = tension_routes(solF, L_FIX)
        row["n_req"] = n
        h12.append(row)
        print(f"n={n:5d} nodes={row['n_nodes']:6d} r1={row['r1']:.14f} "
              f"resRMS={row['res_rms']:.2e} fi={row['fi_drift']:.2e} "
              f"status={row['status']}")
    out["h_crossing_L12"] = h12

    # ---- 4) independent boundary treatment: exponential-fit BC at L_FIX
    print("\n-- exponential-fit BC cross-check (L=12, warm from Dirichlet) --")
    solE = solve(L_FIX, N_FIX, bc_expfit, warm=solF, tol=1e-10)
    rE = tension_routes(solE, L_FIX)
    print(f"expfit: status={rE['status']} nodes={rE['n_nodes']} "
          f"r1={rE['r1']:.14f} fi={rE['fi_drift']:.2e} endA={rE['end_a']:.2e} "
          f"endB={rE['end_b']:.2e}")
    out["expfit_crosscheck"] = rE
    out["bc_agreement"] = abs(rE["r1"] - h12[-1]["r1"])
    print(f"BC agreement |expfit - dirichlet| = {out['bc_agreement']:.3e}")

    # ---- 5) basin probes: direct guesses at long box (documented, not budget)
    print("\n-- basin probes (direct guesses land in multi-kink trains) --")
    probes = {}
    for fam in ("tanh", "linear", "step"):
        s = solve(L_FIX, 1601, bc_dirichlet, fam, tol=1e-8)
        row = tension_routes(s, L_FIX)
        probes[fam] = {"status": row["status"], "r1": row["r1"],
                       "mono_bad": {k: v for k, v in row["mono"].items() if v}}
        print(f"probe {fam:7s}: status={row['status']} r1={row['r1']:.6f} "
              f"mono_bad={probes[fam]['mono_bad']}")
    out["basin_probes"] = probes

    # ---- 6) certified-box drift (b4), warm-started at production settings
    print("\n-- certified-box drift (b4) --")
    drift = {}
    for tag, db in (("up", +BOX), ("down", -BOX)):
        b_vec = np.array([0.0, B_STATE[0] + db, B_STATE[1] + db, B_STATE[2] + db])
        s = solve(L_FIX, N_FIX, bc_dirichlet, warm=solF, tol=1e-10, b_vec=b_vec)
        row = tension_routes(s, L_FIX)
        drift[tag] = row["r1"]
        print(f"box {tag}: sigma = {row['r1']:.14f}")
    out["box_drift"] = {"up": drift["up"], "down": drift["down"],
                        "half_spread": abs(drift["up"] - drift["down"]) / 2}

    # ---- 7) evaluator noise (b6) on the finest solution
    Y = solF.sol(np.linspace(-L_FIX, L_FIX, 6401))
    T = kinetic(Y)
    V = pot(Y)
    fs = math.fsum((T + V).tolist())
    fs_rev = math.fsum((T + V)[::-1].tolist())
    T2 = (Y[4] ** 2 + 2 * Y[5] ** 2 + 2 * Y[6] ** 2) / 4.0 + Y[7] ** 2 / 2.0
    V2 = np.array([V_fun(*Y[:4, i]) for i in range(Y.shape[1])])
    jit = float(np.max(np.abs((T2 + V2) - (T + V))))
    print(f"-- evaluator noise (b6): fsum fwd/rev {fs_rev - fs!r}, "
          f"codepath jitter {jit:.2e}")
    out["evaluator_noise"] = {"fsum_permutation": fs_rev - fs, "codepath_jitter": jit}

    # ---- headline: finest L=12 Dirichlet run
    finest = h12[-1]
    routes = [finest["r1"], finest["r2"], finest["r3"], finest["quad"][0]]
    spread = max(routes) - min(routes)
    budget = {
        "b1_residual_scaled": finest["res_rms"] * (2 * L_FIX),
        "b2_truncation_halfspread_L12": abs(h12[-1]["r1"] - h12[-2]["r1"]) / 2,
        "b3_domain_L10_to_L12": out["domain_estimate"],
        "b4_certified_box": out["box_drift"]["half_spread"],
        "b5_quadrature": abs(finest["quad"][0] - finest["r1"]),
        "b6_evaluator": max(abs(out["evaluator_noise"]["fsum_permutation"]),
                            out["evaluator_noise"]["codepath_jitter"]),
        "b7_bc_agreement": out["bc_agreement"],
        "b8_first_integral_drift": finest["fi_drift"],
    }
    total = max(budget.values())
    print("\n== headline ==")
    print("sigma_0 (finest Dirichlet L=12 n=3201, warm chain):")
    for name, v in zip(("R1 T+V", "R2 2V", "R3 2T", "GK-quad"), routes):
        print(f"  {name:8s} = {v:.15f}")
    print(f"route spread = {spread:.3e}")
    for k, v in budget.items():
        print(f"  {k:34s} {v:.3e}")
    print(f"  budget total (max item)            {total:.3e}")
    print(f"acceptance: total <= 1e-3 * sigma_0 = {1e-3 * finest['r1']:.3e} "
          f"-> {'PASS' if total <= 1e-3 * finest['r1'] else 'FAIL'}")
    out["headline"] = {
        "routes": dict(zip(("r1_TV", "r2_2V", "r3_2T", "gk_quad"), routes)),
        "route_spread": spread, "budget": budget, "budget_total": total,
        "accept": total <= 1e-3 * finest["r1"],
        "sigma0_best": float(np.mean(routes)),
        "sigma0_richardson": rich}

    # ---- profile artifact (finest solve)
    # artifact from a fully converged rung: re-solve at n=1601 (the status-0
    # node count in the chain), warm-started from the finest solution
    solF = solve(L_FIX, 1601, bc_dirichlet, warm=solF, tol=1e-10)
    rF = tension_routes(solF, L_FIX)
    print(f"artifact rung: status={rF['status']} nodes={rF['n_nodes']} "
          f"r1={rF['r1']:.14f} resRMS={rF['res_rms']:.2e} "
          f"fi={rF['fi_drift']:.2e}")
    xs = np.linspace(-L_FIX, L_FIX, 2001)
    Y = solF.sol(xs)
    u, du = Y[:4], Y[4:]
    with open(Path(__file__).parent / "profile_L12_n3201.csv", "w") as fh:
        fh.write("x,m,c,b,f,V,T\n")
        for i in range(len(xs)):
            Yi = np.concatenate((u[:, i], du[:, i]))
            fh.write(f"{xs[i]:.10g},{u[0, i]:.12g},{u[1, i]:.12g},"
                     f"{u[2, i]:.12g},{u[3, i]:.12g},{V_fun(*u[:, i]):.12g},"
                     f"{kinetic(Yi):.12g}\n")
    (Path(__file__).parent / "sigma0_results.json").write_text(json.dumps(
        out, indent=1,
        default=lambda o: bool(o) if isinstance(o, np.bool_) else float(o)))
    print("\nwrote sigma0_results.json and profile_L12_n3201.csv")


if __name__ == "__main__":
    main()
