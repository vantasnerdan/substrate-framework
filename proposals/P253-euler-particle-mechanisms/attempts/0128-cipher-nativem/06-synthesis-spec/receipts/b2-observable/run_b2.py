"""B2 dynamical observable (FROZEN F-bar in ../../06-synthesis-spec/03-b2-fbar.md):
helical filaments L=1,2,3 (plus far reference ring), Magnus test-force,
short-window 3D response vs full Lorentz template (v x B_mutual); slope
ratios vs integer ratios 2/1, 3/1, 3/2.
Usage: python3 run_b2.py [gate|fit] [N]
Self-contained 3D Biot-Savart (Saffman-local self + RM mutual, A3 backbone).
"""
import sys
import time

import numpy as np

G = 1.0
AA = 0.05
AH = 0.1
TWIN = 0.5
EPS = 1e-5
DT = 0.005


def ring_angles(N):
    return np.linspace(0, 2 * np.pi, N, endpoint=False)


def helix(L, N, R0=1.0):
    ph = ring_angles(N)
    X = np.zeros((2, N, 3))
    X[0, :, 0] = (R0 + AH * np.cos(L * ph)) * np.cos(ph)
    X[0, :, 1] = (R0 + AH * np.cos(L * ph)) * np.sin(ph)
    X[0, :, 2] = AH * np.sin(L * ph)
    X[1, :, 0] = R0 * np.cos(ph)
    X[1, :, 1] = R0 * np.sin(ph)
    X[1, :, 2] = -3.0
    return X


def tangent(X):
    T = np.roll(X, -1, axis=1) - np.roll(X, 1, axis=1)
    return T / (np.linalg.norm(T, axis=2, keepdims=True) + 1e-30)


def vel_decomp(X, Gam=1.0, aa=AA):
    """returns (V_total, V_mutual); V_self = axial Saffman (local)."""
    N = X.shape[1]
    ph = ring_angles(N)
    dph = 2 * np.pi / N
    Vt = np.zeros_like(X)
    Vm = np.zeros_like(X)
    for n in range(2):
        Xn = X[n]
        Rn = np.sqrt(Xn[:, 0] ** 2 + Xn[:, 1] ** 2).mean()
        Vt[n, :, 2] = Gam / (4 * np.pi * Rn) * (np.log(8 * Rn / aa) - 0.25)
        Xm = X[1 - n]
        d = Xn[:, None, :] - Xm[None, :, :]
        r2 = (d ** 2).sum(-1) + aa * aa
        Tm = np.zeros_like(Xm)
        Tm[:, 0] = -np.sin(ph)
        Tm[:, 1] = np.cos(ph)
        Rm = np.sqrt(Xm[:, 0] ** 2 + Xm[:, 1] ** 2).mean()
        Vm[n] = Gam / (4 * np.pi) * np.sum(
            np.cross(Tm[None, :, :], d) / (r2 ** 1.5)[:, :, None], axis=1) * (Rm * dph)
        # NOTE: source tangent circular approx (helix corrections O(AH/R0)).
        Vt[n] += Vm[n]
    return Vt, Vm


def rhs(X, eps):
    Vt, Vm = vel_decomp(X)
    return Vt + eps * np.cross(tangent(X), Vt)


def rk4(X, dt, eps):
    f = lambda q: rhs(q, eps)
    k1 = f(X)
    k2 = f(X + dt / 2 * k1)
    k3 = f(X + dt / 2 * k2)
    k4 = f(X + dt * k3)
    return X + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def leg(L, eps, N, dt=DT):
    """returns (slope kappa, residual, response, template) for carrier ring."""
    X0 = helix(L, N)
    n = int(round(TWIN / dt))
    Xb, Xp = X0.copy(), X0.copy()
    Rgr = np.zeros((n + 1,) + X0[0:1].shape)
    Tgr = np.zeros((n + 1,) + X0[0:1].shape)
    for i in range(n + 1):
        Vt, Vm = vel_decomp(Xb)
        Tgr[i, 0] = np.cross(Vt[0], Vm[0])
        if i < n:
            Rnow = (Xp - Xb) / eps
            Rgr[i, 0] = Rnow[0]
            Xb = rk4(Xb, dt, 0.0)
            Xp = rk4(Xp, dt, eps)
    Rg = Rgr[:, 0].reshape(-1)
    Tg = Tgr[:, 0].reshape(-1)
    k = float((Rg * Tg).sum() / (Tg * Tg).sum())
    r = float(np.sqrt(((Rg - k * Tg) ** 2).sum() / (Rg ** 2).sum()))
    return k, r


def stage_gate():
    k1, _ = leg(1, EPS, 64)
    k2, _ = leg(1, EPS / 2, 64)
    d = abs(k1 - k2) / max(abs(k1), 1e-12)
    print(f"eps-linearity (slope halves): {d:.2e} -> {'PASS' if d < 1e-2 else 'FAIL'}")


def stage_fit():
    t0 = time.time()
    N = int(sys.argv[3]) if len(sys.argv) > 3 else 64
    out = {}
    for L in (1, 2, 3):
        k, r = leg(L, EPS, N)
        out[L] = (k, r)
        print(f"L={L}: kappa={k:.4f} resid={r:.4f}", flush=True)
    ratios = {(2, 1): out[2][0] / out[1][0], (3, 1): out[3][0] / out[1][0],
              (3, 2): out[3][0] / out[2][0]}
    ok_r = all(abs(v - a / b) / (a / b) <= 1.0 for (a, b), v in ratios.items())
    res_ok = all(v[1] <= 0.25 for v in out.values())
    res_kill = all(v[1] > 0.50 for v in out.values())
    print("ratios:", " ".join(f"{a}/{b}={v:.3f}(exp {a/b})" for (a, b), v in ratios.items()))
    if res_ok and ok_r:
        vd = "PASS-observable-supplied"
    elif res_kill:
        vd = "KILL-observable-absent"
    else:
        vd = "UNRESOLVED-gray"
    print(f"B2 verdict: {vd} ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    st = sys.argv[1] if len(sys.argv) > 1 else "gate"
    {"gate": stage_gate, "fit": stage_fit}[st]()
