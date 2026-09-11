"""S1 kappa-fit (FROZEN F1 bar in ../../03-prefire.md): reduced axisymmetric
two-ring leapfrog + Magnus-form perturbation (fixed physics coefficient).
Test: period-INTEGRATED response R_P vs INSTANTANEOUS Lorentz template T_P
(t_hat x mutual-field on the BARE orbit); fit scalar slope kappa_P per
passage; F1 bar on (r_A, r_B, kappa ratio). Non-circular: template uses bare
kinematics only; fit tests whether integration preserves single-slope form.
Stages: gate (eps-linearity) | fit (passages A+B)
Usage: python3 run_kappa.py [gate|fit]
"""
import sys
import time

sys.path.insert(0, "/home/dan/substrate-framework/proposals/P253-euler-particle-mechanisms/attempts/0120-cipher-m2b1/receipts/a3-scan")
import numpy as np

from run_a3 import mut_a

G = 1.0
AA = 0.05
NQ = 200
DT = 0.004
RHO = 1.0


def self_V(R):
    return G / (4 * np.pi * R) * (np.log(8 * R / AA) - 0.25)


def rhs2m(s, eps):
    R1, Z1, R2, Z2 = s
    m12 = mut_a(R1, Z1, R2, Z2, G, AA, NQ)
    m21 = mut_a(R2, Z2, R1, Z1, G, AA, NQ)
    V1 = np.array([m12[0], self_V(R1) + m12[2]])
    V2 = np.array([m21[0], self_V(R2) + m21[2]])
    # Magnus slip = MUTUAL field only (self comoves by definition; docstring).
    # t_hat x (Wr, Wz) = (Wz, -Wr).
    a1 = eps * np.array([m12[2], -m12[0]])
    a2 = eps * np.array([m21[2], -m21[0]])
    return np.array([V1[0] + a1[0], V1[1] + a1[1], V2[0] + a2[0], V2[1] + a2[1]])


def rk4(s, dt, eps):
    f = lambda q: rhs2m(q, eps)
    k1 = f(s)
    k2 = f(s + dt / 2 * k1)
    k3 = f(s + dt / 2 * k2)
    k4 = f(s + dt * k3)
    return s + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def passage_data(R1b, R2b, T, eps):
    s0 = np.array([R1b, 0.0, R2b, 0.0])
    n = int(round(T / DT))
    sb, sp = s0.copy(), s0.copy()
    Rg = np.zeros((4, n + 1))
    Tg = np.zeros((4, n + 1))
    for i in range(n + 1):
        R1, Z1, R2, Z2 = sb
        m12 = mut_a(R1, Z1, R2, Z2, G, AA, NQ)
        m21 = mut_a(R2, Z2, R1, Z1, G, AA, NQ)
        # template: t_hat x mutual-field per ring (R,Z comps)
        Tg[:, i] = [m12[2], -m12[0], m21[2], -m21[0]]
        if i < n:
            sb = rk4(sb, DT, 0.0)
            sp = rk4(sp, DT, eps)
            Rg[:, i + 1] = (sp - sb) / eps
    # align: Rg[:,0]=0 (initial diff zero); template full grid
    return Rg, Tg


def stage_gate():
    R1b, R2b, T = 0.773723, 1.185226, 4.08800
    R1, _ = passage_data(R1b, R2b, T, 2.5e-5)
    R2, _ = passage_data(R1b, R2b, T, 1.25e-5)
    d = abs(R1 - R2).max() / max(abs(R1).max(), 1e-12)
    print(f"eps-linearity: rel diff halves {d:.2e} -> {'PASS' if d < 1e-3 else 'FAIL (reduce eps)'}")


def fit_one(Rg, Tg):
    k = float((Rg * Tg).sum() / (Tg * Tg).sum())
    r = float(np.sqrt(((Rg - k * Tg) ** 2).sum() / (Rg ** 2).sum()))
    return k, r


def stage_fit():
    t0 = time.time()
    eps = 2.5e-5
    R1b, R2b, T = 0.773723, 1.185226, 4.08800
    RA, TA = passage_data(R1b, R2b, T, eps)
    s = np.array([R1b, 0.0, R2b, 0.0])
    for _ in range(int(round(T / 4 / DT))):
        s = rk4(s, DT, 0.0)
    R1c, R2c = s[0], s[2]
    RB, TB = passage_data(R1c, R2c, T, eps)
    kA, rA = fit_one(RA, TA)
    kB, rB = fit_one(RB, TB)
    print(f"passage A: kappa={kA:.4f} r={rA:.4f}; passage B ({R1c:.4f},{R2c:.4f}): kappa={kB:.4f} r={rB:.4f}")
    ratio = max(kA, kB) / min(kA, kB) if min(kA, kB) > 0 else float("inf")
    verdict = "FIT-PASS" if max(rA, rB) <= 0.10 and ratio <= 2.0 else (
        "KILL" if max(rA, rB) > 0.25 or ratio > 3.0 or (kA < 0 or kB < 0) else "UNRESOLVED-gray")
    print(f"F1 verdict: {verdict} (ratio={ratio:.3f}) ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    st = sys.argv[1] if len(sys.argv) > 1 else "gate"
    {"gate": stage_gate, "fit": stage_fit}[st]()
