"""A3 archive: 3D filament resonance scan (FROZEN design A3-design.md, C1-C8).
Core model (C4 explicit): Saffman-local self-induction + Rosenhead-Moore mutual.
Scope: reduced filament model; verdicts -in-model; live-field gap uncrossed (C7).
Stages: gate | orbit | mono | ladder | control  (argv[1]; default gate)
"""
import sys, time
import numpy as np

import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "0108-cipher-radical", "receipts", "poc2-filament"))
G = 1.0
AA = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
N = int(sys.argv[3]) if len(sys.argv) > 3 else 64
DT = 0.005
PH = np.linspace(0, 2 * np.pi, N, endpoint=False)
DPH = 2 * np.pi / N
RUNIT = np.stack([np.cos(PH), np.sin(PH), np.zeros(N)], axis=1)   # radial unit
AUNIT = np.stack([-np.sin(PH), np.cos(PH), np.zeros(N)], axis=1)  # azimuthal unit
ZUNIT = np.tile(np.array([0.0, 0.0, 1.0]), (N, 1))


def ring_state(R1, Z1, R2, Z2):
    X = np.zeros((2, N, 3))
    for n, (R, Z) in enumerate([(R1, Z1), (R2, Z2)]):
        X[n, :, 0] = R * np.cos(PH)
        X[n, :, 1] = R * np.sin(PH)
        X[n, :, 2] = Z
    return X


def rhs3(X, Gam=G, aa=AA):
    V = np.zeros_like(X)
    for n in range(2):
        Xn = X[n]
        Rn = np.sqrt(Xn[:, 0] ** 2 + Xn[:, 1] ** 2).mean()
        V[n, :, 2] = Gam / (4 * np.pi * Rn) * (np.log(8 * Rn / aa) - 0.25)
        Xm = X[1 - n]
        d = Xn[:, None, :] - Xm[None, :, :]
        r2 = (d ** 2).sum(-1) + aa * aa
        T = np.zeros_like(Xm)
        T[:, 0] = -np.sin(PH)
        T[:, 1] = np.cos(PH)
        Rm = np.sqrt(Xm[:, 0] ** 2 + Xm[:, 1] ** 2).mean()
        V[n] += Gam / (4 * np.pi) * np.sum(
            np.cross(T[None, :, :], d) / (r2 ** 1.5)[:, :, None], axis=1) * (Rm * DPH)
    return V


def rk4_3(X, dt, **kw):
    k1 = rhs3(X, **kw)
    k2 = rhs3(X + dt / 2 * k1, **kw)
    k3 = rhs3(X + dt / 2 * k2, **kw)
    k4 = rhs3(X + dt * k3, **kw)
    return X + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def flow(X0, T, dt=DT, **kw):
    X = X0.copy()
    n = int(round(T / dt))
    for _ in range(n):
        X = rk4_3(X, dt, **kw)
    return X


# ---- per-m basis (C5): dirs per m; m=0 -> 6 real, m>0 -> 12 real ----
def basis_field(m, n, k, part):
    """unit perturbation field on ring n, direction k (0=rad,1=azi,2=ax), Re/Im part."""
    F = np.zeros((2, N, 3))
    U = [RUNIT, AUNIT, ZUNIT][k]
    if m == 0:
        F[n] = U
    else:
        F[n] = U * (np.cos(m * PH) if part == 0 else np.sin(m * PH))[:, None]
    return F


def ndirs(m):
    return 6 if m == 0 else 12


def dir_index(m, i):
    if m == 0:
        return (i // 3, i % 3, 0)
    return ((i // 6), (i // 2) % 3, i % 2)


def project(Xd):
    """project displacement field onto all (m,n,k,part) comps -> dict m: vec."""
    out = {}
    v0 = np.zeros(6)
    for n in range(2):
        for k, U in enumerate([RUNIT, AUNIT, ZUNIT]):
            v0[n * 3 + k] = (Xd[n] * U).sum() / N
    out[0] = v0
    for m in range(1, 9):
        v = np.zeros(12)
        for n in range(2):
            for k, U in enumerate([RUNIT, AUNIT, ZUNIT]):
                c = (Xd[n] * U)
                v[n * 6 + k * 2 + 0] = (c.sum(axis=1) * np.cos(m * PH)).sum() * 2 / N
                v[n * 6 + k * 2 + 1] = (c.sum(axis=1) * np.sin(m * PH)).sum() * 2 / N
        out[m] = v
    return out


def monodromy(X0, T, mmax=6, eps=1e-6, dt=DT, **kw):
    XB = flow(X0, T, dt=dt, **kw)
    pB = project(XB - axisym_of(XB))
    Mono = {}
    for m in range(0, mmax + 1):
        nd = ndirs(m)
        M = np.zeros((nd, nd))
        for i in range(nd):
            F = basis_field(m, *dir_index(m, i))
            XT = flow(X0 + eps * F, T, dt=dt, **kw)
            col = (project(XT - axisym_of(XT))[m] - pB[m]) / eps
            M[:, i] = col[:nd] if m > 0 else col
        Mono[m] = M
    return Mono, XB


def axisym_of(X):
    """axisymmetric part (m=0 reconstruction) of a state."""
    Y = np.zeros_like(X)
    for n in range(2):
        R = np.sqrt(X[n, :, 0] ** 2 + X[n, :, 1] ** 2).mean()
        Z = X[n, :, 2].mean()
        Y[n, :, 0] = R * np.cos(PH)
        Y[n, :, 1] = R * np.sin(PH)
        Y[n, :, 2] = Z
    return Y


def soft3_vectors(X0):
    """NAMED soft subspace SOFT3 (C6): time-shift + x/y translates as m=0/m=1 fields."""
    t = project(rhs3(X0))[0]          # time-shift tangent, m=0
    ex = np.zeros((2, N, 3)); ex[:, :, 0] = 1.0
    ey = np.zeros((2, N, 3)); ey[:, :, 1] = 1.0
    px = project(ex)  # lives in m=1
    py = project(ey)
    return t, px[1], py[1]


def stage_gate():
    _av, sys.argv = sys.argv, ["run_poc2.py"]
    from run_poc2 import mutual, rhs2
    sys.argv = _av
    V = rhs3(ring_state(1.0, 0.5, 1.0, -0.5))
    m = mutual(1.0, 0.5, 1.0, -0.5, nq=N)
    rad = V[:, :, 0] * np.cos(PH)[None, :] + V[:, :, 1] * np.sin(PH)[None, :]
    azi = -V[:, :, 0] * np.sin(PH)[None, :] + V[:, :, 1] * np.cos(PH)[None, :]
    ok = abs(rad[0, 0] - m[0]) < 1e-9 and abs(azi).max() < 1e-12
    print(f"gate radial {rad[0,0]:.12f} vs {m[0]:.12f}; azim {abs(azi).max():.1e}; GATE {'PASS' if ok else 'FAIL'}")

def stage_orbit():
    # banked PoC-2 section orbit; section return (Z1-Z2) + radii (pair drift is physical)
    R1, R2, T = 0.773723, 1.185226, 4.08800
    X0 = ring_state(R1, 0.0, R2, 0.0)
    t0 = time.time()
    XT = flow(X0, T)
    R = [np.sqrt(XT[n, :, 0] ** 2 + XT[n, :, 1] ** 2).mean() for n in range(2)]
    Z = [XT[n, :, 2].mean() for n in range(2)]
    print(f"orbit return: dR1={R[0]-R1:.2e} dR2={R[1]-R2:.2e} d(Z1-Z2)={Z[0]-Z[1]:.2e} "
          f"pair-drift={(Z[0]+Z[1])/2:.4f} nonaxisym={abs(XT - axisym_of(XT)).max():.1e} ({time.time()-t0:.1f}s)")


def stage_mono():
    R1, R2, T = 0.773723, 1.185226, 4.08800
    X0 = ring_state(R1, 0.0, R2, 0.0)
    t0 = time.time()
    Mono, XB = monodromy(X0, T, mmax=6)
    tS, px, py = soft3_vectors(X0)
    print(f"SOFT3 norms: |t-shift|={np.linalg.norm(tS):.3f} |px|={np.linalg.norm(px):.3f} |py|={np.linalg.norm(py):.3f}")
    for m in range(7):
        ev = np.linalg.eigvals(Mono[m])
        am = abs(ev)
        order = np.argsort(-am)
        top = " ".join(f"{am[j]:.6f}" for j in order[:4])
        grow = am[am > 1 + 1e-6]
        print(f"m={m}: |rho|max4 [{top}] n_grow={len(grow)}")
    # C1: W MEASURED — max relative convergence rate x Rbar^2/Gamma over base period
    X = X0.copy()
    n = int(round(T / DT))
    Wmax = 0.0
    for _ in range(n):
        R = [np.sqrt(X[k, :, 0] ** 2 + X[k, :, 1] ** 2).mean() for k in range(2)]
        Z = [X[k, :, 2].mean() for k in range(2)]
        V = rhs3(X)
        rR = [(V[k, :, 0] * np.cos(PH) + V[k, :, 1] * np.sin(PH)).mean() for k in range(2)]
        if abs(R[0] - R[1]) > 1e-9:
            Wmax = max(Wmax, abs(rR[0] - rR[1]) / abs(R[0] - R[1]) * np.mean(R) ** 2 / G)
        X = rk4_3(X, DT)
    Rm = float(np.mean([R1, R2]))
    print(f"W_MEASURED={Wmax:.4f} Lambda_Saff={np.log(8*Rm/AA):.3f} Lambda_ln={np.log(Rm/AA):.3f} ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1][0].isalpha() else "gate"
    {"gate": stage_gate, "orbit": stage_orbit, "mono": stage_mono}[stage]()
