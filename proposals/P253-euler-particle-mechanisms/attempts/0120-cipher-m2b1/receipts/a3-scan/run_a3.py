"""A3 archive: 3D filament resonance scan (FROZEN design A3-design.md, C1-C8 + D4).
Core model (C4 explicit): Saffman-local self-induction + Rosenhead-Moore mutual.
Scope: reduced filament model; verdicts -in-model; live-field gap uncrossed (C7).
Method history (honest trail): one-sided FD monodromy (R-A floor 1e-4 + eps-leg);
m=0 wrong-subspace block retired (measured non-axisym leak, found all-zero);
GAUGE2 named (azimuthal per-ring exact kernel); crossing-detect Newton retired
(branch jumps) -> smooth 3x3 single shooting (res 2.9e-12); lab-frame fixed-T m=0
Stages: gate | orbit | mono | m0 | newton3 | verdicts
Usage: python3 run_a3.py <stage> [aa [N]]
"""
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "0108-cipher-radical", "receipts", "poc2-filament"))
import numpy as np

G = 1.0
_args = sys.argv[1:]
_stage0 = _args[0] if _args and _args[0][:1].isalpha() else "gate"
_rest = _args[1:] if _args and _args[0] == _stage0 else _args
AA = float(_rest[0]) if len(_rest) > 0 else 0.05
N = int(_rest[1]) if len(_rest) > 1 else 64
DT = 0.005
PH = np.linspace(0, 2 * np.pi, N, endpoint=False)
DPH = 2 * np.pi / N
RUNIT = np.stack([np.cos(PH), np.sin(PH), np.zeros(N)], axis=1)
AUNIT = np.stack([-np.sin(PH), np.cos(PH), np.zeros(N)], axis=1)
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


def basis_field(m, n, k, part):
    """unit perturbation on ring n, dir k (0=rad,1=azi,2=ax), Re/Im part."""
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
    """project displacement field onto (m,n,k,part) comps -> dict m: vec."""
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


def axisym_of(X):
    Y = np.zeros_like(X)
    for n in range(2):
        R = np.sqrt(X[n, :, 0] ** 2 + X[n, :, 1] ** 2).mean()
        Z = X[n, :, 2].mean()
        Y[n, :, 0] = R * np.cos(PH)
        Y[n, :, 1] = R * np.sin(PH)
        Y[n, :, 2] = Z
    return Y


def monodromy(X0, T, mmax=6, eps=1e-6, dt=DT, **kw):
    """per-m FD monodromy. m=0 block RETIRED (wrong subspace) — Quasi-static scope:
    use only m>=1 columns; m=0 covered by stage_m0 (centered section map)."""
    XB = flow(X0, T, dt=dt, **kw)
    pB = project(XB - axisym_of(XB))
    Mono = {}
    for m in range(1, mmax + 1):
        nd = ndirs(m)
        M = np.zeros((nd, nd))
        for i in range(nd):
            F = basis_field(m, *dir_index(m, i))
            XT = flow(X0 + eps * F, T, dt=dt, **kw)
            col = (project(XT - axisym_of(XT))[m] - pB[m]) / eps
            M[:, i] = col[:nd]
        Mono[m] = M
    np.savez("Mono_e%s_m%02d.npz" % (eps, mmax),
             **{"m%d" % m: Mono[m] for m in Mono})
    return Mono, XB


def soft3_vectors(X0):
    t = project(rhs3(X0))[0]
    ex = np.zeros((2, N, 3))
    ex[:, :, 0] = 1.0
    ey = np.zeros((2, N, 3))
    ey[:, :, 1] = 1.0
    px = project(ex)
    py = project(ey)
    return t, px[1], py[1]


def stage_gate():
    _av, sys.argv = sys.argv, ["run_poc2.py"]
    from run_poc2 import mutual
    sys.argv = _av
    V = rhs3(ring_state(1.0, 0.5, 1.0, -0.5))
    m = mutual(1.0, 0.5, 1.0, -0.5, nq=N)
    rad = V[:, :, 0] * np.cos(PH)[None, :] + V[:, :, 1] * np.sin(PH)[None, :]
    azi = -V[:, :, 0] * np.sin(PH)[None, :] + V[:, :, 1] * np.cos(PH)[None, :]
    ok = abs(rad[0, 0] - m[0]) < 1e-9 and abs(azi).max() < 1e-12
    print(f"gate radial {rad[0,0]:.12f} vs {m[0]:.12f}; "
          f"azim {abs(azi).max():.1e}; GATE {'PASS' if ok else 'FAIL'}")


def shape4(X):
    return np.array([np.sqrt(X[0, :, 0] ** 2 + X[0, :, 1] ** 2).mean(), X[0, :, 2].mean(),
                     np.sqrt(X[1, :, 0] ** 2 + X[1, :, 1] ** 2).mean(), X[1, :, 2].mean()])


def shoot(R1g=0.773723, R2g=1.185226, Tg=4.08800, quiet=True):
    """smooth 3x3 single shooting: unknowns (R1,R2,T), section start Z1=Z2=0."""
    y = np.array([R1g, R2g, Tg])

    def res3(yy):
        X0 = ring_state(yy[0], 0.0, yy[1], 0.0)
        XT = flow(X0, yy[2])
        S = shape4(XT)
        return np.array([S[0] - yy[0], S[2] - yy[1], S[1] - S[3]])

    r = res3(y)
    for it in range(12):
        rn = np.linalg.norm(r)
        if not quiet:
            print(f"  shoot it{it}: |res|={rn:.3e} "
                  f"R1={y[0]:.6f} R2={y[1]:.6f} T={y[2]:.5f}", flush=True)
        if rn < 1e-10:
            break
        J = np.zeros((3, 3))
        e = np.array([1e-6, 1e-6, 1e-7])
        for j in range(3):
            dy = np.zeros(3)
            dy[j] = e[j]
            J[:, j] = (res3(y + dy) - r) / e[j]
        dy, *_ = np.linalg.lstsq(J, -r, rcond=None)
        y = y + dy
        r = res3(y)
    if not quiet:
        print(f"SHOOTING converged: R1={y[0]:.6f} R2={y[1]:.6f} "
              f"T={y[2]:.5f} |res|={np.linalg.norm(r):.2e}")
    return y[0], y[1], y[2], np.linalg.norm(r)


def stage_newton3():
    R1, R2, T, rn = shoot(quiet=False)
    print(f"NEWTON3: R1={R1:.6f} R2={R2:.6f} T={T:.5f} |res|={rn:.1e}")


def stage_orbit():
    R1, R2, T, rn = shoot()
    X0 = ring_state(R1, 0.0, R2, 0.0)
    t0 = time.time()
    XT = flow(X0, T)
    R = [np.sqrt(XT[n, :, 0] ** 2 + XT[n, :, 1] ** 2).mean() for n in range(2)]
    Z = [XT[n, :, 2].mean() for n in range(2)]
    print(f"orbit return: dR1={R[0]-R1:.2e} dR2={R[1]-R2:.2e} d(Z1-Z2)={Z[0]-Z[1]:.2e} "
          f"pair-drift={(Z[0]+Z[1])/2:.4f} nonaxisym={abs(XT - axisym_of(XT)).max():.1e} "
          f"({time.time()-t0:.1f}s)")


def stage_mono():
    R1, R2, T, rn = shoot()
    print(f"mono base: shot R1={R1:.6f} R2={R2:.6f} T={T:.5f} |res|={rn:.1e}", flush=True)
    X0 = ring_state(R1, 0.0, R2, 0.0)
    t0 = time.time()
    Mono, XB = monodromy(X0, T, mmax=6)
    tS, px, py = soft3_vectors(X0)
    print(f"SOFT3 norms: |t-shift|={np.linalg.norm(tS):.3f} "
          f"|px|={np.linalg.norm(px):.3f} |py|={np.linalg.norm(py):.3f}")
    for m in range(1, 7):
        ev = np.linalg.eigvals(Mono[m])
        am = abs(ev)
        order = np.argsort(-am)
        top = " ".join(f"{am[j]:.6f}" for j in order[:4])
        # R-B: RAW counts, UNLICENSED (no SOFT3 deflation) — do not consume; see verdicts
        print(f"m={m}: RAW|rho|max4 [{top}] UNLICENSED-nodeflate")
    # C1: W MEASURED — dimensionless velocity-gradient operator norm along orbit
    rng = np.random.default_rng(0)
    X = X0.copy()
    n = int(round(T / DT))
    Wmax = 0.0
    for step in range(n + 1):
        if step % max(n // 8, 1) == 0:
            Rm = float(np.mean([np.sqrt(X[k, :, 0] ** 2 + X[k, :, 1] ** 2).mean()
                                for k in range(2)]))
            V0 = rhs3(X)
            for _ in range(6):
                F = rng.normal(size=X.shape)
                F /= np.sqrt((F ** 2).sum() / F.size)
                dV = (rhs3(X + 1e-7 * F) - V0) / 1e-7
                Wmax = max(Wmax, np.sqrt((dV ** 2).sum() / dV.size) * Rm ** 2 / G)
        X = rk4_3(X, DT)
    Rm = float(np.mean([R1, R2]))
    print(f"W_MEASURED={Wmax:.4f} Lambda_Saff={np.log(8*Rm/AA):.3f} "
          f"Lambda_ln={np.log(Rm/AA):.3f} ({time.time()-t0:.1f}s)")


def stage_m0():
    # CENTERED section-map 3x3 monodromy at SHOT orbit (mirrors PoC-2 flow_shape).
    R1, R2, T, rn = shoot()
    print(f"m0 base: shot R1={R1:.6f} R2={R2:.6f} T={T:.5f} |res|={rn:.1e}", flush=True)
    t0 = time.time()

    def secflow(sh):
        # sh=(R1,R2,Zd): start Z1=+Zd/2,Z2=-Zd/2, flow fixed T, return shape
        X0 = ring_state(sh[0], sh[2] / 2, sh[1], -sh[2] / 2)
        XT = flow(X0, T)
        R = [np.sqrt(XT[k, :, 0] ** 2 + XT[k, :, 1] ** 2).mean() for k in range(2)]
        Z = [XT[k, :, 2].mean() for k in range(2)]
        return np.array([R[0], R[1], Z[0] - Z[1]])

    sh_star = np.array([R1, R2, 0.0])
    F0 = secflow(sh_star)
    print(f"m0 section return res: {np.linalg.norm(F0 - sh_star):.1e}")
    M = np.zeros((3, 3))
    e = 1e-6
    for j in range(3):
        dp = np.zeros(3)
        dp[j] = e
        M[:, j] = (secflow(sh_star + dp) - secflow(sh_star - dp)) / (2 * e)
    ev = np.linalg.eigvals(M)
    print("m=0 section-flow eigs:", " ".join(f"{v.real:.6f}{v.imag:+.6f}j" for v in ev),
          "|.|=", " ".join(f"{abs(v):.6f}" for v in ev))
    print(f"(expect PoC-2 0.9275+-0.3738i + 1; GAUGE2 kernel by construction) "
          f"({time.time()-t0:.1f}s)")


def stage_verdicts():
    # R-A floor 1e-4 + eps-leg; R-B SOFT3 deflation (m=1); m=0 via stage_m0.
    R1, R2, T, rn = shoot()
    print(f"verdicts base: shot R1={R1:.6f} R2={R2:.6f} T={T:.5f} |res|={rn:.1e}", flush=True)
    X0 = ring_state(R1, 0.0, R2, 0.0)
    D = np.load("Mono_e1e-06_m06.npz")
    tS, px, py = soft3_vectors(X0)
    Qp = np.stack([px, py], axis=1)
    Qp, _ = np.linalg.qr(Qp)
    print("SOFT3: m=1 dim2 (x/y-translate), orthonormalized; m=0 via stage_m0")
    for m in range(1, 7):
        M = D["m%d" % m]
        ev, EV = np.linalg.eig(M)
        am = abs(ev)
        if m == 1:
            P = np.eye(M.shape[0]) - Qp @ Qp.T
            amd = abs(np.linalg.eigvals(P @ M @ P))
            ov = [float(np.linalg.norm(Qp.T @ (EV[:, j] / np.linalg.norm(EV[:, j]))))
                  for j in range(len(ev))]
            nsoft = sum(1 for o in ov if o > 0.5)
        else:
            amd = am
            nsoft = 0
        top = " ".join(f"{a:.6f}" for a in sorted(amd)[-4:][::-1])
        ng = int((amd > 1 + 1e-4).sum())
        print(f"m={m}: deflated|rho| [{top}] n_grow(1e-4)={ng} soft-attrib={nsoft}")
    for m in (1, 2):
        M1 = D["m%d" % m]
        Mono2, _ = monodromy(X0, T, mmax=m, eps=5e-7)
        M2 = Mono2[m]
        print(f"m={m} eps-leg: max|M(eps)-M(eps/2)|={abs(M1 - M2).max():.2e}")


if __name__ == "__main__":
    {"gate": stage_gate, "orbit": stage_orbit, "mono": stage_mono, "m0": stage_m0,
     "newton3": stage_newton3, "verdicts": stage_verdicts}[_stage0]()
