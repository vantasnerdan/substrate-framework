"""PoC-2 archive: reduced thin-ring leapfrog — periodic orbit + Floquet + reduced Hessian.

Reproduces receipts/poc2-filament verdicts from eval cells "PoC-2 leapfrog demo"
and "PoC-2 Newton plus Floquet" (2026-09-10). Run: python3 run_poc2.py
Scope: REDUCED model (fixed core, axisymmetric). Full-3D filament + Euler shadowing OUT.
"""
import time
import numpy as np

G, a = 1.0, 0.05
NQ = 200
DT = 0.004


def self_V(R):
    return G / (4 * np.pi * R) * (np.log(8 * R / a) - 0.25)


def mutual(Rt, Zt, Rs, Zs, nq=NQ):
    ph = np.linspace(0, 2 * np.pi, nq, endpoint=False)
    dph = 2 * np.pi / nq
    sx, sy = Rs * np.cos(ph), Rs * np.sin(ph)
    rx, ry, rz = Rt - sx, -sy, Zt - Zs
    r = np.sqrt(rx * rx + ry * ry + rz * rz + a * a)
    dlx, dly = -Rs * np.sin(ph) * dph, Rs * np.cos(ph) * dph
    f = G / (4 * np.pi * r ** 3)
    return np.array([np.sum(f * dly * rz),
                     np.sum(f * (-dlx * rz)),
                     np.sum(f * (dlx * ry - dly * rx))])


def rhs2(s):
    R1, Z1, R2, Z2 = s
    m12 = mutual(R1, Z1, R2, Z2)
    m21 = mutual(R2, Z2, R1, Z1)
    return np.array([m12[0], self_V(R1) + m12[2],
                     m21[0], self_V(R2) + m21[2]])


def rk4(s, dt):
    k1 = rhs2(s)
    k2 = rhs2(s + 0.5 * dt * k1)
    k3 = rhs2(s + 0.5 * dt * k2)
    k4 = rhs2(s + dt * k3)
    return s + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def shape(s):
    R1, Z1, R2, Z2 = s
    return np.array([R1, R2, Z1 - Z2])


def full_map(x, dt=DT, nmax=4000):
    R1, R2 = x
    s = np.array([R1, -1e-9, R2, 1e-9])
    prevD = s[1] - s[3]
    t = 0.0
    ncr = 0
    for _ in range(nmax):
        s = rk4(s, dt)
        t += dt
        D = s[1] - s[3]
        if prevD < 0 and D >= 0:
            ncr += 1
            if ncr == 2:
                return shape(s), t
        prevD = D
    return None, None


def residual(x):
    sh, T = full_map(x)
    return sh - np.array([x[0], x[1], 0.0]), T


def flow_shape(sh0, T, dt=DT):
    s = np.array([sh0[0], sh0[2] / 2, sh0[1], -sh0[2] / 2])
    for _ in range(int(T / dt)):
        s = rk4(s, dt)
    return shape(s)


def main():
    t0 = time.time()
    # demo: section crossings
    s = np.array([1.0, -0.25, 1.0, 0.25])
    prev = s.copy()
    ncross = 0
    for i in range(int(40.0 / 0.005)):
        s = rk4(s, 0.005)
        if (prev[1] - prev[3]) * (s[1] - s[3]) < 0:
            ncross += 1
            if ncross <= 4:
                print(f"crossing t={i*0.005:.3f} R1={s[0]:.4f} R2={s[2]:.4f} Z={s[1]:.4f}")
        prev = s.copy()
    print(f"crossings in T=40: {ncross}")
    # Newton
    x = np.array([0.7729, 1.1843])
    T = None
    for it in range(12):
        r, T = residual(x)
        rn = np.linalg.norm(r)
        print(f"it{it}: |res|={rn:.3e} T={T:.5f} R1={x[0]:.6f} R2={x[1]:.6f}")
        if rn < 1e-9:
            break
        J = np.zeros((3, 2))
        e = 1e-6
        for j in range(2):
            dx = np.zeros(2)
            dx[j] = e
            r2, _ = residual(x + dx)
            J[:, j] = (r2 - r) / e
        dx, *_ = np.linalg.lstsq(J, -r, rcond=None)
        x = x + dx
    # monodromy (fixed-T flow on shape vars)
    sh_star = np.array([x[0], x[1], 0.0])
    F0 = flow_shape(sh_star, T)
    print("fixed-T return res:", np.linalg.norm(F0 - sh_star))
    M = np.zeros((3, 3))
    e = 1e-6
    for j in range(3):
        dp = np.zeros(3)
        dp[j] = e
        M[:, j] = (flow_shape(sh_star + dp, T) - flow_shape(sh_star - dp, T)) / (2 * e)
    ev = np.linalg.eigvals(M)
    print("monodromy eigs:", ev, "|.|=", np.abs(ev))
    # single-ring reduced Hessian
    R = 1.0
    h = 1e-5
    E = lambda R: 0.5 * G * G * R * (np.log(8 * R / a) - 1.75)
    I = lambda R: np.pi * G * R * R
    V = (E(R + h) - E(R - h)) / (I(R + h) - I(R - h))
    H2 = (E(R + h) - V * I(R + h) - 2 * (E(R) - V * I(R)) + E(R - h) - V * I(R - h)) / h ** 2
    print(f"single-ring: V={V:.5f} d2(E-VI)/dR2={H2:.5f} (Z: exact zero mode)")
    print(f"elapsed {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
