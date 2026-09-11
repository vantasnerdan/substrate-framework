"""S4 tangle pilot (FROZEN F-bar in ../../05-s4/00-fbar.md): random-ring vacuum
tangle + unit-circle carrier; Gauss linking statistics over seeds; reports
(mu, sigma, m1, CV, Rb-scaling, Lambda-exponent, rotation spread) + verdict.
Usage: python3 run_tangle.py  (S=32 seeds, ~seconds)
"""
import time

import numpy as np

RNG = np.random.default_rng(7)
MC, K = 64, 16
RB = 3.0
RT = 0.3


def ring_pts(center, R, n, u, v):
    ph = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return center[None, :] + R * (np.cos(ph)[:, None] * u[None, :] + np.sin(ph)[:, None] * v[None, :])


def rand_frame(rng):
    A = rng.normal(size=(3, 3))
    Q, _ = np.linalg.qr(A)
    return Q[0] / np.linalg.norm(Q[0]), Q[1] / np.linalg.norm(Q[1])


def carrier(Rot=np.eye(3)):
    ph = np.linspace(0, 2 * np.pi, MC, endpoint=False)
    P = np.stack([np.cos(ph), np.sin(ph), np.zeros(MC)], axis=1) @ Rot.T
    return P


def gauss_link(P, Q):
    """Discrete Gauss integral over closed polylines P (M segs), Q (K segs)."""
    Pm = (P + np.roll(P, -1, axis=0)) / 2
    dl = np.roll(P, -1, axis=0) - P
    Qm = (Q + np.roll(Q, -1, axis=0)) / 2
    dm = np.roll(Q, -1, axis=0) - Q
    R = Pm[:, None, :] - Qm[None, :, :]
    r3 = np.linalg.norm(R, axis=2) ** 3 + 1e-30
    cross = np.cross(dl[:, None, :], dm[None, :, :])
    return float(((R * cross).sum(axis=2) / r3).sum() / (4 * np.pi))


def tangle(rng, Nt, Rb, Rt):
    rings = []
    for _ in range(Nt):
        c = rng.normal(size=3)
        c *= Rb * rng.random() ** (1 / 3) / max(np.linalg.norm(c), 1e-12)
        u, v = rand_frame(rng)
        rings.append(ring_pts(c, Rt, K, u, v))
    return rings


def stats(C, rings):
    return sum(gauss_link(C, Q) for Q in rings)


def run_case(Rb, Nt, S, Rot=np.eye(3)):
    C = carrier(Rot)
    Lks = np.zeros(S)
    for s in range(S):
        rng = np.random.default_rng(1000 + s)
        Lks[s] = stats(C, tangle(rng, Nt, Rb, RT))
    return Lks


def main():
    t0 = time.time()
    S = 32
    Nt = 40  # Lambda ~ 40*2pi*0.3/(4pi*27/3) ~= 0.67
    L = run_case(RB, Nt, S)
    mu, sd, m1 = L.mean(), L.std(), np.abs(L).mean()
    cv = L.std() / m1 if m1 > 0 else float("inf")
    print(f"base Rb={RB} Nt={Nt}: mu={mu:+.4f} sd={sd:.4f} m1={m1:.4f} CV={cv:.3f}")
    print(f"F-S4a (generator): {'PASS' if abs(mu) <= 2 * sd / np.sqrt(S) else 'FAIL-biased'}")
    L5 = run_case(5.0, int(Nt * (5 / RB) ** 3), S)
    m15 = np.abs(L5).mean()
    print(f"Rb=5: m1={m15:.4f} ratio={m15 / m1:.3f} -> {'PASS' if 0.5 <= m15 / m1 <= 2.0 else 'KILL-bnd'}")
    L2 = run_case(RB, 2 * Nt, S)
    m12 = np.abs(L2).mean()
    import math
    alpha = math.log(m12 / m1) / math.log(2) if m1 > 0 and m12 > 0 else float("nan")
    print(f"2xLambda: m1={m12:.4f} alpha={alpha:.3f} -> {'PASS-lean' if alpha <= 0.25 else ('KILL-lean' if alpha > 0.5 else 'gray')}")
    R_rots = [np.eye(3)]
    Rx = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
    Ry = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
    R_rots += [Rx, Ry, Rx @ Ry]
    qs = []
    for R in R_rots:
        LL = run_case(RB, Nt, S, R)
        qs.append(np.abs(LL).mean())
    qs = np.array(qs)
    print(f"rotations: {np.round(qs, 4)} spread={qs.std():.4f} (2sig={2 * sd / np.sqrt(S):.4f}) -> {'PASS' if qs.std() <= 2 * sd / np.sqrt(S) else 'KILL-frame'}")
    kills = [cv > 0.50, not (0.5 <= m15 / m1 <= 2.0), alpha > 0.5,
             qs.std() > 2 * sd / np.sqrt(S)]
    passes = [cv <= 0.25, 0.5 <= m15 / m1 <= 2.0, alpha <= 0.25,
              qs.std() <= 2 * sd / np.sqrt(S)]
    print(f"F-bar verdict: {'KILL' if any(kills[1:]) or kills[0] else ('FIT-PASS' if all(passes) else 'UNRESOLVED-gray')} ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
