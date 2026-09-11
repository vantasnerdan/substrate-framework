"""S4b ordered-vacuum pilot (FROZEN F-bar in ../../10-s4b/00-fbar.md, committed
pre-compute): polarized ring tangle (cone half-angle th0) + carrier; ladder
th0 in {pi, pi/2, pi/4, pi/8} + frame-rotation leg; verdict per bar.
Usage: python3 run_s4b.py
"""
import os
import sys
import time

sys.path.insert(0, "/home/dan/substrate-framework/proposals/P253-euler-particle-mechanisms/attempts/0128-cipher-nativem/05-s4/receipts/tangle-pilot")
import numpy as np
import run_tangle as T

S = 32
RB, NT, RT = 3.0, 40, 0.3


def polarized_tangle(rng, Nt, Rb, Rt, th0):
    rings = []
    for _ in range(Nt):
        c = rng.normal(size=3)
        c *= Rb * rng.random() ** (1 / 3) / max(np.linalg.norm(c), 1e-12)
        # normal uniform within cone th0 around +z
        cos_t = 1 - rng.random() * (1 - np.cos(th0))
        sin_t = np.sqrt(max(0.0, 1 - cos_t ** 2))
        phi = rng.uniform(0, 2 * np.pi)
        n = np.array([sin_t * np.cos(phi), sin_t * np.sin(phi), cos_t])
        a = np.array([1.0, 0, 0]) if abs(n[0]) < 0.9 else np.array([0, 1.0, 0])
        u = a - a.dot(n) * n
        u /= np.linalg.norm(u)
        v = np.cross(n, u)
        rings.append(T.ring_pts(c, Rt, T.K, u, v))
    return rings


def run_case(th0, S, Rot=np.eye(3)):
    C = T.carrier(Rot)
    L = np.zeros(S)
    for s in range(S):
        rng = np.random.default_rng(3000 + s)
        L[s] = T.stats(C, polarized_tangle(rng, NT, RB, RT, th0))
    return L


def main():
    t0 = time.time()
    rows = {}
    for name, th0 in [("pi", np.pi), ("pi/2", np.pi / 2), ("pi/4", np.pi / 4),
                      ("pi/8", np.pi / 8)]:
        L = run_case(th0, S)
        mu, sd, m1 = L.mean(), L.std(), np.abs(L).mean()
        cv = sd / m1 if m1 > 0 else float("inf")
        rows[name] = (mu, sd, m1, cv)
        print(f"th0={name}: mu={mu:+.4f} sd={sd:.4f} m1={m1:.4f} CV={cv:.3f}", flush=True)
    m1s = [rows[k][2] for k in ("pi", "pi/2", "pi/4", "pi/8")]
    mono = all(b >= a for a, b in zip(m1s, m1s[1:]))
    print(f"monotone-coherent: {mono}")
    # frame leg at th0=pi/4: rotate carrier out of alignment plane
    Rx = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
    Lf = run_case(np.pi / 4, S, Rx @ np.eye(3))
    print(f"frame-rotated: m1={np.abs(Lf).mean():.4f} mu={Lf.mean():+.4f}")
    cv8 = rows["pi/8"][3]
    mu8, sd8 = rows["pi/8"][0], rows["pi/8"][1]
    stable = cv8 <= 0.25 and abs(mu8) > 3 * sd8 / np.sqrt(S)
    kills = [all(rows[k][3] > 0.50 for k in rows), not mono]
    print(f"S4b verdict: {'KILL' if any(kills) else ('PASS-lean (F4-conditional)' if stable else 'UNRESOLVED-gray')} ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
