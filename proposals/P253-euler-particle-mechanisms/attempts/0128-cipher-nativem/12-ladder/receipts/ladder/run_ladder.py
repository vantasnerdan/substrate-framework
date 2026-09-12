"""L-ladder build (FROZEN F-bar 00-fbar.md incl pre-run amendment):
Neumann mutual inductance M(L) between ribbon edges E+- of a framed circular
ring (R=1, half-width w), framing L twists edge offset around centerline.
Ladder Q(L) = M(L) - M(0); verdict on Q(L)/Q(1) = L (5% band) + sign + leakage.
Usage: python3 run_ladder.py
"""
import time

import numpy as np

R = 1.0
W = 0.02
N = 512


def edge(sign, L, ph):
    # centerline ring in xy-plane + offset along framing-twisted normal
    # Frenet of circle: normal n̂ = radial-in, binormal b̂ = ẑ
    a = L * ph
    off = sign * (W / 2.0)
    nx, nz = np.cos(a), np.sin(a)  # twist L times around centerline
    cx, cy = R * np.cos(ph), R * np.sin(ph)
    # offset dir = nx * (-radial) + nz * z
    px = cx - off * nx * np.cos(ph)
    py = cy - off * nx * np.sin(ph)
    pz = off * nz
    return np.stack([px, py, pz], axis=1)


def tangent(curve):
    d = np.roll(curve, -1, axis=0) - np.roll(curve, 1, axis=0)
    return d / np.linalg.norm(d, axis=1, keepdims=True)


def neumann(L):
    ph = np.linspace(0, 2 * np.pi, N, endpoint=False)
    A = edge(+1.0, L, ph)
    B = edge(-1.0, L, ph)
    dA = tangent(A)
    dB = tangent(B)
    ds = 2 * np.pi * R / N
    D = A[:, None, :] - B[None, :, :]
    r = np.sqrt((D ** 2).sum(axis=2))
    np.fill_diagonal(r, np.inf)  # exclude coincident (same-ph) pairs
    dot = dA @ dB.T
    return float((dot / r).sum()) * ds * ds / (4 * np.pi)


def main():
    t0 = time.time()
    Ms = {}
    for L in [0, 1, -1, 2, -2]:
        Ms[L] = neumann(L)
        print(f"L={L:+d} M={Ms[L]:.8f}", flush=True)
    Q = {L: Ms[L] - Ms[0] for L in Ms}
    print(f"Q: { {L: round(v, 8) for L, v in Q.items()} }")
    ok, why = True, []
    for L in [1, -1, 2, -2]:
        rat = Q[L] / Q[1] if Q[1] != 0 else float("nan")
        band = abs(rat - L) <= 0.05 * abs(L)
        print(f"L={L:+d} Q/Q1={rat:.5f} expect {L:+d} {'IN' if band else 'OUT'}")
        if not band:
            ok = False
            why.append(f"L={L}: {rat:.4f} vs {L}")
    # even-channel leakage: [Q(L)+Q(-L)]/2 should vanish vs odd scale
    for L in [1, 2]:
        leak = abs(Q[L] + Q[-L]) / (2 * abs(Q[1]))
        print(f"leak L={L}: {leak:.4f} {'OK' if leak <= 0.05 else 'LEAK'}")
        if leak > 0.05:
            ok = False
            why.append(f"leak L={L}: {leak:.3f}")
    vd = ("CONFIRMED (integer ladder; SYN P1 opens)" if ok
          else f"P1 DEAD ({'; '.join(why)})")
    print(f"L-ladder verdict: {vd} ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
