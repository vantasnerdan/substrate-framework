"""B1 phase-charge lock (FROZEN F-bar in ../../06-synthesis-spec/01-b1-fbar.md):
framed-unknot edge source (linking L by construction) -> BS vector potential
-> phase around core circuit; control far circuit; N-leg 128 vs 256.
Usage: python3 run_ab.py
"""
import time

import numpy as np

N = 256
R0 = 1.0


def frame_circle(n):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    C = np.stack([R0 * np.cos(t), R0 * np.sin(t), np.zeros(n)], axis=1)
    T = np.stack([-np.sin(t), np.cos(t), np.zeros(n)], axis=1)
    Nn = np.stack([-np.cos(t), -np.sin(t), np.zeros(n)], axis=1)
    B = np.zeros((n, 3))
    B[:, 2] = 1.0
    return C, T, Nn, B, t


def edge_curve(L, n, eps=0.05):
    C, _, Nn, B, t = frame_circle(n)
    return C + eps * (np.cos(L * t)[:, None] * Nn + np.sin(L * t)[:, None] * B)


def bs_A(src, obs):
    """A(obs) = (1/4pi) sum dl x r / |r|^3 (unit current)."""
    dl = np.roll(src, -1, axis=0) - src
    mid = (src + np.roll(src, -1, axis=0)) / 2
    R = obs[:, None, :] - mid[None, :, :]
    r3 = np.linalg.norm(R, axis=2) ** 3 + 1e-30
    return (np.cross(dl[None, :, :], R) / r3[:, :, None]).sum(axis=1) / (4 * np.pi)


def phase(test, src):
    A = bs_A(src, test)
    dl = np.roll(test, -1, axis=0) - test
    return float((A * dl).sum(axis=1).sum() )


def main():
    t0 = time.time()
    C, _, _, _, _ = frame_circle(N)
    Phis = []
    for L in (1, 2, 3):
        E = edge_curve(L, N)
        Phis.append(phase(C, E))
    ctrl = C + np.array([5.0, 0, 0])
    Pc = phase(ctrl, edge_curve(1, N))
    print("Phis:", " ".join(f"{p:.4f}" for p in Phis), f"ctrl={Pc:.4f}")
    # N-leg
    C2, _, _, _, _ = frame_circle(128)
    P1 = phase(C2, edge_curve(1, 128))
    print(f"N-leg: Phi1(N=256)={Phis[0]:.4f} Phi1(N=128)={P1:.4f} dPhi={abs(Phis[0]-P1):.2e}")
    ok = all(abs(p - L) <= 0.05 * max(L, 1) for p, L in zip(Phis, (1, 2, 3))) and abs(Pc) <= 0.05
    kill = any(abs(p - L) > 0.25 for p, L in zip(Phis, (1, 2, 3))) or abs(Pc) > 0.25
    print(f"B1 verdict: {'HOLD' if ok else ('KILL' if kill else 'UNRESOLVED-gray')} ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
