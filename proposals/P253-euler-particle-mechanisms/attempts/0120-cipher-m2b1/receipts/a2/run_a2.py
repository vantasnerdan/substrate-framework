"""A2 archive: Rankine-core rigidity (repaired attempt). Run: python3 run_a2.py. CPython 3.12.2, numpy.

Banked intermediate failures: (1) naive material loop in exterior-model field -> ell 4.2/16.3
(void inside r<a: model has no core interior); (2) explicit-Euler rotation blowup 8e13x;
(3) 3-vector broadcast bug. Repair: Rankine solid rotation (exact matrix split) + other-ring
strain only. Pyright/lint static notices only; runtime exit 0.
"""
import numpy as np
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', '0108-cipher-radical', 'receipts', 'poc2-filament'))
from run_poc2 import rk4, mutual, self_V

G, a, T, dt = 1.0, 0.05, 4.088, 0.004
Om = G / (2 * np.pi * a * a)
c0, s0 = np.cos(Om * dt), np.sin(Om * dt)
Rot = np.array([[c0, -s0], [s0, c0]])
s = np.array([0.773723, -1e-9, 1.185226, 1e-9])
th = np.linspace(0, 2 * np.pi, 64, endpoint=False)
C1 = np.stack([np.full(64, s[0]) + a * np.cos(th), np.full(64, s[1]) + a * np.sin(th)], 1)
C2 = np.stack([np.full(64, s[2]) + a * np.cos(th), np.full(64, s[3]) + a * np.sin(th)], 1)


def M2(p, Ro, Zo):
    m = mutual(p[0], p[1], Ro, Zo)
    return np.array([m[0], m[2]])


def step(C, Rc, Zc, Ro, Zo, Vring):
    cen = np.array([Rc, Zc])
    rel = (C - cen) @ Rot.T
    m = np.array([M2(p, Ro, Zo) for p in C])
    base = M2(cen, Ro, Zo)
    return cen + rel + dt * (Vring + (m - base))


def ell(C):
    c = C.mean(0)
    r = np.linalg.norm(C - c, axis=1)
    return r.max() / r.min(), r.mean() / a


for _ in range(int(T / dt)):
    s = rk4(s, dt)
    V1 = np.array([M2(np.array([s[0], s[1]]), s[2], s[3])[0],
                   self_V(s[0]) + M2(np.array([s[0], s[1]]), s[2], s[3])[1]])
    V2 = np.array([M2(np.array([s[2], s[3]]), s[0], s[1])[0],
                   self_V(s[2]) + M2(np.array([s[2], s[3]]), s[0], s[1])[1]])
    C1 = step(C1, s[0], s[1], s[2], s[3], V1)
    C2 = step(C2, s[2], s[3], s[0], s[1], V2)
print("t=T ell=", [round(x, 4) for x in ell(C1)], [round(x, 4) for x in ell(C2)])
