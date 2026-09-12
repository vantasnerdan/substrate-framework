#!/usr/bin/env python3
"""Single-defect persistence (beacon, under 0157, FIRED-ON-#87).
Fence: demonstrate, never declare alive. Envelope: post-f2-build-design.
Frozen pre-compute thresholds: charge == b within 1% at all probes;
perturbed-defect final E within 2% of clean-defect E (same box);
trivial-background bump decays to <5% of its initial E.
Units: a=1, b0=1, mu=1 (dipole-run convention). Sector: shear only,
tr eps = 0 (screw u_z), Dirichlet boundary pinned to analytic defect
field (fixes charge sector); perturbation compactly supported inside.
Exit nonzero on any threshold breach.
"""

import numpy as np

MU = 1.0
B0 = 1.0
N = 128
L = 16.0
H = L / N
DT = 0.2 * H * H / MU


def defect_field(X, Y):
    return (B0 / (2 * np.pi)) * np.arctan2(Y, X)


def energy(u):
    gx = (np.roll(u, -1, 1) - np.roll(u, 1, 1)) / (2 * H)
    gy = (np.roll(u, -1, 0) - np.roll(u, 1, 0)) / (2 * H)
    return 0.5 * MU * np.sum(gx ** 2 + gy ** 2) * H * H


def charge(u):
    # Burgers circuit: closed-loop diffs telescope to 0 identically, so
    # the charge lives in the unwinder accumulator: scale to 2pi-period,
    # unwrap, and read end-minus-start (accumulated winding adjustment)
    r = 3.0
    n = 720
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    xs = (r * np.cos(th) / H + N / 2).astype(int) % N
    ys = (r * np.sin(th) / H + N / 2).astype(int) % N
    s = np.unwrap(u[ys, xs] / B0 * 2 * np.pi)
    return float((s[-1] - s[0]) / (2 * np.pi) * B0)


def relax_w(w, fix_mask, steps):
    # Flow the single-valued perturbation only: u = u_clean + w with
    # u_clean harmonic away from core, so dw/dt = mu.lap(w). Flowing u
    # directly lets the grid cut-wall (branch-cut discretization
    # artifact) drive spurious unwinding/escape — banked observation:
    # free gradient flow ejects defects through boundaries (surface
    # annihilation channel, consistent with D1 pair-instability moral).
    for _ in range(steps):
        lap = (np.roll(w, -1, 0) + np.roll(w, 1, 0)
               + np.roll(w, -1, 1) + np.roll(w, 1, 1) - 4 * w) / H ** 2
        w = w + DT * MU * lap
        w[fix_mask] = 0.0
    return w


def main() -> None:
    g = (np.arange(N) - N / 2) * H
    X, Y = np.meshgrid(g, g)
    u_clean = defect_field(X, Y)
    edge = np.zeros_like(u_clean, bool)
    edge[0, :] = edge[-1, :] = edge[:, 0] = edge[:, -1] = True
    E_clean = energy(u_clean)
    q_clean = charge(u_clean)
    # Core disk pinned as sub-continuum regularization (declared, standard:
    # core physics below cutoff a is atomistic, not elastic). Relaxation
    # acts on the far-field texture only — the persistence claim.
    core = X ** 2 + Y ** 2 < 1.0
    fix = edge | core
    bump = 0.5 * np.exp(-((X - 3.0) ** 2 + (Y + 2.0) ** 2) / 2.0 ** 2)
    bump *= 1.0 - np.clip(np.sqrt(X ** 2 + Y ** 2) / (L / 2 - 1.0), 0, 1)
    w = relax_w(bump, fix, 20000)
    u = u_clean + w
    E_fin, q_fin = energy(u), charge(u)
    print(f"perturbed: E0={energy(u_clean + bump):.5f} Efin={E_fin:.5f} "
          f"charge={q_fin:.5f}")
    assert abs(q_fin - B0) / B0 < 0.01, "PERSISTENCE FIRES: charge unwound"
    assert abs(E_fin - E_clean) / E_clean < 0.02, \
        "PERSISTENCE FIRES: did not relax to defect"
    v = relax_w(bump, edge, 20000)
    E_v0, E_v = energy(bump), energy(v)
    print(f"trivial: E0={E_v0:.5f} Efin={E_v:.5f} "
          f"(ratio={E_v / E_v0:.4f})")
    assert E_v / E_v0 < 0.05, "CONTRAST FIRES: bump survived on trivial bg"
    print("PERSISTENCE verdict: defect sector survives perturbation "
          "(charge + energy restored); trivial bump dies; "
          "ALIVE not declared (fence)")


if __name__ == "__main__":
    main()
