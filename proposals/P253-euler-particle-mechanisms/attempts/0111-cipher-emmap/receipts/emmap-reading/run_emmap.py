"""EM-map archive: reading selection (current-loop vs flux-tube) + motional E read-off.

Reproduces receipts/emmap-reading verdicts (eval "EM-map reading selection", 2026-09-10).
Run: python3 run_emmap.py. Env: CPython 3.12.2, numpy 1.26.4 (np.trapz; no trapezoid).
Scope: one-way kinematic map only (R-EM1). No back-reaction (R-EM2 open), no charge value (R-EM3).
"""
import time
import numpy as np

R, Gam, V, a, b, Phi = 1.0, 1.0, 0.3, 0.05, 0.1, 1.0


def B_loop(sig, z, nq=720):
    ph = np.linspace(0, 2 * np.pi, nq, endpoint=False)
    dph = 2 * np.pi / nq
    sx, sy = R * np.cos(ph), R * np.sin(ph)
    rx, ry, rz = sig - sx, -sy, z + 0 * sig
    r = np.sqrt(rx * rx + ry * ry + rz * rz + a * a)
    dlx, dly = -R * np.sin(ph) * dph, R * np.cos(ph) * dph
    f = Gam / (4 * np.pi * r ** 3)
    return np.sum(f * dly * rz), np.sum(f * (dlx * ry - dly * rx))


def main():
    t0 = time.time()
    for rt in [0.15, 0.2, 0.3]:
        n = 60
        u = np.linspace(-rt, rt, n)
        du = u[1] - u[0]
        inside = u[:, None] ** 2 + u[None, :] ** 2 < b * b
        flux_tube = np.sum(inside) * (Phi / (np.pi * b * b)) * du * du
        print(f"test r={rt}: flux_tube/Phi={flux_tube / Phi:.5f} | current-loop toroidal flux=0")
    print(f"self-disk inductance L=R[ln(8R/a)-2]={np.log(8 * R / a) - 2:.5f} (geometry, not charge)")
    for s, z in [(1.5, 0.0), (0.0, 1.5), (2.0, 1.0)]:
        Br, Bzz = B_loop(s, z)
        print(f"(sig={s},z={z}): B=({Br:.4f},{Bzz:.4f}) E_phi_motional={-V * Br:.4f}")
    print(f"elapsed {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
