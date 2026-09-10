#!/usr/bin/env python3
"""P253/0114 S9 exposing probe (beacon, attempt-local, frozen design.md).

Advects passive tracers in the FROZEN exact Hill comoving field (S9 (2.5),
cross-checked in design.md) and evaluates the observable pair (D, F)
against frozen criteria. Scope: observable sensitivity only.
"""

from __future__ import annotations

import numpy as np

U = 2.0 / 15.0  # W_H, S9 section 2.2


def hill_velocity(r, z):
    r = np.asarray(r, dtype=float)
    z = np.asarray(z, dtype=float)
    rho2 = r**2 + z**2
    rho = np.sqrt(np.maximum(rho2, 1e-300))
    inside = rho2 <= 1.0
    ur = np.empty_like(rho)
    uz = np.empty_like(rho)
    ur[inside] = 1.5 * U * r[inside] * z[inside]
    uz[inside] = 1.5 * U * (1.0 - z[inside]**2 - 2.0 * r[inside]**2)
    out = ~inside
    rho5 = rho[out] ** 5
    ur[out] = 1.5 * U * r[out] * z[out] / rho5
    uz[out] = -U * (1.0 - 1.0 / rho[out]**3 + 1.5 * r[out]**2 / rho5)
    return ur, uz


def run(n: int, dt: float, T: float, delta: float = 0.05, seed: int = 0):
    rng = np.random.default_rng(seed)
    # behind-hemisphere shell set + full-shell control set
    th = rng.uniform(np.pi / 2, 3 * np.pi / 2, n)  # polar angle from +z... use spherical
    rr = rng.uniform(1.0, 1.0 + delta, n)
    # spherical: r = rr sin(th), z = rr cos(th); behind = z<0 -> th in (pi/2, 3pi/2)
    r = rr * np.sin(th)
    z = rr * np.cos(th)
    r = np.abs(r) + 1e-9
    thc = rng.uniform(0.0, 2 * np.pi, n)
    rrc = rng.uniform(1.0, 1.0 + delta, n)
    rc = np.abs(rrc * np.sin(thc)) + 1e-9
    zc = rrc * np.cos(thc)
    steps = int(T / dt)
    D, F = np.empty(steps + 1), np.empty(steps + 1)
    D[0] = float(np.max(np.sqrt(r**2 + z**2)))
    F[0] = 1.0
    for k in range(steps):
        r, z = rk4_step(r, z, dt)
        rc, zc = rk4_step(rc, zc, dt)
        D[k + 1] = float(np.max(np.sqrt(r**2 + z**2)))
        F[k + 1] = float(np.mean(np.sqrt(r**2 + z**2) < 2.0))
    t = np.arange(steps + 1) * dt
    return t, D, F

def rk4_step(r, z, dt):
    k1 = hill_velocity(r, z)
    k2 = hill_velocity(r + 0.5 * dt * k1[0], z + 0.5 * dt * k1[1])
    k3 = hill_velocity(r + 0.5 * dt * k2[0], z + 0.5 * dt * k2[1])
    k4 = hill_velocity(r + dt * k3[0], z + dt * k3[1])
    r2 = np.abs(r + dt * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6) + 1e-12
    z2 = z + dt * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6
    return r2, z2


def verdict(t, D, F):
    m = (t >= 10.0) & (t <= 30.0)
    A = np.vstack([t[m], np.ones_like(t[m])]).T
    s, d0 = np.linalg.lstsq(A, D[m], rcond=None)[0]
    pred = s * t[m] + d0
    ss = 1 - np.sum((D[m] - pred) ** 2) / np.sum((D[m] - D[m].mean()) ** 2)
    return s, ss, F[-1], bool(s > 0.05 and ss > 0.95 and F[-1] > 0.9)


def main() -> None:
    import json
    out = {}
    for tag, n, dt in (("base", 4000, 0.02), ("N2", 8000, 0.02), ("dt2", 4000, 0.01)):
        t, D, F = run(n, dt, 30.0)
        s, ss, f30, ok = verdict(t, D, F)
        print(f"{tag}: slope={s:.4f} R2={ss:.4f} F30={f30:.4f} -> "
              f"{'EXPOSED' if ok else 'BLIND'}")
        out[tag] = {"slope": s, "R2": ss, "F30": f30,
                    "verdict": "EXPOSED" if ok else "BLIND"}
    with open("proposals/P253-euler-particle-mechanisms/attempts/0114-beacon-s9/"
              "probe-result.json", "w") as fh:
        json.dump(out, fh, indent=2)
    base, n2, d2 = out["base"], out["N2"], out["dt2"]
    conv = (abs(n2["slope"] - base["slope"]) / base["slope"] < 0.2
            and abs(d2["slope"] - base["slope"]) / base["slope"] < 0.2)
    print("convergence:", "PASS" if conv else "FAIL")
    if not (base["verdict"] == "EXPOSED" and conv):
        raise SystemExit("probe did not meet frozen EXPOSED criteria (recorded as-is)")


if __name__ == "__main__":
    main()
