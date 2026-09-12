#!/usr/bin/env python3
"""Defect dynamics intermediates (beacon, under 0157, FIRED-ON-#87).
Fence: demonstrate, never declare alive. Envelope: post-f2-build-design.
Units: a=1, b0=1, mu=1 (dipole-run convention); rho, tau stated imports.

D1 stability: E''(d) analytic + numeric -> pair configurationally
  UNSTABLE (annihilation channel); single defect topologically
  protected (I-Sing charge conservation) -> persistence moral banked.
  (standard dislocation effective mass, model-level import, stated) ->
  analytic t_c (erf-form, core-contact stop) vs RK4 integration.
D3 drift: imposed shear -> Peach-Koehler glide force, overdamped
  estimate v ~ F.tau/m* (tau = F-A tangle relaxation, stated).
Exit nonzero on any internal contradiction.
"""

import math

import numpy as np

MU = 1.0
B0 = 1.0
A = 1.0
C = MU * B0 * B0 / (2 * np.pi)
RHO = 1.0   # stated import: background mass density, model-level
TAU = 10.0  # stated import: F-A tangle relaxation time, model-level


def E_dip(d):
    return C * np.log(d / A)


def main() -> None:
    d0 = 10.0
    # D1: second variation, analytic + numeric second difference
    e2ana = -C / d0 ** 2
    h = 1e-3 * d0
    e2num = (E_dip(d0 + h) - 2 * E_dip(d0) + E_dip(d0 - h)) / h ** 2
    print(f"D1: E''(d0) analytic={e2ana:.6f} numeric={e2num:.6f} (<0: pair "
          "unstable to collapse = annihilation channel)")
    assert e2ana < 0 and abs(e2num - e2ana) / abs(e2ana) < 1e-4
    print("D1: single defect persists (charge conserved, I-Sing 0147); "
          "only pairs annihilate (particle-like persistence moral)")
    # D2: collapse trajectory, m* frozen at d0 (stated approx)
    ms = RHO * B0 * B0 * np.log(d0 / A) / 2
    # t_c to d=2a (not 0): sqrt(pi).erf(sqrt(ln(d0/2a))) — stopping at
    # core contact, matching the integrator
    tc_ana = d0 * np.sqrt(np.pi * ms / (2 * C)) * math.erf(
        np.sqrt(np.log(d0 / (2 * A))))
    def acc(x):
        return -(C / (ms * x))
    dt, d, v, t = 1e-3, d0, 0.0, 0.0
    while d > 2 * A and t < 10 * tc_ana:  # RK4 (explicit Euler overshoots)
        k1v, k1x = acc(d), v
        k2v, k2x = acc(d + k1x * dt / 2), v + k1v * dt / 2
        k3v, k3x = acc(d + k2x * dt / 2), v + k2v * dt / 2
        k4v, k4x = acc(d + k3x * dt), v + k3v * dt
        v += (k1v + 2 * k2v + 2 * k3v + k4v) * dt / 6
        d += (k1x + 2 * k2x + 2 * k3x + k4x) * dt / 6
        t += dt
    print(f"D2: collapse d0={d0}->2a: analytic t_c={tc_ana:.4f}, "
          f"integrated t={t:.4f} (m*={ms:.4f} stated-import)")
    assert d <= 2 * A, "D2 contradiction: pair did not collapse"
    assert abs(t - tc_ana) / tc_ana < 0.05, "D2 contradiction: t_c mismatch"
    F = B0 * 0.01 * MU  # Peach-Koehler glide force, imposed shear 0.01mu
    v_drift = F * TAU / ms
    print(f"D3: glide force F={F:.5f} -> drift v~{v_drift:.5f} "
          "(overdamped estimate, tau-import tiered; propagation direction "
          "for SYN missing-2)")
    assert v_drift > 0
    print("DYNAMICS verdict: pair-unstable + collapse + drift demonstrated; "
          "single-charge persistence banked; ALIVE not declared (fence)")


if __name__ == "__main__":
    main()
