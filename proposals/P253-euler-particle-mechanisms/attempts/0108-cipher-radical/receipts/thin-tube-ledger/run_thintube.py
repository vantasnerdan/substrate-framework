"""Thin-tube ledger archive: scale/label separation bookkeeping (M1 input).

Reproduces receipts/thin-tube-ledger eval "Thin-tube ledger PoC" (2026-09-10).
Run: python3 run_thintube.py. Env: CPython 3.12.2, numpy 1.26.4.
Scope: analytic thin-ring formulas + helicity-link integers; exploratory only.
"""
import numpy as np

rho, Gamma, R, a = 1.0, 1.0, 1.0, 0.05
L = np.log(8 * R / a)
I = np.pi * rho * Gamma * R ** 2
E = 0.5 * rho * Gamma ** 2 * R * (L - 1.75)
V = Gamma / (4 * np.pi * R) * (L - 0.25)
m_star = I / V
G1 = G2 = 1.0
for n in [0, 1, 2]:
    print(f"linking n={n}: H_link={2 * n * G1 * G2:.1f} (discrete jump 2.0), "
          f"I_pair={2 * I:.4f} (continuous in R,Gamma)")
print(f"single ring: I={I:.4f} E={E:.4f} V={V:.4f} m*=I/V={m_star:.4f} "
      f"(Gamma-free, scales rho*R^3*log-ratio)")
print(f"Hill ball added mass: m={(2 / 3) * np.pi * rho * R ** 3:.4f}")
print("ledger: R,Gamma continuous knobs; H-linking integer; m from I-V not from Gamma declaration")
