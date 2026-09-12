#!/usr/bin/env python3
"""Exact algebra checks for the P253/0083 scale and hierarchy ledger."""

import sympy as sp


p, R, kappa, Lambda = sp.symbols("p R kappa Lambda", positive=True)
epsilon, d = sp.symbols("epsilon d", positive=True)

C_s = (Lambda / kappa) ** ((p - 1) / 2) * R ** (-(p + 1) / 2)
H_0 = (R * C_s) ** (-sp.Rational(2, 1) / (p - 1))

# Local PDE balance and physical circulation determine the same coefficient.
pde_balance = R**2 * C_s**2 * H_0 ** (p - 1)
mass_balance = H_0**p * R * C_s**2 * Lambda
# Positivity makes logarithms injective and avoids SymPy's deliberately
# conservative cancellation of symbolic noninteger powers.
assert sp.simplify(sp.expand_log(sp.log(pde_balance), force=True)) == 0
assert sp.simplify(sp.expand_log(sp.log(mass_balance / kappa), force=True)) == 0

# A center/auxiliary-radius error d*epsilon^2 produces only an epsilon^3
# physical-scale error through the exact second equation of Cao (3.36).
r = sp.symbols("r", positive=True)
C_of_r = (Lambda / kappa) ** ((p - 1) / 2) * r ** (-(p + 1) / 2)
center_to_scale = sp.simplify(
    sp.limit(
        (epsilon * C_of_r.subs(r, R + d * epsilon**2) - epsilon * C_s)
        / epsilon**3,
        epsilon,
        0,
    )
)
assert sp.simplify(center_to_scale + d * (p + 1) * C_s / (2 * R)) == 0

# The safe area-radius estimate is O(epsilon^2).  At epsilon~1/N this becomes
# O(1/N) after multiplication by n~N, and is o(h_N) for the exposing choice
# h_N=sqrt(1/N).
N = sp.symbols("N", positive=True)
wave_error = 1 / N
h_N = sp.sqrt(1 / N)
assert sp.limit(wave_error / h_N, N, sp.oo) == 0

# The Kelvin--Hicks row on delta=k/(N*P) separates into log N plus a fixed
# rational-ray constant.
P_ray, k_ray = sp.symbols("P_ray k_ray", positive=True)
delta_N = k_ray / (N * P_ray)
log_split = sp.expand_log(sp.log(8 / delta_N), force=True)
expected_log_split = sp.log(N) + sp.log(8 * P_ray / k_ray)
assert sp.simplify(log_split - expected_log_split) == 0

# The sufficient finite-speed ceiling saturates the exposed upper bound.
alpha, c_em, margin, C_rem, C_ray = sp.symbols(
    "alpha c_em margin C_rem C_ray", positive=True
)
N_em = sp.exp((c_em - margin - C_rem) / alpha - C_ray)
speed_upper_at_ceiling = sp.simplify(alpha * (sp.log(N_em) + C_ray) + C_rem)
assert speed_upper_at_ceiling == c_em - margin

print("PASS scale coefficient from PDE balance")
print("PASS scale coefficient from circulation balance")
print("PASS O(epsilon^2) center error gives O(epsilon^3) scale error")
print("PASS wave-number safe area-scale error is o(h_N)")
print("PASS massive Kelvin-Hicks logarithm split")
print("PASS finite subluminal ceiling algebra")
