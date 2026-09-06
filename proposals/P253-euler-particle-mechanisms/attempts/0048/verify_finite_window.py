#!/usr/bin/env python3
"""Exact algebra checks for the P253/0048 finite-window construction."""

import sympy as sp


q, x, p = sp.symbols("q x p", positive=True)
U, V, Z = sp.symbols("U V Z", positive=True)
Ux, Vx = sp.symbols("Ux Vx")
LapU, LapV, LapZ = sp.symbols("LapU LapV LapZ")

# Expand the exact local Cao equation
# -Delta w + q/(1+q*x) w_x - (1+q*x)^2 w^p = 0.
w_power = (
    U**p
    + q * p * U ** (p - 1) * V
    + q**2
    * (p * U ** (p - 1) * Z + p * (p - 1) * U ** (p - 2) * V**2 / 2)
)
transport = sp.series(q / (1 + q * x), q, 0, 3).removeO() * (
    Ux + q * Vx
)
geometry = (1 + q * x) ** 2 * w_power
residual = -(LapU + q * LapV + q**2 * LapZ) + transport - geometry
c1 = sp.expand(residual).coeff(q, 1)
c2 = sp.expand(residual).coeff(q, 2)

expected_c1 = -LapV + Ux - p * U ** (p - 1) * V - 2 * x * U**p
expected_c2 = (
    -LapZ
    + Vx
    - x * Ux
    - p * U ** (p - 1) * Z
    - p * (p - 1) * U ** (p - 2) * V**2 / 2
    - 2 * p * x * U ** (p - 1) * V
    - x**2 * U**p
)
assert sp.simplify(c1 - expected_c1) == 0
assert sp.simplify(c2 - expected_c2) == 0

# Physical m=0 normalization: h=lambda and sigma=k/sqrt(lambda).
k, lam, M = sp.symbols("k lam M", positive=True)
sigma = k / sp.sqrt(lam)
energy = M / sigma**2
energy_scale = sp.simplify(1 / sp.sqrt(energy))
assert energy_scale == k / (sp.sqrt(M) * sp.sqrt(lam))

# Resonance and spacing scales, with L=log(1/delta).
delta, L, C = sp.symbols("delta L C", positive=True)
nstar = 1 / (delta * L)
kelvin = delta / sp.Symbol("n", positive=True)
n = list(kelvin.free_symbols - {delta})[0]
spacing = delta / n**2
coupling = delta**2 / n
assert sp.simplify(coupling.subs(n, nstar)) == delta**3 * L
assert sp.simplify(spacing.subs(n, nstar)) == delta**3 * L**2
assert sp.simplify((coupling / spacing).subs(n, nstar)) == 1 / L

# Spectral window / crossing slope gives one bad interval O(delta^2).
crossing_slope = delta * L
bad_interval = sp.simplify(coupling.subs(n, nstar) / crossing_slope)
crossing_count = 1 / (delta * L)
assert bad_interval == delta**2
assert sp.simplify(bad_interval * crossing_count) == delta / L

# Tail estimates: below resonance is O(delta^3 log n*) and above is O(delta^3).
# The exact summands reduce to these powers before elementary series bounds.
low_summand = sp.simplify((delta**2 / n) ** 2 / (delta / n))
high_summand = sp.simplify((delta**2 / n) ** 2 / (delta**2 * L))
assert low_summand == delta**3 / n
assert high_summand == delta**2 / (L * n**2)
assert sp.simplify((delta**2 / L) / nstar) == delta**3

# The nonlinear transport critical harmonic is beyond every fixed thin-ring
# algebraic order but finite for each fixed nonzero delta.
jcrit = 1 / (delta**2 * L)
pattern_speed = delta**2 * L
assert sp.simplify(jcrit * pattern_speed) == 1
fixed_order = sp.symbols("N", positive=True)
assert sp.limit(fixed_order * pattern_speed, delta, 0, dir="+") == 0

# Same-sign KKS oscillators have an avoided crossing, not a Hamiltonian-Hopf.
tau, sig, h = sp.symbols("tau sig h", positive=True)
disc_same = sp.expand(((tau - sig) / 2) ** 2 + h**2)
disc_opposite = sp.expand(((tau - sig) / 2) ** 2 - h**2)
assert sp.simplify(disc_same - disc_opposite) == 2 * h**2
assert disc_same.subs(sig, tau) == h**2
assert disc_opposite.subs(sig, tau) == -h**2

checks = [
    "cao_first_cell",
    "cao_second_cell",
    "physical_energy_scale",
    "resonant_coupling_scale",
    "resonant_spacing_scale",
    "window_to_spacing_ratio",
    "single_crossing_bad_width",
    "dyadic_bad_measure",
    "low_tail_summand",
    "high_tail_summand",
    "high_tail_bound_scale",
    "critical_harmonic_scale",
    "fixed_order_noncritical_limit",
    "same_krein_avoided_crossing",
    "opposite_krein_exposing_sign",
]
for name in checks:
    print(f"PASS {name}")
print(f"PASS total={len(checks)}")
