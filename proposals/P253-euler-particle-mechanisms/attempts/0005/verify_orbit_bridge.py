#!/usr/bin/env python3
"""Exact algebra checks for the P253/0005 rotation-orbit bridge."""

import sympy as sp

r, theta, phi, rho, I1, I2, j, hbar = sp.symbols(
    "r theta phi rho I1 I2 j hbar", positive=True
)

# Angular factors for |n x y|^2 and y x (n x y).
angular_energy = sp.integrate(
    sp.sin(theta) ** 3, (theta, 0, sp.pi), (phi, 0, 2 * sp.pi)
)
angular_Lz = sp.integrate(
    sp.sin(theta) ** 3, (theta, 0, sp.pi), (phi, 0, 2 * sp.pi)
)
assert sp.simplify(angular_energy - 8 * sp.pi / 3) == 0
assert sp.simplify(angular_Lz - 8 * sp.pi / 3) == 0

H = sp.simplify(rho * angular_energy * I2 / 2)
Lz = sp.simplify(rho * angular_Lz * I1)
assert H == 4 * sp.pi * rho * I2 / 3
assert Lz == 8 * sp.pi * rho * I1 / 3

# KKS sphere period and north/south chart jump.
period = sp.integrate(
    j * sp.sin(theta), (theta, 0, sp.pi), (phi, 0, 2 * sp.pi)
)
assert sp.simplify(period - 4 * sp.pi * j) == 0
A_N_phi = j * (1 - sp.cos(theta))
A_S_phi = -j * (1 + sp.cos(theta))
assert sp.simplify(A_N_phi - A_S_phi - 2 * j) == 0
equator_phase_exponent = sp.simplify(
    sp.integrate(A_N_phi.subs(theta, sp.pi / 2), (phi, 0, 2 * sp.pi)) / hbar
)
assert equator_phase_exponent == 2 * sp.pi * j / hbar

# U,L scaling: energy exponent U^2 L^3 and action/angular-momentum U L^4.
U, L, c = sp.symbols("U L c", positive=True)
energy_scale = rho * U**2 * L**3
time_scale = L / U
action_scale = sp.simplify(energy_scale * time_scale)
assert action_scale == rho * U * L**4
assert sp.simplify((c**2 * energy_scale) * (time_scale / c) - c * action_scale) == 0

print("PASS angular energy coefficient: 8*pi/3")
print("PASS H=(4*pi*rho/3) I2")
print("PASS j=(8*pi*rho/3) I1")
print("PASS KKS period=4*pi*j and A_N-A_S=2*j*dphi")
print("PASS equator holonomy exponent=2*pi*j/hbar")
print("PASS energy*time=rho*U*L^4 and Euler amplitude scaling is linear in action")
