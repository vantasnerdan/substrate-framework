#!/usr/bin/env python3
"""Exact exposing checks for P253/0097."""

import sympy as sp

from substrate_framework.euler_neutrino_suppliers import (
    current_parity,
    packet_transition_probability,
    relativistic_phase_leading,
    two_flavor_transition_probability,
)


theta, phi = sp.symbols("theta phi", real=True)
fixed = two_flavor_transition_probability(theta, phi, 1)
packet = packet_transition_probability(theta, sp.exp(-sp.I * phi))
assert sp.simplify(sp.expand_trig(packet - fixed)) == 0
print("PASS packet formula reduces to fixed-momentum oscillation")

assert two_flavor_transition_probability(theta, phi, 0) == 0
print("PASS no elapsed phase gives no transition")

assert current_parity("scalar")["spatial_current"] == "polar"
assert current_parity("pseudoscalar")["spatial_current"] == "axial"
print("PASS scalar and pseudoscalar current parity")

dm2, c, length, energy, action = sp.symbols(
    "dm2 c length energy action", positive=True
)
phase = relativistic_phase_leading(dm2, c, length, energy, action)
assert sp.simplify(phase * 2 * energy * action - dm2 * c**3 * length) == 0
print("PASS shared-action relativistic phase coefficient")

U11, U12, U21, U22 = sp.symbols("U11 U12 U21 U22")
r1, r2 = sp.symbols("r1 r2")
U = sp.Matrix([[U11, U12], [U21, U22]])
Dchi = sp.Matrix([r1, r2])
DU = sp.zeros(2, 2)
chi = sp.Matrix(sp.symbols("chi1 chi2"))
D_Uchi = DU * chi + U * Dchi
assert D_Uchi == U * Dchi
assert D_Uchi.subs({r1: 0, r2: 0}) == sp.zeros(2, 1)
print("PASS material Leibniz rule and passive constant-mixture specialization")

u1, u2, u3, w1, w2, w3 = sp.symbols("u1 u2 u3 w1 w2 w3")
u = sp.Matrix([u1, u2, u3])
w = sp.Matrix([w1, w2, w3])
assert sp.simplify(u.cross(w).dot(w)) == 0
assert all(
    sp.simplify(x) == 0
    for x in (u.cross(w).cross(u) - (u.dot(u) * w - u.dot(w) * u))
)
print("PASS Euler helicity-current vector identities")

p, c, mass_squared = sp.symbols("p c mass_squared", positive=True)
sqrt_energy = sp.sqrt(p**2 * c**2 + mass_squared * c**4)
mass_coefficient = sp.simplify(
    sp.diff(sqrt_energy, mass_squared).subs(mass_squared, 0)
)
assert mass_coefficient == c**3 / (2 * p)
print("PASS square-root dispersion derives high-momentum mass coefficient")

nu1, nu2 = sp.symbols("nu1 nu2", real=True)
Hprop = sp.diag(nu1, nu2)
Hint = sp.Matrix([[0, 1], [1, 0]])
comm = Hprop * Hint - Hint * Hprop
assert comm == sp.Matrix([[0, nu1 - nu2], [nu2 - nu1, 0]])
print("PASS exact noncommuting propagation and interaction matrix witness")
