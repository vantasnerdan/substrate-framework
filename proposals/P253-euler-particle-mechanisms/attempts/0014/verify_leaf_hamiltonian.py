#!/usr/bin/env python3
"""Exact symbolic checks for the 0014 leaf/Hamiltonian calculation."""

import sympy as sp

r, e, xi = sp.symbols("r e xi", nonzero=True, real=True)
b, kz, az, lam = sp.symbols("b kz az lam", real=True)
Aprime, Fprime, apsi, D = sp.symbols(
    "Aprime Fprime apsi D", nonzero=True, real=True
)

# Cao active-core differentiation: zeta=const+xi/(e*r**2).
az_cao = b / (e * r**2) - 2 * xi * kz / (e * r**4)
J0 = sp.Matrix([[az, b], [b, 0]])
L0 = sp.Matrix([[0, -e], [-e, r**-2]])
cao_product = sp.simplify((J0 * L0).subs(az, az_cao))
cao_expected = sp.Matrix([[-e * b, 2 * xi * kz / r**4], [0, -e * b]])
assert sp.simplify(cao_product - cao_expected) == sp.zeros(2)

# Double convective characteristic and nonzero nilpotent square.
charpoly = sp.factor((sp.I * cao_product - lam * sp.eye(2)).det())
assert sp.simplify(charpoly - (lam + sp.I * e * b) ** 2) == 0
nilpotent = sp.simplify(cao_product + e * b * sp.eye(2))
assert nilpotent**2 == sp.zeros(2)

# General regular-swirl patch: A'=-1/F', b=F' a_psi.
general_product = sp.Matrix([[az, b], [b, 0]]) * sp.Matrix(
    [[0, Aprime], [Aprime, D]]
)
general_diag = sp.simplify(
    general_product[0, 0].subs({Aprime: -1 / Fprime, b: Fprime * apsi})
)
assert general_diag == -apsi
assert general_product[1, 0] == 0

# WKB Cao local cross term: eta=chi/(e*r**2) gives the negative of swirl norm.
chi = sp.symbols("chi", real=True)
eta_lead = chi / (e * r**2)
local_q = sp.simplify(chi**2 / r**2 - 2 * e * eta_lead * chi)
assert local_q == -chi**2 / r**2

# Exposing mutation: the wrong cross sign would be positive, not the asserted result.
wrong_local_q = sp.simplify(chi**2 / r**2 + 2 * e * eta_lead * chi)
assert wrong_local_q != -chi**2 / r**2
assert wrong_local_q == 3 * chi**2 / r**2

print("cao_product=", cao_product)
print("cao_charpoly=", charpoly)
print("cao_nilpotent_square=", nilpotent**2)
print("general_convective_diagonal=", general_diag)
print("wkb_local_q=", local_q)
print("wrong_cross_sign_local_q=", wrong_local_q)
print("checks=7 passed")
