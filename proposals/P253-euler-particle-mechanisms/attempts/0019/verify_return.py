#!/usr/bin/env python3
"""Exact algebra for the P253/0019 solvable Gavrilov return."""

import sympy as sp

R, t = sp.symbols("R t", positive=True, real=True)
pi = sp.pi
Omega = R
gamma = R / sp.sqrt(2)
T = 2 * pi / Omega

# Moving basis (n,t,Y): rigid meridional rotation plus helical radial shear.
Arot = sp.Matrix([[0, -Omega, 0], [Omega, 0, 0], [0, 0, 0]])
Shear = sp.Matrix([[0, 0, 0], [0, 0, 0], [gamma, 0, 0]])
F = sp.eye(3) + sp.Matrix([[0, 0, 0], [0, 0, 0], [gamma * t, 0, 0]])
assert F.det() == 1
assert sp.simplify(F.inv().T - sp.Matrix([[1, 0, -gamma * t], [0, 1, 0], [0, 0, 1]])) == sp.zeros(3)

# Full Leray geometric-optics amplitude generator for the returning k=e_t ray.
k = sp.Matrix([0, 1, 0])
P_k = k * k.T
amplitude_generator = sp.simplify(-2 * Arot - Shear + 2 * P_k * (Arot + Shear))
expected_generator = sp.Matrix([[0, 2 * Omega, 0], [0, 0, 0], [-gamma, 0, 0]])
assert amplitude_generator == expected_generator

# The divergence-free subspace A_t=0 contains the exact normal-to-swirl shear.
M = sp.Matrix([[1, 0], [-sp.simplify(gamma * T), 1]])
M_expected = sp.Matrix([[1, 0], [-sp.sqrt(2) * pi, 1]])
assert sp.simplify(M - M_expected) == sp.zeros(2)
J2 = sp.Matrix([[0, 1], [-1, 0]])
assert sp.simplify(M.T * J2 * M - J2) == sp.zeros(2)
assert M.det() == 1
assert M.eigenvals() == {sp.Integer(1): 2}
assert (M - sp.eye(2)).rank() == 1

witness = sp.Matrix([1, 0])
gain = sp.simplify((M * witness).dot(M * witness) / witness.dot(witness))
assert gain == 1 + 2 * pi**2

# Energy-weighted K/shear pairing retained by the 0020 correction.
xi, kz, bxi, r, kappa = sp.symbols("xi kz bxi r kappa", nonzero=True, real=True)
upper = sp.I * 2 * xi * kz / (r**4 * kappa)
lower = sp.I * r**2 * bxi / kappa
pair_product = sp.simplify(upper * lower)
assert pair_product == -2 * xi * bxi * kz / (r**2 * kappa**2)

# Canonical action sign: -Omega(q,qdot)/2-Q/2 gives qdot=J grad Q.
L11, L12, L22 = sp.symbols("L11 L12 L22", real=True)
L = sp.Matrix([[L11, L12], [L12, L22]])
# Euler--Lagrange gives (-J) qdot-Lq=0 for the displayed minus sign.
assert sp.simplify((-J2).inv() * L - J2 * L) == sp.zeros(2)
# Flipping only the kinetic sign instead gives J qdot-Lq=0, hence -J Lq.
assert sp.simplify((J2).inv() * L + J2 * L) == sp.zeros(2)

print("deformation_det=", F.det())
print("cotangent_return=", sp.simplify(F.subs(t, T).inv().T))
print("amplitude_generator=", amplitude_generator)
print("axisymmetric_monodromy=", M)
print("monodromy_symplectic_residual=", sp.simplify(M.T * J2 * M - J2))
print("monodromy_jordan_rank=", (M - sp.eye(2)).rank())
print("normal_witness_energy_gain=", gain)
print("weighted_pressure_shear_product=", pair_product)
print("checks=12 passed")
