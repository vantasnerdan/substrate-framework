#!/usr/bin/env python3
"""Exact exposing checks for the P253/0058 nonlinear range calculation.

This verifier checks finite-dimensional algebra only.  It does not claim a
Cao branch, a graph-domain inverse, or convergence of a nonlinear iteration.
"""

import sympy as sp


def vec(x, y, z):
    return sp.Matrix([x, y, z])


# A minimal exact quadratic Euler interaction.
k1 = vec(-1, 0, 0)
eta1 = vec(0, -1, 0)
k2 = vec(0, -1, 0)
eta2 = vec(0, 0, -1)
k = k1 + k2
v1 = k1.cross(eta1) / k1.dot(k1)
v2 = k2.cross(eta2) / k2.dot(k2)
numerator = sp.simplify(k.cross(v1.cross(eta2) + v2.cross(eta1)))

assert k1.dot(eta1) == 0
assert k2.dot(eta2) == 0
assert k.dot(numerator) == 0
assert v1 == vec(0, 0, 1)
assert v2 == vec(1, 0, 0)
assert numerator == vec(1, -1, 0)

# Choose a frozen relative velocity for which the output transport divisor is
# exactly zero.  The nonzero numerator proves that the quadratic tensor is not
# universally divisible by that divisor.
W = vec(1, -1, 0)
divisor = sp.expand(W.dot(k))
assert divisor == 0
assert numerator != sp.zeros(3, 1)

# Dynamically accessible principal symbol:
# i k x (xi x omega) = i[(k.omega)xi-(k.xi)omega].
kx, ky, kz, xx, xy, xz, ox, oy, oz = sp.symbols(
    "kx ky kz xx xy xz ox oy oz", real=True
)
ks = vec(kx, ky, kz)
xis = vec(xx, xy, xz)
oms = vec(ox, oy, oz)
triple = sp.simplify(ks.cross(xis.cross(oms)))
rhs = sp.simplify(xis * ks.dot(oms) - oms * ks.dot(xis))
assert sp.simplify(triple - rhs) == sp.zeros(3, 1)

# On k.xi=0 and k.omega != 0 this is a scalar multiple of xi, so the DA map
# spans the transverse vorticity plane at symbol level.  This is not a range
# theorem for the nonnormal critical-layer operator.

# Centralizer of a positive axisymmetric toroidal vorticity zeta*partial_theta.
# In one nonzero theta harmonic ell, [Y,zeta*partial_theta]=0 has the displayed
# triangular coefficient matrix.  Its determinant is nonzero.
ell, zeta, zr, zz = sp.symbols("ell zeta zr zz", nonzero=True, real=True)
I = sp.I
centralizer_matrix = sp.Matrix(
    [
        [-I * ell * zeta, 0, 0],
        [zr, -I * ell * zeta, zz],
        [0, 0, -I * ell * zeta],
    ]
)
centralizer_det = sp.factor(centralizer_matrix.det())
assert centralizer_det == I * ell**3 * zeta**3

# Full frozen linearized vorticity symbol retains the Hodge term.  On the
# transverse plane (k.eta=0), (k cross)^2 eta=-|k|^2 eta.
ex, ey, ez = sp.symbols("ex ey ez", real=True)
etas = vec(ex, ey, ez)
cross_square = sp.expand(ks.cross(ks.cross(etas)))
lagrange = sp.expand(ks * ks.dot(etas) - ks.dot(ks) * etas)
assert sp.simplify(cross_square - lagrange) == sp.zeros(3, 1)

# Supervisor-corrected physical KKS rotation sign.
a, sigma = sp.symbols("a sigma", real=True, nonzero=True)
# X_J=-a*ell*q_s and Omega(q_s,q_c)=-sigma.
Jprime = sp.simplify((-a * ell) * (-sigma))
j2 = sp.simplify(Jprime / (2 * a))
assert Jprime == a * ell * sigma
assert j2 == ell * sigma / 2

print("quadratic_input_1:", tuple(k1), tuple(eta1))
print("quadratic_input_2:", tuple(k2), tuple(eta2))
print("quadratic_output_k:", tuple(k))
print("transport_divisor:", divisor)
print("nondivisible_numerator:", tuple(numerator))
print("da_symbol_triple_product: PASS")
print("nonzero_harmonic_centralizer_det:", centralizer_det)
print("frozen_hodge_cross_square: PASS")
print("physical_kks_Jprime:", Jprime)
print("physical_kks_j2:", j2)
print("checks: 13 passed")
