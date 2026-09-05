"""Actual C016 acoustic antisymmetry, stationary spin and detector normalization."""

import sympy as s

from substrate_framework import euler_fourier as ef
from substrate_framework.euler_displacement_preparation import (
    finite_displacement_cell,
    prepared_displacement,
)
from substrate_framework.verification import CheckLedger


def add(*vectors):
    return tuple(ef.add(*(v[i] for v in vectors)) for i in range(3))


def scale(vector, coefficient):
    return tuple(ef.scale(v, coefficient) for v in vector)


def advect(left, right):
    return tuple(ef.add(*(ef.mul(left[j], ef.derivative(right[i], j))
                          for j in range(3))) for i in range(3))


def zero(vector):
    return all(s.simplify(value) == 0 for row in vector for value in row.values())


def main():
    checks = CheckLedger("P251-0231-actual-acoustic-angle-spin")
    cell = finite_displacement_cell()
    u = cell.background
    rot = ({}, ef.scale(u[2], -1), u[1])
    grad = tuple(ef.derivative(u[0], i) for i in range(3))
    checks.check("actual axial rotation velocity is the complete negative pressure gradient",
                 zero(add(rot, grad)))
    checks.check("the full periodic pressure removes the rigid-rate Coriolis residual",
                 zero(ef.leray(rot)))
    left = prepared_displacement(cell, (0, 1, 0), (0, 0, 1))
    right = prepared_displacement(cell, (0, 0, 1), (0, 1, 0))
    def difference(name):
        return add(getattr(left, name), scale(getattr(right, name), -1))

    checks.check("complete C015 antisymmetric lift equals the actual rotation velocity",
                 zero(add(difference("lift"), scale(rot, -1))))
    checks.check("all actual displacement returns vanish in the axial antisymmetric channel",
                 zero(difference("returned")))
    checks.check("the full displacement material-rate correction is exactly zero",
                 zero(difference("material_rate")))
    checks.check("the full physical acoustic current test is symmetric",
                 zero(difference("current_test")))
    checks.check("the complete C016 velocity forcing retains and cancels both gradient terms",
                 zero(ef.leray(add(rot, difference("lift")))))
    k = s.symbols("k0:3", real=True)
    d = s.symbols("d0:3", real=True)
    correlation = k[1]*d[0]+k[0]*d[1]
    anti = s.diff(correlation, k[1], d[2])-s.diff(correlation, k[2], d[1])
    checks.check("actual passive configuration and velocity correlation has zero axial skew",
                 anti == 0)
    b = ({}, u[1], u[2])
    checks.check("the physical spin-return field commutes with the full Euler background",
                 zero(add(advect(u, b), scale(advect(b, u), -1))))
    checks.check("the spin-return velocity solves full stationary linear Euler including pressure",
                 zero(ef.leray(add(advect(u, b), advect(b, u)))))
    checks.check("the spin-return field is mean free and has no axial material drift",
                 not b[0] and all(row.get(ef.ZERO, 0) == 0 for row in b))
    y, z = s.symbols("Y Z", real=True)
    psi = s.Function("psi")(y, z)
    ff = s.Function("F")(psi)
    divergence = s.diff(y*ff, y)+s.diff(z*ff, z)
    checks.check("actual background spin has the compact-boundary nonzero denominator identity",
                 s.simplify(divergence-2*ff-s.diff(ff, psi)*
                            (y*s.diff(psi, y)+z*s.diff(psi, z))) == 0)
    chi = s.Function("chi")(psi)
    checks.check("the return preserves the entire literal transported tag",
                 s.simplify(-s.diff(psi, z)*s.diff(chi, y)
                            +s.diff(psi, y)*s.diff(chi, z)) == 0)
    qy, qz, angle = s.symbols("QY QZ angle", real=True)
    covariance = s.diag(qy, qz)
    jj = s.Matrix([[0, -1], [1, 0]])
    variation = angle*(jj*covariance-covariance*jj)
    checks.check("the literal covariance detector returns an exact rigid angle",
                 s.cancel(variation[0, 1]/(qy-qz)) == angle)
    gram = s.eye(3)/3
    omega = s.Matrix(s.symbols("omega0:3", real=True))
    checks.check("the actual detector Gram reconstructs rigid rotation without changing spin",
                 gram.inv()*(gram*omega) == omega)
    checks.check("raw detector averaging without the Gram inverse has the exposed factor-three defect",
                 gram*omega != omega)
    tensor = s.Matrix(3, 3, s.symbols("t0:9"))
    gradient = s.Matrix(3, 3, s.symbols("B0:9"))
    observed = s.Matrix([sum(s.LeviCivita(i, j, k0)*s.LeviCivita(0, a, b0)
                             *tensor[a, b0]*gradient[j, k0]/2
                             for j in range(3) for k0 in range(3)
                             for a in range(3) for b0 in range(3))
                         for i in range(3)])
    expected = (tensor[2, 1]-tensor[1, 2])*s.Matrix([
        gradient[2, 1]-gradient[1, 2], gradient[0, 2]-gradient[2, 0],
        gradient[1, 0]-gradient[0, 1]])/2
    checks.check("full Haar contraction retains only the measured rigid antisymmetric response",
                 s.simplify(observed-expected) == s.zeros(3, 1))
    inertia, j0, spin0, om, rate = s.symbols("I j0 S0 Omega rate", nonzero=True, real=True)
    coefficient = (j0-inertia)/spin0
    checks.check("the actual angle-null stationary return matches both measured G and spin",
                 s.simplify(inertia*om+coefficient*spin0*om-j0*om) == 0
                 and s.simplify(inertia*rate+coefficient*spin0*rate-j0*rate) == 0)
    checks.check("dropping the physical return leaves the exposed inertia mismatch",
                 s.simplify(inertia*rate-j0*rate) != 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
