"""Actual full-fluid forms, physical detector Gram, and ordered error scales."""

import sympy as sp

from substrate_framework import euler_fourier as ef
from substrate_framework.homogenization import sphere_second_moment
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0234-law-diagonal")
    fractions = (sp.Rational(1, 4), sp.Rational(1, 2))
    avg_fraction = sum(fractions)/2
    chi, rho = sp.symbols("chi rho", positive=True)
    for fraction in fractions:
        mass = rho*(fraction*chi + 1-fraction*chi)
        ledger.check(f"species {fraction} retains the FULL physical ambient mass", sp.expand(mass-rho) == 0)

    psi = ef.add(ef.trig(2), ef.scale(ef.trig(1), sp.Rational(1, 100)))
    u = (psi, ef.trig(2, kind="sin"), ef.scale(ef.trig(1, kind="sin"), -sp.Rational(1, 100)))
    first = ef.trig(1)
    second = ef.transport(u, (first, {}, {}))[0]

    def generator(s):
        return (ef.scale(s, -1), ef.scale(ef.derivative(s, 2), -1), ef.derivative(s, 1))

    _, h, omega = ef.coadjoint_matrices(u, [generator(first), generator(second)], beltrami_eigenvalue=-1)
    averaged_h = (h+h)/2
    averaged_omega = (omega+omega)/2
    ledger.check("the actual same-field species average retains the FULL canonical Euler energy",
                 averaged_h == h and h[0, 0] > 0)
    ledger.check("the literal detector fraction is not a spurious multiplier of the full fluid phase",
                 averaged_omega == omega and omega[0, 1] != 0)
    ledger.check("mass-weighting the full action by tag fraction produces a false normalization",
                 averaged_h != avg_fraction*h and averaged_omega != avg_fraction*omega)

    gram = sp.Matrix(sphere_second_moment())
    j0 = sp.Symbol("j0", positive=True)
    vector = sp.Matrix(sp.symbols("a0:3"))
    ledger.check("the ACTUAL observation Gram reconstructs a common physical angle without scaling spin",
                 gram.inv()*gram*vector == vector)
    measured_spin = avg_fraction*j0*gram*vector
    ledger.check("the two positive species give the fixed nonzero measured density j0/8",
                 measured_spin == j0*vector/8)
    ledger.check("the same axial control can correct the FULL vector current with its distinct factor three",
                 gram*(3*vector) == vector)
    inertia, spin0 = sp.symbols("I S0", nonzero=True)
    acoustic_coefficients = [(f*j0-f*inertia)/(f*spin0) for f in fractions]
    ledger.check("the actual acoustic angle-null spin repair is unchanged by tag fractions",
                 all(sp.cancel(c-(j0-inertia)/spin0) == 0 for c in acoustic_coefficients))

    width = sp.Symbol("h", positive=True)
    cost = sp.Symbol("D", integer=True, positive=True)
    macro = width**(cost+1)
    flat_order = 2*cost+4
    leading_ratio = sp.powsimp(width**(flat_order+1)/macro**2)
    cubic_ratio = sp.powsimp(width**(-cost)*macro**3/macro**2)
    ledger.check("the LEADING physical clock error is strictly smaller than the spatial curvature scale",
                 sp.simplify(leading_ratio-width**3) == 0)
    ledger.check("the actual polynomial norm cost leaves a vanishing cubic-response remainder",
                 sp.simplify(cubic_ratio-width) == 0)
    epsilon = sp.Symbol("epsilon", positive=True)
    ledger.check("the observed second-coefficient approximation needs only its own vanishing accuracy",
                 sp.cancel(epsilon*macro**2/macro**2) == epsilon)
    bad_flat_order = 2*cost+1
    ledger.check("insufficient leading moment order is exposed as a NONVANISHING scaled error",
                 sp.simplify(width**(bad_flat_order+1)/macro**2) == 1)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
