"""Correct the crossing functional and execute the finite-segment repair.

Exact algebra: geometry and second variation are computed before integration.
Finite open segments give a local interaction contribution, not a closed
divergence-free vortex configuration. Its endpoint completion remains explicit.
"""
import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0033-mutual-energy")
    theta, r, q = s.symbols("theta r q", real=True)
    d = s.Symbol("d", positive=True)
    t1 = s.Matrix([0, 0, 1])
    t2 = s.Matrix([s.sin(theta), 0, s.cos(theta)])
    separation = s.Matrix([0, d, 0])+q*t2-r*t1
    distance2 = s.simplify(separation.dot(separation))
    ledger.check("declared crossing geometry gives the actual squared distance",
                 s.simplify(distance2-(d**2+r**2+q**2-2*r*q*s.cos(theta))) == 0)
    energy_kernel = t1.dot(t2)/s.sqrt(distance2)
    ledger.check("mutual-energy kernel is even under reflected tilt",
                 s.simplify(energy_kernel.subs(theta, -theta)-energy_kernel) == 0)
    circulation_numerator = s.simplify(t1.cross(separation).dot(t2))
    ledger.check("velocity circulation has a different, parity-odd numerator",
                 circulation_numerator == -d*s.sin(theta))
    ledger.check("velocity line integral is not the kinetic-energy kernel",
                 s.simplify(circulation_numerator/distance2**s.Rational(3, 2)-energy_kernel) != 0)
    second = s.simplify(s.diff(energy_kernel, theta, 2).subs(theta, 0))
    target = -(d**2+r**2+q**2-r*q)/(d**2+(r-q)**2)**s.Rational(3, 2)
    ledger.check("finite-segment tilt Hessian derives directly from mutual energy",
                 s.simplify(second-target) == 0)
    positive_numerator = d**2+(r-q)**2/2+(r**2+q**2)/2
    ledger.check("tilt Hessian sign has an exact positive-square decomposition",
                 s.expand(positive_numerator-(d**2+r**2+q**2-r*q)) == 0)
    print("E12'' kernel / (rho Gamma1 Gamma2 / 4pi) =", second)
    print("For Gamma1*Gamma2<0 this finite-segment contribution is strictly positive.")
    print("Closed-loop endpoint completion and the dynamical action remain required.")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
