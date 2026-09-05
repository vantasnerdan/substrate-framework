"""Helicity insertion, exact periodic Euler and normalized defect-energy identities."""

import sympy as s

from substrate_framework import euler_fourier as fourier
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0153-same-field")
    a, b = s.symbols("a b", real=True)
    n = s.Matrix([2*a, 2*b, 1-a*a-b*b])/(1+a*a+b*b)
    cross = s.Matrix([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    projection = (s.eye(3)-n*n.T+s.I*cross)/2
    checks.check("rational-direction positive helicity is an exact Hermitian projector",
                 s.simplify(projection.H-projection) == s.zeros(3)
                 and s.simplify(projection*projection-projection) == s.zeros(3))
    checks.check("projection enforces the actual divergence and curl equations",
                 s.simplify(n.T*projection) == s.zeros(1, 3)
                 and s.simplify(s.I*cross*projection-projection) == s.zeros(3))
    checks.check("antipodal helicity projection preserves real Fourier fields",
                 s.simplify((s.eye(3)-n*n.T-s.I*cross)/2-s.conjugate(projection))
                 == s.zeros(3))

    reference = ({}, {}, {})
    north, south = (0, 0, 1), (0, 0, -1)
    ref_coeff = s.Matrix([s.Rational(1, 2), s.I/2, 0])
    for i in range(3):
        if ref_coeff[i]:
            reference[i][north] = ref_coeff[i]
            reference[i][south] = s.conjugate(ref_coeff[i])
    defect = ({}, {}, {})
    energies = []
    e1 = s.Matrix([1, 0, 0])
    for aa, bb in ((s.Rational(1, 10), 0),
                   (s.Rational(1, 7), s.Rational(1, 11)),
                   (-s.Rational(1, 9), s.Rational(1, 8))):
        wave = tuple(n.subs({a: aa, b: bb}))
        coeff = s.simplify((projection*e1).subs({a: aa, b: bb}))/3
        energies.append(2*s.simplify(coeff.H.dot(coeff)))
        for i in range(3):
            if coeff[i]:
                defect[i][wave] = coeff[i]
                defect[i][tuple(-v for v in wave)] = s.conjugate(coeff[i])
    joined = tuple(fourier.add(reference[i], defect[i]) for i in range(3))
    actual_curl = fourier.curl(joined)
    checks.check("one joined periodic field retains exactly the same curl",
                 not fourier.divergence(joined)
                 and all(not fourier.add(actual_curl[i], fourier.scale(joined[i], -1))
                         for i in range(3)))
    norm2 = fourier.add(*(fourier.mul(v, v) for v in joined))
    pressure = fourier.scale(norm2, -s.Rational(1, 2))
    convection = fourier.transport(joined, joined)
    checks.check("all nonlinear cross modes cancel against the full Euler pressure",
                 all(not fourier.add(convection[i], fourier.derivative(pressure, i))
                     for i in range(3)))
    mean = norm2.get((0, 0, 0), 0)
    checks.check("actual normalized energy is reference plus distinct-mode defect energy",
                 s.simplify(mean-1-sum(energies)) == 0)
    gradient_norm = fourier.add(*(
        fourier.mul(fourier.derivative(v, j), fourier.derivative(v, j))
        for v in defect for j in range(3)))
    checks.check("fixed-shell derivative energy has no growing period factor",
                 s.simplify(gradient_norm.get((0, 0, 0), 0)-sum(energies)) == 0)
    reference_pressure = fourier.scale(fourier.add(*(
        fourier.mul(v, v) for v in reference)), -s.Rational(1, 2))
    checks.check("omitting insertion cross-pressure is exposed by the actual Euler field",
                 any(fourier.add(convection[i], fourier.derivative(reference_pressure, i))
                     for i in range(3)))
    # The defect is locally coherent even though its normalized energy is
    # a sum of squared distinct Fourier coefficients, not their square sum.
    at_origin = s.Matrix([sum(component.values()) for component in defect])
    checks.check("coherent local value and averaged energy are genuinely different norms",
                 at_origin != s.zeros(3, 1)
                 and s.simplify(at_origin.dot(at_origin)-sum(energies)) != 0)
    print("Scope: exact helicity/periodic Euler/energy identities supporting the analytic join;")
    print("prescribed tubes use local approximation; uniform long-wave closure remains separate")
    raise SystemExit(checks.finish())


if __name__ == "__main__":
    main()
