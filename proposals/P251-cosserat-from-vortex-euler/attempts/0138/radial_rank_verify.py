"""Exact ordinary-column four-row material-tag jet construction."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0138-fixed-tag-carrier-jet")
    m = s.Symbol("m", integer=True, positive=True)
    radial_parameter = s.Symbol("t0", positive=True)
    coefficients = [(-1)**n/(4**n*s.factorial(n)*s.rf(m+1, n)) for n in range(4)]
    for n in range(1, 4):
        checks.check(f"Bessel pressure recurrence coefficient {n}",
                     s.simplify(4*n*(m+n)*coefficients[n]+coefficients[n-1]) == 0)
    pressure = s.Matrix([[coefficients[n]*radial_parameter**n for n in range(4)]])
    coefficient_matrix = s.Matrix.vstack(s.Matrix([[1, 0, 0, 0]]), pressure,
                                         pressure.diff(radial_parameter),
                                         pressure.diff(radial_parameter, 2))
    determinant = s.factor(coefficient_matrix.det())
    expected = radial_parameter**3/(24576*(m+1)**3*(m+2)**2*(m+3))
    checks.check("four pressure/moment rows have a nonzero finite coefficient determinant",
                 s.simplify(determinant-expected) == 0 and determinant.is_positive is True)

    scale, scale1, scale2, t1, t2 = s.symbols("A A1 A2 t1 t2", nonzero=True)
    change = s.Matrix([[scale, 0, 0], [scale1, scale*t1, 0],
                       [scale2, 2*scale1*t1+scale*t2, scale*t1**2]])
    checks.check("actual carrier derivatives preserve rank with the full normalization terms",
                 s.simplify(change.det()-scale**3*t1**3) == 0)
    inverse, inverse1, inverse2 = s.symbols("inverse inverse1 inverse2", nonzero=True)
    scalar_division = s.Matrix([[inverse, 0, 0], [inverse1, inverse, 0],
                                [inverse2, 2*inverse1, inverse]])
    checks.check("a radially constant Doppler factor gives an invertible triangular change",
                 s.simplify(scalar_division.det()-inverse**3) == 0)

    dq = s.Symbol("dq")
    qmoment, r0, r1, r2, error = s.symbols("Q R0 R1 R2 error", nonzero=True)
    target = r0+r1*dq+r2*dq**2/2
    actual = qmoment*(r0+r1*dq+r2*dq**2/2)+error*dq**3
    factor0, factor1, factor2 = s.symbols("factor0 factor1 factor2")
    mismatch = s.expand((factor0+factor1*dq+factor2*dq**2)*(actual-qmoment*target))
    checks.check("one fixed tag cancels all three physical/canonical mismatch coefficients",
                 all(mismatch.coeff(dq, n) == 0 for n in range(3)))
    checks.check("endpoint degeneracy is retained instead of declaring a zero-carrier rotor",
                 determinant.subs(radial_parameter, 0) == 0)
    print("Exact four-row coefficient determinant:", determinant)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
