"""Exact finite-cylinder Neumann repair of the polynomial pressure limit."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0133-finite-cylinder-pressure")
    r, radius, axial = s.symbols("r R q", positive=True)
    denominator = axial*(s.besseli(1, axial*radius)+s.besseli(3, axial*radius))/2
    correction = -2*radius*s.besseli(2, axial*r)/denominator
    pressure = r*r+correction
    checks.check("full Neumann correction removes the polynomial wall-normal return",
                 s.simplify(s.diff(pressure, r).subs(r, radius)) == 0)
    small = s.series(pressure, axial, 0, 4).removeO()
    checks.check("finite-cylinder pressure has a regular zero-axial-wave-number limit",
                 s.simplify(small-axial**2*r*r*(2*radius**2-r*r)/12) == 0)
    radial_operator = s.diff(r*r, r, 2)+s.diff(r*r, r)/r-(4/r**2+axial**2)*r*r
    checks.check("polynomial particular pressure solves the correct azimuth-two source",
                 s.simplify(radial_operator+axial**2*r*r) == 0)
    print("Leading Neumann-corrected pressure:", s.factor(small))
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
