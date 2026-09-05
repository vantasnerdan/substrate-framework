"""Exact ABC Euler response; no magnetic or mean-field dynamic closure."""

import sympy as s

from substrate_framework.euler_fourier import (
    add, coadjoint_matrices, cross, curl, derivative, inner, leray, mul, scale, trig,
)
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0132-Euler-magnetic-reconciliation")
    # Length unit alpha^-1, velocity amplitudes U_i; 0 < |k| < 1/3.
    # Full real Fourier fields, bilinear volume average, density one.
    k = s.Symbol("k", real=True, nonzero=True)
    u1, u2, u3 = s.symbols("U1 U2 U3", real=True)
    background = (
        add(scale(trig(2), u3), scale(trig(1, kind="sin"), u2)),
        add(scale(trig(0), u1), scale(trig(2, kind="sin"), u3)),
        add(scale(trig(1), u2), scale(trig(0, kind="sin"), u1)),
    )
    plus, minus = (0, 0, k), (0, 0, -k)
    eta = ({plus: s.Rational(1, 2), minus: s.Rational(1, 2)},
           {plus: 1/(2*s.I), minus: -1/(2*s.I)}, {})
    eta_quadrature = (scale(eta[1], -1), eta[0], {})
    velocity, hessian, kks = coadjoint_matrices(
        background, [eta, eta_quadrature], beltrami_eigenvalue=-1)
    v = velocity[0]
    expected_h = -(u1**2+u2**2)*k**3/(2*(1+k**2))
    print("Actual coadjoint Hessian:", hessian[0, 0])
    checks.check("full Euler ABC Hessian reproduces its signed cubic response",
                 s.simplify(hessian[0, 0]-expected_h) == 0)
    checks.check("axial background component does not enter this displacement energy",
                 s.diff(hessian[0, 0], u3) == 0)
    checks.check("opposite spatial helicity reverses the exact Hessian",
                 s.simplify(hessian[0, 0].subs(k, -k)+hessian[0, 0]) == 0)
    checks.check("this bare displacement has no quadratic long-wave modulus",
                 s.limit(hessian[0, 0]/k**2, k, 0) == 0)

    stress = [[add(mul(background[i], v[j]), mul(v[i], background[j]))
               for j in range(3)] for i in range(3)]
    force_full = leray(tuple(scale(add(*(derivative(stress[i][j], j)
                                        for j in range(3))), -1)
                             for i in range(3)))
    force = tuple({wave: s.factor(value) for wave, value in component.items()
                   if wave in (plus, minus)} for component in force_full)
    print("Actual instantaneous slow Euler acceleration:", force)
    target = (scale(eta[0], k**2*((u1**2+u2**2)*k+u1**2-u2**2)/(2*(1+k**2))),
              scale(eta[1], k**2*((u1**2+u2**2)*k-u1**2+u2**2)/(2*(1+k**2))), {})
    checks.check("complete pressure gives the exact slow Reynolds acceleration",
                 all(s.simplify(value) == 0 for i in range(3)
                     for value in add(force[i], scale(target[i], -1)).values()))
    # Independent Cartesian pressure calculation for U1 only.
    x, z = s.symbols("x z", real=True)
    pressure_potential = -(1-k)/(1+k**2)*s.sin(k*z)*s.cos(x)
    ax, ay, az = s.sin(k*z)*s.sin(x), -s.cos(k*z)*s.sin(x), s.cos(k*z)*s.cos(x)
    vx = -ax+s.diff(pressure_potential, x)
    vy = -ay
    vz = -az+s.diff(pressure_potential, z)
    checks.check("Cartesian Poisson correction satisfies the full divergence equation",
                 s.simplify(s.diff(vx, x)+s.diff(vz, z)) == 0)
    t13 = s.integrate(s.sin(x)*vx, (x, 0, 2*s.pi))/(2*s.pi)
    t23 = s.integrate(s.cos(x)*vz+s.sin(x)*vy, (x, 0, 2*s.pi))/(2*s.pi)
    checks.check("independent physical pressure reproduces anisotropic acceleration",
                 s.simplify(-s.diff(t13, z)-k**2*(k+1)*s.cos(k*z)/(2*(1+k**2))) == 0
                 and s.simplify(-s.diff(t23, z)-k**2*(k-1)*s.sin(k*z)/(2*(1+k**2))) == 0)
    checks.check("instantaneous virtual work equals minus the Euler Hessian",
                 s.simplify(inner(eta, force)+hessian[0, 0]) == 0)

    # Same equilibrium field, distinct transported field under perturbation.
    advected_velocity = curl(cross(eta, background))
    magnetic_h = s.factor(inner(advected_velocity, advected_velocity)
                         - inner(advected_velocity, v))
    # Magnetic Hessian is |curl(eta cross u)|² - v.curl(eta cross u).
    difference = tuple(add(advected_velocity[i], scale(v[i], -1)) for i in range(3))
    checks.check("sum of distinct Hessians is the squared generator mismatch",
                 s.simplify(magnetic_h+hessian[0, 0]-inner(difference, difference)) == 0)
    print("Magnetic Hessian with same base field:", magnetic_h)
    checks.check("magnetic field transport differs from Kelvin velocity transport",
                 s.simplify(inner(difference, difference)) != 0)
    checks.check("the two physical macro displacement columns have degenerate KKS",
                 kks == s.zeros(2) and hessian.det() != 0)
    checks.check("importing the magnetic Hessian would change the actual Euler force",
                 s.simplify(inner(eta, force)+magnetic_h) != 0)
    print("Scope: exact energy and initial Euler acceleration, not closed mean dynamics")
    raise SystemExit(checks.finish())


if __name__ == "__main__":
    main()
