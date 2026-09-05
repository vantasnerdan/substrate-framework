"""C-CST-002: conditional Rankine Euler mode equations and branch asymptotics.

Cartesian-derived velocity forms, Poincare pressure reduction and exact Bessel
boundary equation use exp(i*(m theta+k z-omega t)). At x=|k|a -> 0:
m=1: omega/Omega = -x^2/2*(log(2/x)+1/4-EulerGamma)+o(x^2);
m=2: omega/Omega = 1-x^2/6+O(x^4 log x).
The logarithmic m=2 remainder is not a constant times x^4.

Attempt 0029 restores the original velocity signs: 0019's alleged correction
reversed Coriolis terms. Attempt 0031 differentiates K1 exactly, removing
the spurious 1/4 in the old bending constant. Signed roots and 40/60-digit
refinement replace coarse grid minima. Analytic specification and error
budget: attempts/0031/README.md. These modes do not establish rod bend/twist
energy, inertia, or an Euler-to-Cosserat action.
"""
import mpmath as mp
import sympy as s

from substrate_framework.rankine_modes import boundary_determinant, core_velocity, rankine_residual
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("C-CST-002")
    r = s.Symbol("r", positive=True)
    Om, rho, wt, axial_k = s.symbols("Omega rho wt k", nonzero=True)
    azimuthal_m, lam2 = s.symbols("m lambda_squared")
    pressure = s.Function("P")(r)
    vr, vt, vz = core_velocity(pressure, r, axial_k, azimuthal_m, wt, Om, rho)
    divergence = s.diff(r*vr, r)/r+s.I*azimuthal_m*vt/r+s.I*axial_k*vz
    reduced = s.simplify(divergence.subs(
        s.diff(pressure, r, 2),
        -s.diff(pressure, r)/r+(azimuthal_m**2/r**2-lam2)*pressure))
    ledger.check("Poincare reduction from solved Euler velocity",
                 s.simplify(reduced.subs(lam2, axial_k**2*(4*Om**2-wt**2)/wt**2)) == 0)
    wrong_reduced = reduced.subs(lam2, -axial_k**2*(4*Om**2-wt**2)/wt**2)
    ledger.check("wrong pressure species breaks incompressibility",
                 s.simplify(wrong_reduced) != 0)
    x = s.Symbol("x", positive=True)
    ell, gamma, c = s.symbols("L gamma c", real=True)
    K1 = 1/x+x/2*(s.log(x/2)+gamma-s.Rational(1, 2))
    derivative = s.diff(K1, x)
    ledger.check("K1 derivative retains the missing quarter",
                 s.simplify(derivative-(-1/x**2+(s.log(x/2)+gamma+s.Rational(1, 2))/2)) == 0)
    Kratio = s.series(x*derivative/K1, x, 0, 3).removeO().expand()
    ledger.check("K1 logarithmic derivative has no spurious constant",
                 s.simplify(Kratio-(-1+x**2*(s.log(x/2)+gamma))) == 0)
    Jratio = 1-3*x**2/4
    w = -x**2/2*(ell+c)
    F = boundary_determinant(-1+w, 1, 1, Jratio, Kratio.subs(s.log(x), s.log(2)-ell))
    leading = s.expand(s.series(F, x, 0, 3).removeO()).coeff(x, 2)
    constant = s.solve(leading, c)[0]
    ledger.check("m1 boundary equation derives c=1/4-gamma",
                 s.simplify(constant-(s.Rational(1, 4)-gamma)) == 0)
    ledger.check("old c=1/2-gamma fails the defining boundary equation",
                 s.simplify(leading.subs(c, s.Rational(1, 2)-gamma)) != 0)
    F2 = boundary_determinant(-1+c*x**2, 1, 2, 2-x**2/2, -2-x**2/2)
    m2_coefficient = s.solve(s.expand(s.series(F2, x, 0, 3).removeO()).coeff(x, 2), c)[0]
    ledger.check("m2 branch retains its independently derived -1/6 coefficient",
                 m2_coefficient == -s.Rational(1, 6))
    ledger.check("exact translation neutrality at k=0",
                 boundary_determinant(-1, 1, 1, 1, -1) == 0)
    print("m1 constant:", constant, "; m2 coefficient:", m2_coefficient)
    # Irreducible remainder: exact Bessel root vs asymptotic series.
    for m in (1, 2):
        previous_error = None
        for exponent in (2, 3, 4):
            roots = []
            for digits in (40, 60):
                ctx = mp.mp.clone()
                ctx.dps = digits
                xx = ctx.mpf(10)**(-exponent)
                shift = (-xx**2/2*(ctx.log(2/xx)+ctx.mpf(1)/4-ctx.euler)
                         if m == 1 else -xx**2/6)
                def residual(ss):
                    return rankine_residual(ctx, xx, m, ss)
                root = ctx.findroot(residual, (-1+shift/2, -1+shift*2),
                                    solver="anderson", tol=ctx.mpf(10)**(-digits+8),
                                    maxsteps=200)
                frequency = root+m
                remainder = abs((frequency-(m-1)-shift)/xx**2)
                roots.append(ctx.nstr(root, digits))
                ledger.check(f"m={m} x=1e-{exponent} dps={digits} boundary residual",
                             abs(residual(root)) < ctx.mpf(10)**(-digits+10))
            precision_error = abs(ctx.mpf(roots[0])-ctx.mpf(roots[1]))
            ledger.check(f"m={m} x=1e-{exponent} precision refinement",
                         precision_error < ctx.mpf("1e-30"))
            if previous_error is not None:
                ledger.check(f"m={m} x=1e-{exponent} asymptotic remainder converges",
                             remainder < previous_error/10)
            previous_error = remainder
            print(f"m={m} x=1e-{exponent}: scaled remainder={ctx.nstr(remainder, 8)} "
                  f"precision discrepancy={ctx.nstr(precision_error, 5)}")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
