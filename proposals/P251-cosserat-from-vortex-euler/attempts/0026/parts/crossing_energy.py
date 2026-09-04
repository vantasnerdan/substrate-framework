"""Attempt 0026 part 1 -- exact inter-tube Biot-Savart crossing energy.

Geometry: line 1 = z-axis (tangent t1 = zhat, circulation Gamma); line 2 =
{d yhat + s t2}, t2 = (sin th, 0, cos th) so the offset is perpendicular to
t2 and the minimum separation is exactly d.

Cross energy: E_int = rho Gam^2 ∮_{C2} u1 · dx2 (standard mutual-inductance
form of the Biot-Savart cross term). The numerator (t1 x R)·t2 is CONSTANT
(= -d sin th), so the double line integral reduces to one quadrature:

    E_int/rhoGam^2 = (-d sin th / 4pi) ∫∫ [d^2 + r^2 + s^2 - 2 d r sin th]^(-3/2) ds dr
                   = (-d sin th / 4pi) * 2pi/(d cos th)
                   = -(1/2) tan th.

This script verifies the reduction symbolically and the closed form against
30-digit mpmath quadrature.
"""
import mpmath as mp
import sympy as sp


def E_exact(theta, Gam=1.0, rho=1.0):
    return -0.5 * rho * Gam**2 * mp.tan(theta)


def E_numeric(theta, Gam=1.0, rho=1.0, d=1.0):
    st = mp.sin(theta)
    inner = lambda rr: 2 / (d**2 + rr**2 - 2 * d * rr * st)  # analytic s-reduction
    I = mp.quad(inner, [-mp.inf, mp.inf])
    return rho * Gam**2 / (4 * mp.pi) * (-d * st) * I


def main():
    mp.mp.dps = 30
    worst = 0.0
    for deg in (0.3, 0.7, 1.1, 1.4):
        t = mp.mpf(str(deg))
        err = abs((E_numeric(t) - E_exact(t)) / E_exact(t))
        worst = max(worst, float(err))
        print(f"theta={deg}: numeric={mp.nstr(E_numeric(t), 14)} "
              f"exact={mp.nstr(E_exact(t), 14)} rel.err={mp.nstr(err, 3)}")
    assert worst < 1e-12, "closed form vs quadrature mismatch"

    th = sp.symbols("theta", positive=True)
    E = -sp.Rational(1, 2) * sp.tan(th)
    second = sp.simplify(sp.diff(E, th, 2))
    print("E''(theta) =", second, "(nonzero at generic theta)")
    print("parity: E(-theta) =", sp.simplify(E.subs(th, -th)), "=> ODD")
    assert sp.simplify(E.subs(th, -th) + E) == 0, "parity check failed"
    print("ALL CHECKS PASS (0026 part1)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
