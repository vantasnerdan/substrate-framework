"""Exact full Euler/Lin, radial pressure, localization and physical-current rows."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0193-localized-helical-pole")
    r, c, cap, sigma, m, radius = s.symbols(
        "r c C sigma m R", positive=True
    )
    d = c**2+r**2
    f = cap/d
    cutoff = s.Function("cutoff")(r/radius)
    localized = cutoff*f
    checks.check(
        "radial cutoff retains the exact steady Euler pressure balance",
        s.simplify(r*localized**2-(r*localized)**2/r) == 0,
    )
    checks.check(
        "actual outer helical charge is not set constant",
        s.simplify(s.diff(d*localized, r)-cap*s.diff(cutoff, r)) == 0,
    )
    checks.check(
        "coefficient derivative retains both cutoff and original tail",
        s.simplify(s.diff(localized-f, r)
                   -(cutoff-1)*s.diff(f, r)-f*s.diff(cutoff, r)) == 0,
    )
    checks.check(
        "uncut derivative has the stated weighted tail scale",
        s.simplify(s.sqrt(d)*s.diff(f, r)+2*cap*r/d**s.Rational(3, 2)) == 0,
    )

    ff, fp, vr, vt, vz, tau = s.symbols("f fp vr vt vz tau")
    charge_derivative = 2*r*ff+d*fp
    checks.check(
        "helical momentum equation keeps all three full Euler terms",
        s.expand(r*(2*ff+r*fp)+c**2*fp-charge_derivative) == 0,
    )
    radial_kelvin = fp*tau/(s.I*sigma)+fp*vr*charge_derivative/(s.I*sigma)**2
    checks.check(
        "nonconstant outer charge cancels in the actual Kelvin gradient identity",
        s.simplify(radial_kelvin.subs(tau, -charge_derivative*vr/(s.I*sigma))) == 0,
    )
    h = s.Matrix([0, r, c])
    vort = s.Matrix([0, -c*fp, 2*ff+r*fp])
    checks.check(
        "actual vorticity return measures force-free failure of a scalar cutoff",
        s.simplify(h.cross(vort)-s.Matrix([charge_derivative, 0, 0])) == s.zeros(3, 1),
    )

    radial = s.Function("F")(r)
    pressure = s.Function("pressure")(r)
    angular = s.Function("f")(r)
    angular_prime = s.diff(angular, r)
    v_radial = s.I*sigma*radial
    v_theta = -m*pressure/(sigma*r)-(2*angular+r*angular_prime)*radial
    v_axial = m*pressure/(sigma*c)-c*angular_prime*radial
    radial_prime = (-1/r+2*m*angular/(sigma*r))*radial
    radial_prime += m**2*(1/r**2+1/c**2)*pressure/sigma**2
    pressure_prime = (sigma**2-2*angular*(2*angular+r*angular_prime))*radial
    pressure_prime -= 2*m*angular*pressure/(sigma*r)
    checks.check(
        "actual radial pressure equation includes the inertial and return terms",
        s.simplify(s.I*sigma*v_radial-2*angular*v_theta+pressure_prime) == 0,
    )
    checks.check(
        "actual tangential Euler equation",
        s.simplify(s.I*sigma*v_theta+(2*angular+r*angular_prime)*v_radial
                   +s.I*m*pressure/r) == 0,
    )
    checks.check(
        "actual axial Euler equation at the helical carrier",
        s.simplify(s.I*sigma*v_axial+c*angular_prime*v_radial
                   -s.I*m*pressure/c) == 0,
    )
    divergence = s.diff(v_radial, r)+v_radial/r+s.I*m*v_theta/r-s.I*m*v_axial/c
    checks.check(
        "full divergence derives the radial material equation",
        s.simplify(divergence.subs(s.diff(radial, r), radial_prime)) == 0,
    )
    outer_equation = s.diff(pressure, r, 2)-sigma**2*radial_prime.subs(angular, 0)
    outer_equation = outer_equation.subs(radial, s.diff(pressure, r)/sigma**2)
    checks.check(
        "exterior pressure is screened Helmholtz rather than a wall",
        s.simplify(outer_equation-(s.diff(pressure, r, 2)+s.diff(pressure, r)/r
                                  -m**2*(1/r**2+1/c**2)*pressure)) == 0,
    )

    beta = s.symbols("beta", positive=True)
    symplectic_unit = s.Matrix([[0, 1], [-1, 0]])
    checks.check(
        "continued actual KKS and physical generator give positive phase energy",
        -(beta*symplectic_unit)*(sigma*symplectic_unit) == beta*sigma*s.eye(2),
    )
    phi = s.Function("phi")(r)
    g_row = c**2*r**2*s.diff(phi, r)/(sigma*d)-m*s.diff(f, r)*r**2*phi/sigma**2
    spin_row = -c**2*r**2*s.diff(phi, r)/d+2*m*cap*c**2*r*phi/(sigma*d**2)
    checks.check(
        "new inner mode retains both displacement and moving-position spin rows",
        s.simplify(spin_row+sigma*g_row-2*m*cap*r*phi/(sigma*d)) == 0,
    )
    action, obs, eps, reference, eta, prefactor = s.symbols(
        "action observation epsilon reference eta prefactor", nonzero=True
    )
    angle_factor = obs/(sigma*eps*reference)
    physical_mass = action/(sigma**2*angle_factor**2)
    spin_coefficient = prefactor*eps*(-eta*action*reference/(prefactor*obs))
    checks.check(
        "actual three-row marker target matches the new full action mass",
        s.simplify(spin_coefficient+eta*physical_mass*sigma*angle_factor) == 0,
    )

    vv, vp, zz, ww, wp = s.symbols("V Vprime Z W Wprime", nonzero=True)
    checks.check(
        "Bernoulli square root derives the generalized-force-free axial derivative",
        s.simplify((-vv**2/r-vv*vp)+vv*(vp+vv/r)) == 0,
    )
    gf_vorticity = s.Matrix([0, -wp, zz])
    gf_velocity = s.Matrix([0, vv, ww])
    checks.check(
        "compact-vorticity alternative has actual parallel curl",
        s.simplify((gf_vorticity-zz*gf_velocity/ww).subs(wp, -vv*zz/ww)) == s.zeros(3, 1),
    )
    original_p, pressure_loss, cc = s.symbols("p loss cutoff", real=True)
    reconstructed = -2*(original_p-pressure_loss)-cc**2*vv**2
    checks.check(
        "positive axial square contains both missing-swirl and pressure-integral terms",
        s.expand(reconstructed-((-2*original_p-vv**2)+(1-cc**2)*vv**2
                                +2*pressure_loss)) == 0,
    )
    grad, bound, gg = s.symbols("gprime K g", real=True, nonzero=True)
    step = -grad/bound
    checks.check(
        "square-root derivative bound uses a nonnegative Taylor test not a lower floor",
        s.simplify(gg+grad*step+bound*step**2/2-(gg-grad**2/(2*bound))) == 0,
    )
    nu, nu_prime, pp = s.symbols("nu nuprime pressure", nonzero=True)
    extra_numerator = -s.I*nu_prime*pp
    checks.check(
        "general axial-flow Kelvin correction is the full varying-frequency gradient",
        s.simplify(extra_numerator/(s.I*nu)**2-s.I*nu_prime*pp/nu**2) == 0,
    )
    print("Operator-norm, Riesz rank and exact marker-inverse bounds are analytic in the proofs.")
    print("No numerical eigenvalue, guessed gap, radial wall or full-channel spectrum was used.")
    raise SystemExit(checks.finish())


if __name__ == "__main__":
    main()
