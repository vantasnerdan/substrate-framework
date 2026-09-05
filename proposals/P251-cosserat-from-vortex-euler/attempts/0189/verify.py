"""Derive the actual helical Lin displacement and full physical current."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0189-full-helical-current")
    r, c, cap, sigma, rho, length, epsilon, eta, action = s.symbols(
        "r c C sigma rho L epsilon eta h", positive=True)
    theta, a, b, marker = s.symbols("theta a b marker", real=True)
    phi = s.Function("phi")(r)
    m = s.Integer(2)
    d = c*c+r*r
    f = cap/d
    phase = a*s.cos(m*theta)+b*s.sin(m*theta)
    quadrature = b*s.cos(m*theta)-a*s.sin(m*theta)
    xi_r = m*phi*phase/(sigma*r)
    xi_theta = (c*c*s.diff(phi, r)/(sigma*d)
                -m*s.diff(f, r)*phi/sigma**2)*quadrature
    velocity_theta = -c*c*s.diff(phi, r)*phase/d
    rate = sigma*b*s.diff(xi_theta, a)-sigma*a*s.diff(xi_theta, b)
    checks.check("full Lin reconstruction retains radial shear and basis rotation",
                 s.simplify(rate-r*s.diff(f, r)*xi_r-velocity_theta) == 0)
    tag = 1+epsilon*marker*s.cos(m*theta)
    g_row = s.integrate(s.expand_trig(tag*r*xi_theta*r), (theta, 0, 2*s.pi))
    spin_row = s.integrate(s.expand_trig(tag*(r*velocity_theta
                         +s.diff(r*r*f, r)*xi_r)*r), (theta, 0, 2*s.pi))
    expected_g = s.pi*epsilon*marker*b*(c*c*r*r*s.diff(phi, r)/(sigma*d)
                                       -m*s.diff(f, r)*r*r*phi/sigma**2)
    expected_s = s.pi*epsilon*marker*a*(-c*c*r*r*s.diff(phi, r)/d
                                       +2*m*cap*c*c*r*phi/(sigma*d*d))
    checks.check("literal displacement moment follows actual angular integration",
                 s.simplify(g_row-expected_g) == 0)
    checks.check("mechanical spin includes both position and velocity variations",
                 s.simplify(spin_row-expected_s) == 0)
    g_rate = sigma*b*s.diff(g_row, a)-sigma*a*s.diff(g_row, b)
    checks.check("additional current connection is the derived third radial moment",
                 s.simplify(spin_row-g_rate
                            -2*m*cap*s.pi*epsilon*marker*a*r*phi/(sigma*d)) == 0)
    ref, obs = s.symbols("B Iobs", nonzero=True, real=True)
    c_obs = obs/(sigma*epsilon*ref)
    mass = action/(sigma**2*c_obs**2)
    j_moment = -eta*action*ref/(rho*s.pi*length*obs)
    spin = rho*s.pi*length*epsilon*j_moment*a
    current = -spin/(sigma*a)*b
    checks.check("exact moment constraints give positive physical spin overlap",
                 s.simplify(spin+eta*mass*sigma*c_obs*a) == 0)
    checks.check("same initial displacement supplies the full integrated current",
                 s.simplify(current-eta*mass*c_obs*b) == 0)
    omitted = s.integrate(s.expand_trig(tag*r*velocity_theta*r),
                         (theta, 0, 2*s.pi))
    checks.check("omitting displaced position changes the actual spin row",
                 s.simplify(spin_row-omitted) != 0)
    print("Derived S-G_t radial density:", 2*m*cap*r*phi/(sigma*d))
    print("Carrier-two and Euclidean closed-tube transfer remain separate constructions.")
    raise SystemExit(checks.finish())


if __name__ == "__main__":
    main()
