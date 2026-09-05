"""Exact full homogeneous column and physical ring-tilt anchors."""

import sympy as s

from substrate_framework.euler_phase import physical_scalar_chart
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0217-closed-ring-material-clock")
    th = s.symbols("theta", real=True)
    r, q, z, wp, op, m, aa, ap = s.symbols("r q Z Wprime Oprime m A Aprime", real=True)
    alpha = m*aa/r
    x = s.Matrix([alpha*s.sin(th), ap*s.cos(th), 0])
    xt = s.Matrix([-q*alpha*s.cos(th),
                   (q*ap+r*op*alpha)*s.sin(th), wp*alpha*s.sin(th)])
    omega = s.Matrix([0, -wp, z])
    avg = s.integrate(omega.dot(x.cross(xt)), (th, 0, 2*s.pi))/(2*s.pi)
    checks.check("full homogeneous energy retains axial shear and radial strain",
                 s.simplify(avg-z*q*alpha*ap-(wp**2+z*r*op)*alpha**2/2) == 0)
    xs = s.Matrix([-alpha*s.cos(th), ap*s.sin(th), 0])
    checks.check("real quadrature phase has the actual vorticity cross sign",
                 s.simplify(-omega.dot(x.cross(xs))+z*alpha*ap) == 0)
    zp, k, w = s.symbols("Zprime k W", real=True)
    reduced = -m*(zp*q+z*(m*op+k*wp))/2+m**2*(wp**2+z*r*op)/(2*r)
    checks.check("integration by parts cancels only the actual planar shear term",
                 s.simplify(reduced-(-m*zp*q-m*k*z*wp+m**2*wp**2/r)/2) == 0)
    lam, om = s.symbols("lambda Omega", positive=True)
    beltrami = reduced.subs({z: lam*w, zp: lam*wp, q: m*om+k*w})
    checks.check("literal Beltrami relation yields complete energy minus beta times q",
                 s.simplify((beltrami+m*lam*wp*(m*om+k*w)).subs(wp, -lam*r*om)) == 0)
    checks.check("dropping the axial rate loses a genuinely nonzero energy term",
                 s.simplify(avg-avg.subs(wp, 0)-wp**2*alpha**2/2) == 0)

    rr = s.symbols("R", positive=True)
    radial = rr+r*s.cos(th)
    # The exact meridional preparation has R/r_cyl; it cancels the
    # cylindrical measure, leaving the following angular integrals.
    checks.check("literal tilt integral retains the physical radius factor",
                 s.simplify(s.integrate(radial*s.sin(th)*s.exp(-s.I*th),
                                         (th, 0, 2*s.pi))+s.I*s.pi*rr) == 0)
    spinweight = s.I*om*r**2+s.I*rr*om*r*s.cos(th)+r*w*s.sin(th)
    checks.check("actual Euclidean transverse spin retains both circulation clocks",
                 s.simplify(s.integrate(spinweight*s.exp(-s.I*th),
                                         (th, 0, 2*s.pi))-s.I*s.pi*r*(rr*om-w)) == 0)
    checks.check("linked physical centroid is nonzero and carries its separate moment",
                 s.simplify(s.integrate(radial*s.exp(-s.I*th),
                                         (th, 0, 2*s.pi))-s.pi*r) == 0)
    ix, iz, ax, ay = s.symbols("Ix Iz ax ay", real=True)
    rot = s.Matrix([[0, 0, ay], [0, 0, -ax], [-ay, ax, 0]])
    cov = s.diag(ix, ix, iz)
    dc = rot*cov-cov*rot
    checks.check("global covariance tilt has exact unit rigid-rotation response",
                 s.simplify(dc[0, 2]+s.I*dc[1, 2]-s.I*(ix-iz)*(ax+s.I*ay)) == 0)
    beta, nu, ct, time = s.symbols("beta nu c t", positive=True)
    row = s.Matrix([[ct*s.sin(nu*time), -ct*s.cos(nu*time)]])
    chart = physical_scalar_chart(
        beta*s.Matrix([[0, 1], [-1, 0]]), s.zeros(2), row,
        angle_rate=row.diff(time), angle_acceleration=row.diff(time, 2),
        generator_rate=s.zeros(2), spin=beta/(nu*ct**2)*row.diff(time))
    checks.check("actual angle history gives positive phase mass with its full time map",
                 s.simplify(chart.mass-beta/(nu*ct**2)) == 0)
    checks.check("positive scalar stiffness follows the actual observed rotation",
                 s.simplify(chart.stiffness-beta*nu/ct**2) == 0)
    rho, ic, ia, bb = s.symbols("rho Ic IA B", positive=True)
    phase = 2*s.pi**2*rho*rr*lam**2*om*ia
    inertia = 2*s.pi**2*rho*rr**3*ic
    cangle = bb/(2*rr*ic)
    checks.check("tag fraction is derived from physical phase and spin integrals",
                 s.simplify(phase/(om*cangle**2*inertia)-4*lam**2*ic*ia/bb**2) == 0)
    a, ci, ai, bi = s.symbols("a c0 a0 b0", positive=True)
    checks.check("fixed smooth profile scaling makes the positive tag fraction admissible",
                 s.simplify((4*lam**2*ic*ia/bb**2).subs(
                     {ic: a**2*ci, ia: a**4*ai, bb: a**2*bi})
                     -4*(lam*a)**2*ci*ai/bi**2) == 0)
    inverse = chart.coordinates.inv()
    conserved = inverse.T*(2*beta*nu*s.eye(2))*inverse
    checks.check("complete conserved energy retains the exposed factor two",
                 s.simplify(conserved-2*s.diag(chart.stiffness, chart.mass)) == s.zeros(2))
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
