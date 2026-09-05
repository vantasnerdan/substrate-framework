"""Exact intrinsic dispersion, complete material rows, and reaction algebra."""

import sympy as s

from substrate_framework.rankine_modes import boundary_determinant, core_velocity
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0135-intrinsic-surface-response")
    m = s.Symbol("m", integer=True, positive=True)
    x, slope = s.symbols("x slope", real=True)
    doppler = -1+slope*x**2
    jr = m-3*x**2/(2*(m+1))
    kr = -m-x**2/(2*(m-1))
    residual = boundary_determinant(doppler, 1, m, jr, kr)
    coefficient = s.expand(residual).coeff(x, 2)
    root_slope = s.solve(coefficient, slope)[0]
    checks.check("surface root coefficient follows from the pressure matching determinant",
                 s.simplify(root_slope+1/(2*(m*m-1))) == 0)
    checks.check("opposite surface curvature fails the actual boundary residual",
                 s.simplify(coefficient.subs(slope, -root_slope)) != 0)
    lab2 = s.expand((m-1+root_slope*x*x)**2).coeff(x, 2)
    material2 = s.expand((-1+root_slope*x*x)**2).coeff(x, 2)
    checks.check("laboratory squared-frequency curvature is negative",
                 s.simplify(lab2+1/(m+1)) == 0)
    checks.check("transported-material squared-frequency curvature is positive",
                 s.simplify(material2-1/(m*m-1)) == 0)
    checks.check("m2 material curvature is one third", material2.subs(m, 2) == s.Rational(1, 3))

    r, om, rho = s.symbols("r Omega rho", positive=True)
    sigma, k = s.symbols("sigma k", nonzero=True, real=True)
    pressure = s.Function("P")(r)
    velocity = core_velocity(pressure, r, k, m, sigma, om, 1)
    eta = s.I*velocity/sigma
    spin = s.simplify(r*(2*om*eta[0]+velocity[1]))
    checks.check("moving-material spin retains both deformation and velocity",
                 s.simplify(spin-m*pressure/sigma) == 0)
    checks.check("velocity-only spin is exposed by the exact material identity",
                 s.simplify(r*velocity[1]-m*pressure/sigma) != 0)
    shape_kernel = s.simplify(eta[0]+eta[1]/s.I)
    checks.check("material shape numerator follows from the same full Euler field",
                 s.simplify(shape_kernel-(s.diff(pressure, r)+m*pressure/r)
                            /(sigma*(2*om+sigma))) == 0)

    beta, row = s.symbols("beta c", positive=True)
    j2 = s.Matrix([[0, -1], [1, 0]])
    symplectic = -beta*j2
    lab_frequency = m*om+sigma
    h_lab = -symplectic*(lab_frequency*j2)
    momentum = -symplectic*(m*j2)
    h_rot = s.simplify(h_lab-om*momentum)
    checks.check("physical SO2 momentum retains m and its material connection",
                 h_rot == -beta*sigma*s.eye(2))
    canonical_mass = -beta/(sigma*row**2)
    checks.check("intrinsic scalar action has positive mass on sigma negative",
                 canonical_mass.subs(sigma, -om).is_positive is True)

    mu, eps, z0, zc, q, t, g = s.symbols("mu eps Z0 Zc Q T G", nonzero=True, real=True)
    angle = zc*g/(eps*z0*q*sigma*(2*om+sigma))
    measured_spin = rho*mu*m*s.pi*eps*zc*t/sigma
    measured_mass = s.simplify(measured_spin/(sigma*angle))
    expected_mass = rho*mu*m*s.pi*eps**2*z0*q*t*(2*om+sigma)/(sigma*g)
    checks.check("signed radial tag inertia follows from physical shape and spin",
                 s.simplify(measured_mass-expected_mass) == 0)
    actual_mass = canonical_mass.subs(row, angle)
    target = s.solve(measured_mass-actual_mass, t)[0]/q
    expected_target = -beta*z0*sigma**2*(2*om+sigma)/(rho*mu*m*s.pi*zc**2*g)
    checks.check("moment matching ratio is solved from mechanical and KKS masses",
                 s.simplify(target-expected_target) == 0)
    checks.check("the derived ratio exactly matches both masses",
                 s.simplify((measured_mass-actual_mass).subs(t, q*target)) == 0)
    p = s.Symbol("p", positive=True)
    endpoint = s.factor(measured_mass.subs({sigma: -om, t: p*q}))
    checks.check("endpoint signed lobes still give a negative square",
                 s.simplify(endpoint+rho*mu*m*s.pi*eps**2*z0*p*q*q/g) == 0)
    q1, q2, t1, t2 = s.symbols("Q1 Q2 T1 T2", real=True)
    moment_matrix = s.Matrix([[q1, q2], [t1, t2]])
    prescribed = s.Matrix([q, q*target])
    coefficients = moment_matrix.inv()*prescribed
    checks.check("two radial bump coefficients prescribe both actual moment rows",
                 (moment_matrix*coefficients-prescribed).applyfunc(s.simplify) == s.zeros(2, 1))
    ell, r1, r2 = s.symbols("lambda r1 r2", positive=True)
    ratio1 = p*(1-ell**2*r1**2/(2*(m+1)))
    ratio2 = p*(1-ell**2*r2**2/(2*(m+1)))
    checks.check("radial moment independence has explicit second-order conditioning",
                 s.simplify(q1*q2*(ratio2-ratio1)
                            -p*q1*q2*ell**2*(r1*r1-r2*r2)/(2*(m+1))) == 0)

    # Exact frequency-dependent elimination, not a static stiffness inverse.
    nu, nu0, stiff, grad, coupling, wave = s.symbols("nu nu0 L a c k", real=True)
    pencil = s.Matrix([[nu0**2-nu**2+grad*wave**2, coupling*wave],
                       [coupling*wave, stiff-nu**2]])
    schur = pencil[0, 0]-pencil[0, 1]*pencil[1, 0]/pencil[1, 1]
    checks.check("dynamical Schur complement keeps the full fast denominator",
                 s.simplify(pencil.det()-(stiff-nu**2)*schur) == 0)
    optical_slope = -s.diff(schur, wave, 2).subs({wave: 0, nu: nu0})/2 \
        /s.diff(schur, nu).subs({wave: 0, nu: nu0})
    checks.check("actual optical shift includes the frequency-dependent reaction",
                 s.simplify(optical_slope-(grad-coupling**2/(stiff-nu0**2))/(2*nu0)) == 0)
    checks.check("a positive retained gradient can coexist with negative actual curvature",
                 optical_slope.subs({nu0: 1, stiff: 2, grad: 1, coupling: 2}) == -s.Rational(3, 2))
    checks.check("replacing the dynamical denominator by static stiffness is exposed",
                 s.simplify(optical_slope-(grad-coupling**2/stiff)/(2*nu0)) != 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
