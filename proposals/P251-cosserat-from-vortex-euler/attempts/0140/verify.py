"""Full force-free KKS and exact time-dependent material-angle action."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0140-forcefree-material-action")
    f, b, c, vort, axial_derivative = s.symbols("f B C Z Wprime", real=True)
    real_eta = s.Matrix([f, 0, 0])
    imag_eta = s.Matrix([0, b, c])
    background_vorticity = s.Matrix([0, -axial_derivative, vort])
    kks = background_vorticity.dot(real_eta.cross(imag_eta))
    checks.check("force-free KKS includes the azimuthal vorticity contribution",
                 s.expand(kks-f*(vort*b+axial_derivative*c)) == 0)
    checks.check("discarding axial-shear vorticity changes the exact action",
                 s.expand(kks-f*vort*b) != 0)

    om, radius, mode, pressure, sigma = s.symbols("O r m p s", nonzero=True, real=True)
    vt = mode*pressure/(radius*sigma)-2*om*f
    checks.check("same full material spin includes deformation and velocity",
                 s.simplify(radius*(2*om*f+vt)-mode*pressure/sigma) == 0)

    t = s.Symbol("t", real=True)
    beta, frequency = s.symbols("beta sigma_ref", nonzero=True, real=True)
    phase = s.Function("delta")(t)
    z1, z2 = s.Function("z1")(t), s.Function("z2")(t)
    rotation = s.Matrix([[s.cos(phase), -s.sin(phase)],
                         [s.sin(phase), s.cos(phase)]])
    z = s.Matrix([z1, z2])
    y = rotation.T*z
    cross_y = y[0]*s.diff(y[1], t)-y[1]*s.diff(y[0], t)
    cross_z = z1*s.diff(z2, t)-z2*s.diff(z1, t)
    checks.check("time-dependent physical row rotates the complete symplectic one-form",
                 s.trigsimp(cross_y-cross_z+s.diff(phase, t)*(z1*z1+z2*z2)) == 0)
    gamma = frequency+s.diff(phase, t)
    lagrangian = -beta*(z1*s.diff(z2, t)-z2*s.diff(z1, t))/2 \
        +beta*gamma*(z1*z1+z2*z2)/2
    # Subtract a total derivative to make z1 algebraic.
    algebraic = s.simplify(lagrangian-s.diff(beta*z1*z2/2, t))
    solved = s.solve(s.diff(algebraic, z1), z1)[0]
    reduced = s.simplify(algebraic.subs(z1, solved))
    checks.check("conjugate quadrature is eliminated from the same time-dependent action",
                 s.simplify(reduced+beta*s.diff(z2, t)**2/(2*gamma)
                            -beta*gamma*z2*z2/2) == 0)

    theta = s.Function("theta")(t)
    row = s.Function("c")(t)
    pulled_back = s.simplify(reduced.subs({s.diff(z2, t): s.diff(theta/row, t), z2: theta/row}))
    connection = s.diff(row, t)/row
    mass = -beta/(gamma*row**2)
    expected = mass*((s.diff(theta, t)-connection*theta)**2-gamma**2*theta**2)/2
    checks.check("physical scalar pullback retains both temporal connection terms",
                 s.simplify(pulled_back-expected) == 0)
    momentum = s.diff(expected, s.diff(theta, t))
    checks.check("physical canonical momentum is not silently replaced by M theta_dot",
                 s.simplify(momentum-mass*(s.diff(theta, t)-connection*theta)) == 0)
    psi = frequency*t+phase
    phase_solution = row*s.sin(psi)
    exact_momentum = s.simplify(momentum.subs({s.diff(theta, t): s.diff(phase_solution, t),
                                             theta: phase_solution}))
    checks.check("prepared physical phase has the exact momentum minus beta cos psi over c",
                 s.simplify(exact_momentum+beta*s.cos(psi)/row) == 0)

    rho, mu, azimuth, mark, z0, zc, q, a0, tag_t = s.symbols(
        "rho mu m eps Z0 Zc Q A0 Ttag", nonzero=True, real=True)
    row0 = zc*a0/(mark*z0*q)
    physical_spin0 = rho*mu*azimuth*s.pi*mark*zc*tag_t
    target = s.solve(physical_spin0+beta/row0, tag_t)[0]/q
    expected_target = -beta*z0/(rho*mu*azimuth*s.pi*zc**2*a0)
    checks.check("tag target is derived from full physical-action momentum",
                 s.simplify(target-expected_target) == 0)
    checks.check("actual initial spin exactly matches canonical momentum",
                 s.simplify((physical_spin0+beta/row0).subs(tag_t, q*target)) == 0)

    actual_doppler, radial_g = s.symbols("sigma G", nonzero=True, real=True)
    ordinary_target = s.simplify(target.subs(a0, radial_g/(actual_doppler*(2*om+actual_doppler)))
                                *actual_doppler)
    checks.check("ordinary-column limit agrees after converting the pressure-over-Doppler row",
                 s.simplify(ordinary_target+beta*z0*actual_doppler**2*(2*om+actual_doppler)
                            /(rho*mu*azimuth*s.pi*zc**2*radial_g)) == 0)

    k, drift = s.symbols("k h", real=True)
    ref = s.Function("sigma")(k)
    exact_phase = s.exp(s.I*(ref-k*drift)*t)
    reference_phase = s.exp(s.I*ref*t)
    for order in range(3):
        coefficient = s.simplify(s.diff(exact_phase-reference_phase, k, order).subs(drift, 0))
        checks.check(f"phase-error carrier derivative {order} vanishes at zero axial drift",
                     coefficient == 0)
    checks.check("carrier differentiation retains the axial drift source",
                 s.simplify(s.diff(exact_phase, k)-s.I*t*(s.diff(ref, k)-drift)*exact_phase) == 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
