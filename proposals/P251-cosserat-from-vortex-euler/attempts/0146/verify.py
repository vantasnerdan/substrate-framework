"""Exact full-Euler mean/Hodge/normal-form identities; no spectral sampling."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0146-acoustic-normal-form")
    x, y = s.symbols("x y", real=True)
    k = s.Symbol("k", nonzero=True, real=True)
    coords = (x, y)
    psi = s.cos(x)*s.cos(y)
    v = s.Matrix([-s.diff(psi, y), s.diff(psi, x)])
    pressure = -(s.cos(2*x)+s.cos(2*y))/4
    dv = v.jacobian(coords)
    av = dv*v
    gradp = s.Matrix([s.diff(pressure, q) for q in coords])
    checks.check("smooth fixture satisfies stationary Euler with physical pressure",
                 all(s.trigsimp(q) == 0 for q in av+gradp))
    c1, c2 = s.symbols("c1 c2", real=True)
    c = s.Matrix([c1, c2])
    translation = -dv*c
    pressure_t = -(c.dot(gradp))
    translated_euler = translation.jacobian(coords)*v+dv*translation \
        +s.Matrix([s.diff(pressure_t, q) for q in coords])
    checks.check("translation kernel follows the differentiated full Euler equation",
                 all(s.trigsimp(q) == 0 for q in translated_euler))
    poisson_source = 2*sum(s.diff(v[j], coords[i])*s.diff(translation[i], coords[j])
                           for i in range(2) for j in range(2))
    checks.check("translation pressure is the full periodic Poisson solution",
                 s.trigsimp(-sum(s.diff(pressure_t, q, 2) for q in coords)
                            -poisson_source) == 0)
    checks.check("translation pressure retains its transport derivative",
                 s.trigsimp(pressure_t-c.dot(av)) == 0)

    wz = v[0]
    wh = s.I*k*s.Matrix([s.diff(wz, q)/2 for q in coords])
    checks.check("nonzero-divergence Hodge return satisfies the axial constraint",
                 s.simplify(sum(s.diff(wh[j], coords[j]) for j in range(2))+s.I*k*wz) == 0)
    horizontal_rhs = -wh.jacobian(coords)*v-dv*wh
    mean_rhs = horizontal_rhs.applyfunc(
        lambda f: s.integrate(s.integrate(s.expand_trig(f), (x, -s.pi, s.pi)),
                              (y, -s.pi, s.pi))/(4*s.pi**2))
    mean_flux = (v*wz).applyfunc(
        lambda f: s.integrate(s.integrate(f, (x, -s.pi, s.pi)),
                              (y, -s.pi, s.pi))/(4*s.pi**2))
    checks.check("actual Euler mean retains the nonzero axial momentum flux",
                 mean_rhs == -s.I*k*mean_flux and mean_flux[0] == s.Rational(1, 4))

    aa, vv, xx, mm, zz, residual = s.symbols("Av v X m Z pi_r")
    awz = s.Symbol("A_wz")
    # Differentiate Z=w_z+i k v.X in time and apply the stationary A.
    z_lhs = -awz-s.I*k*(aa*xx+residual)+s.I*k*vv*mm+awz+s.I*k*aa*xx
    checks.check("compensated vertical variable keeps the cell pressure remainder",
                 s.expand(z_lhs-s.I*k*(vv*mm-residual)) == 0)
    cv, vz = s.symbols("C_v VZ")
    physical_mean_rhs = -s.I*k*(vz-s.I*k*cv*xx)
    checks.check("mean restoring sign is derived from the exact physical flux",
                 s.expand(physical_mean_rhs+k*k*cv*xx+s.I*k*vz) == 0)
    rzdot, crv, rpi = s.symbols("rZdot C_rv rpi_r")
    paired_vz = rzdot-s.I*k*crv*mm+s.I*k*rpi
    corrected_rhs = physical_mean_rhs.subs(vz, paired_vz)+s.I*k*rzdot
    target = -k*k*cv*xx-k*k*crv*mm+k*k*rpi
    checks.check("actual current normal form retains gyro and pressure rows",
                 s.expand(corrected_rhs-target) == 0)
    checks.check("omitting pressure is exposed by a nonzero mutation residual",
                 s.expand(corrected_rhs-target.subs(rpi, 0)) == k*k*rpi)

    t, omega = s.symbols("t omega", positive=True)
    rotation = s.Matrix([[s.cos(omega*t), -s.sin(omega*t)],
                         [s.sin(omega*t), s.cos(omega*t)]])
    cr = rotation
    rv = -s.diff(cr, t, 2)
    tau = s.Symbol("tau", real=True)
    memory = rv.applyfunc(lambda f: s.integrate((t-tau)*f.subs(t, tau), (tau, 0, t)))
    primitive = cr.subs(t, 0)-cr+t*s.diff(cr, t).subs(t, 0)
    checks.check("twice-integrated exact transport correlation has a bounded primitive",
                 all(s.simplify(q) == 0 for q in memory-primitive))
    initial_slope = s.diff(cr, t).subs(t, 0)
    checks.check("the remaining linear correlation term is antisymmetric",
                 initial_slope+initial_slope.T == s.zeros(2))
    paired_memory = (memory+memory.subs(omega, -omega))/2
    checks.check("actual sign-paired transport removes only the odd linear part",
                 all(s.simplify(q) == 0 for q in paired_memory
                     -s.eye(2)*(1-s.cos(omega*t))))

    lap, transport, qop, dop, fop, top, vop = s.symbols("s A Q d F T V")
    exact_matrix = s.Matrix([[lap, 0, s.I*k*vop],
                             [-top, lap, -s.I*k*fop],
                             [0, s.I*k*qop, lap+transport-k*k*qop*dop]])
    wop = lap+transport-k*k*qop*(dop+fop/lap)
    checks.check("full scalar block Schur retains both horizontal pressure couplings",
                 s.factor(exact_matrix.det()-(lap*lap*wop+k*k*vop*qop*top)) == 0)
    checks.check("transport resolvent uses the primitive rather than a spectral gap",
                 s.factor(transport/(lap+transport)
                          -(1-lap/(lap+transport))) == 0)
    rr = s.Symbol("r")
    checks.check("bounded primitive resolves the small-frequency transport row",
                 s.factor(transport*rr/(lap+transport)
                          -(rr-lap*rr/(lap+transport))) == 0)
    epsilon, speed, initial = s.symbols("epsilon speed V0", positive=True)
    leading_m = initial*s.cos(epsilon*speed*t)
    leading_x = initial*s.sin(epsilon*speed*t)/speed
    checks.check("the leading acoustic flow solves the rescaled physical mean system",
                 s.simplify(s.diff(leading_x, t)-epsilon*leading_m) == 0
                 and s.simplify(s.diff(leading_m, t)+epsilon*speed**2*leading_x) == 0)
    h1, h2, a1, a2, b1, b2 = s.symbols("h1 h2 a1 a2 b1 b2", real=True)
    h = s.diag(h1, h2)
    j = s.Matrix([[0, -1], [1, 0]])
    f = s.Matrix([a1, a2])
    forcing = s.Matrix([b1, b2])
    energy_rate = s.expand(2*(h*f).dot(j*h*f+forcing))
    checks.check("positive core energy retains actual ambient forcing work",
                 s.expand(energy_rate-2*(h*f).dot(forcing)) == 0
                 and energy_rate != 0)
    print("Scope: exact calculus/block identities; no bounded full-cell group or acoustic pole is inferred.")
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
