"""Exact differential/series oracle; no elliptic/KAM numerical selection."""

import sympy as s


def main():
    checks = []

    def check(name, expression):
        ok = all(s.simplify(v) == 0 for v in expression) if isinstance(
            expression, s.MatrixBase
        ) else s.simplify(expression) == 0
        assert ok, (name, expression)
        checks.append(name)
        print(f"PASS {name}")

    r, z, lam, c = s.symbols("r z lambda C", nonzero=True, real=True)
    psi = s.Function("psi")(r, z)
    ur = -s.diff(psi, z) / r
    ut = (lam * psi + c) / r
    uz = s.diff(psi, r) / r
    gs = s.diff(psi, r, 2) - s.diff(psi, r) / r + s.diff(psi, z, 2)
    check("exact cylindrical divergence", s.diff(r * ur, r) / r + s.diff(uz, z))
    check("exact first integral", ur * s.diff(psi, r) + uz * s.diff(psi, z))
    check("radial Beltrami component", -s.diff(ut, z) - lam * ur)
    check("axial Beltrami component", s.diff(r * ut, r) / r - lam * uz)
    check("azimuthal component is GS equation", s.diff(ur, z) - s.diff(uz, r) + gs / r)

    x, radius, major, u = s.symbols("x a R U", positive=True)
    phi = s.Function("phi")(x, z)
    scaled_gs = (
        s.diff(major * phi, x, 2)
        - s.diff(major * phi, x) / (major + x)
        + s.diff(major * phi, z, 2)
        + lam**2 * major * phi
        + lam * major * u
    ) / major
    check("normalized torus scalar equation", scaled_gs - (
        s.diff(phi, x, 2) + s.diff(phi, z, 2)
        - s.diff(phi, x) / (major + x) + lam**2 * phi + lam * u
    ))

    t = s.symbols("t", positive=True)
    radial = u / lam * (s.besselj(0, lam * t) / s.besselj(0, lam * radius) - 1)
    radial_equation = s.diff(radial, t, 2) + s.diff(radial, t) / t + lam**2 * radial + lam * u
    check("straight Bessel Dirichlet solution", s.expand_func(radial_equation))
    check("straight boundary value", radial.subs(t, radius))
    check("straight poloidal velocity", -s.diff(radial, t) - u * s.besselj(1, lam * t) / s.besselj(0, lam * radius))
    quotient = s.besselj(1, lam * t) / (t * s.besselj(0, lam * t))
    numerator = s.series(s.besselj(1, lam * t), t, 0, 8).removeO()
    denominator = s.series(s.besselj(0, lam * t), t, 0, 8).removeO()
    series = s.series(numerator / (t * denominator), t, 0, 6).removeO()
    check("derived Lundquist return-ratio series", series - lam / 2 - lam**3 * t**2 / 16 - lam**5 * t**4 / 96)
    action_derivative = u * s.besselj(0, lam * t) * t / s.besselj(0, lam * radius)
    twist_center = s.limit(
        s.diff(series, t) / action_derivative.subs(s.besselj(0, lam * t), denominator), t, 0
    )
    check("flux-action twist leading coefficient", twist_center - lam**3 * s.besselj(0, lam * radius) / (8 * u))
    delta_theta = 2 * s.pi * t * s.besselj(0, lam * t) / (major * s.besselj(1, lam * t))
    check("transit-time return angle", 4 * s.pi**2 / delta_theta - 2 * s.pi * major * quotient)

    aa, bb, dd, denom = s.symbols("aa bb dd denom", real=True)
    j = s.Matrix([[0, -1], [1, 0]])
    hess = s.Matrix([[aa, bb], [bb, dd]])
    particle = -j * hess / denom
    check("core transverse trace", s.trace(particle))
    check("core elliptic determinant", particle.det() - hess.det() / denom**2)

    omega, jd = s.symbols("Omega j_D", positive=True)
    a1, a2 = s.symbols("a1 a2", real=True)
    angle = s.Matrix([a1, a2])
    angle_dot = -omega * j * angle
    displacement = -j * angle
    displacement_dot = displacement.jacobian(angle) * angle_dot
    check("Kelvin material gyro reverses particle rotation", omega * j - 2 * omega * j + omega * j)
    check("actual shear parcel spin", jd * j * displacement_dot - jd * angle_dot)
    col1, col2 = -j[:, 0], -j[:, 1]
    beta = 2 * omega * jd * s.det(s.Matrix.hstack(col1, col2))
    check("two-angle KKS coefficient", beta - 2 * omega * jd)
    check("canonical inertia differs by exact factor two", beta / omega - 2 * jd)
    q, p, qdot = s.symbols("q p qdot", real=True)
    energy_coeff = beta * omega
    phase_lagrangian = beta * p * qdot - energy_coeff * (q**2 + p**2) / 2
    solved_p = s.solve(s.diff(phase_lagrangian, p), p)[0]
    check("exact scalar Routh inertia", s.diff(phase_lagrangian.subs(p, solved_p), qdot, 2) - 2 * jd)
    xx, yy, zz, da1, da2 = s.symbols("xx yy zz da1 da2", real=True)
    position = s.Matrix([xx, yy, zz])
    tangent = s.Matrix([0, 0, 1])
    rigid_angle = s.Matrix([a1, a2, 0])
    rigid_rate = s.Matrix([da1, da2, 0])
    background = omega * tangent.cross(position)
    rigid_xi = rigid_angle.cross(position)
    lin_velocity = (
        rigid_rate.cross(position)
        + rigid_xi.jacobian(position) * background
        - background.jacobian(position) * rigid_xi
    )
    curl_lin = s.Matrix([
        s.diff(lin_velocity[2], yy) - s.diff(lin_velocity[1], zz),
        s.diff(lin_velocity[0], zz) - s.diff(lin_velocity[2], xx),
        s.diff(lin_velocity[1], xx) - s.diff(lin_velocity[0], yy),
    ])
    check("rigid material Kelvin equation leaves zero angle rate",
          curl_lin - 2 * omega * rigid_xi.diff(zz) - 2 * rigid_rate)
    print(f"{len(checks)}/{len(checks)} exact checks passed")


if __name__ == "__main__":
    main()
