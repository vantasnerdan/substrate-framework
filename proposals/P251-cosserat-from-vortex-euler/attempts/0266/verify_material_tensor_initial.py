"""Integrate the full Cartesian material tensor and isotropic contraction."""

import math

import sympy as sp


def circle_integral(expression, cosine, sine):
    """Exact uniform angular integral of a polynomial in cos and sin."""
    result = 0
    for (cos_power, sin_power), coefficient in sp.Poly(
        sp.expand(expression), cosine, sine
    ).terms():
        if cos_power % 2 or sin_power % 2:
            continue
        a, b = cos_power//2, sin_power//2
        moment = 2*sp.pi*sp.Rational(
            math.factorial(2*a)*math.factorial(2*b),
            4**(a+b)*math.factorial(a)*math.factorial(b)*math.factorial(a+b),
        )
        result += coefficient*moment
    return sp.expand(result)


def main():
    ct, st, cp, snp = sp.symbols("ct st cp snp", real=True)
    radius, minor, r = sp.symbols("R s r", positive=True)
    speed, swirl, coefficient, gamma = sp.symbols("V W C gamma", real=True)
    er, ep, ez = sp.Matrix([cp, snp, 0]), sp.Matrix([-snp, cp, 0]), sp.Matrix([0, 0, 1])
    et = -st*er-ct*ez
    position = r*er-minor*st*ez
    velocity = speed*et+swirl*ep
    wa = swirl-gamma*r  # m=+1; gamma is not fitted to a result.
    displacement = coefficient*(speed*et+wa*ep)
    material_velocity = coefficient*(
        speed**2/minor*(-ct*er+st*ez)+wa/r*(-speed*st*ep-swirl*er)
    )

    tensor = {}
    for i in range(3):
        for j in range(3):
            for ell in range(j, 3):
                integrand = (
                    material_velocity[i]*position[j]*position[ell]
                    +velocity[i]*(displacement[j]*position[ell]+position[j]*displacement[ell])
                )
                phi_integral = circle_integral(integrand, cp, snp)
                theta_integrand = sp.cancel(phi_integral.subs(r, radius+minor*ct))*(ct+sp.I*st)
                value = sp.simplify(circle_integral(theta_integrand, ct, st))
                tensor[i, j, ell] = tensor[i, ell, j] = value
    parallel = sp.I*sp.pi**2*coefficient*speed**2*minor/2
    perpendicular = sp.I*sp.pi**2*coefficient*speed**2*(radius**2/minor+3*minor/4)
    assert tensor[2, 2, 2] == parallel
    assert tensor[2, 0, 0] == tensor[2, 1, 1] == perpendicular
    assert tensor[0, 2, 0] == tensor[1, 2, 1] == -parallel/2
    assert tensor[1, 2, 0] == -sp.pi**2*coefficient*speed*swirl*radius

    kx, ky, kz = sp.symbols("kx ky kz", real=True)
    wave = sp.Matrix([kx, ky, kz])
    full_output = sp.Matrix([
        sum(tensor[i, j, ell]*wave[j]*wave[ell] for j in range(3) for ell in range(3))
        for i in range(3)
    ])
    # The parity-even part is derived by reversing the swirl in the full tensor.
    even_output = sp.simplify((full_output+full_output.subs(swirl, -swirl))/2)
    expected = sp.Matrix([-parallel*kx*kz, -parallel*ky*kz,
                          perpendicular*(kx**2+ky**2)+parallel*kz**2])
    assert sp.simplify(even_output-expected) == sp.zeros(3, 1)
    transverse = sp.eye(3)-wave*wave.T/wave.dot(wave)
    assert sp.simplify(
        transverse*even_output
        -transverse*ez*(perpendicular*(kx**2+ky**2)+2*parallel*kz**2)
    ) == sp.zeros(3, 1)

    # Derive sphere moments by integrating azimuth and mu=cos(polar angle).
    mu = sp.symbols("mu", real=True)
    nx_squared_phi_integral = circle_integral((1-mu**2)*cp**2, cp, snp)
    moment2 = sp.integrate(nx_squared_phi_integral, (mu, -1, 1))/(4*sp.pi)
    moment4 = sp.integrate(nx_squared_phi_integral*mu**2, (mu, -1, 1))/(4*sp.pi)
    assert moment2 == sp.Rational(1, 3)
    assert moment4 == sp.Rational(1, 15)
    averaged = sp.expand(perpendicular*(moment2-moment4)+2*parallel*moment4)
    assert sp.simplify(averaged-(4*perpendicular+2*parallel)/15) == 0
    # Deleting the oblique mixed entries changes the measured transverse row.
    axial_only = sp.Matrix([0, 0, expected[2]])
    assert sp.simplify(transverse*(even_output-axial_only)) != sp.zeros(3, 1)
    print("Full Cartesian material tensor (all 18 symmetric entries), oblique "
          "projection, swirl parity, sphere moments, and mixed-row mutation: PASS")
    print("b_parallel =", parallel)
    print("b_perp =", perpendicular)
    print("transverse whole-law coefficient =", sp.factor(averaged))


if __name__ == "__main__":
    main()
