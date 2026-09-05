"""Actual nonuniform Euler operator and full-pressure polynomial response jet.

The polynomial pressure is the unique polynomial Helmholtz inverse in
each nonzero axial Fourier sector. This is a local comparison algebra,
not a finite-energy replacement for the separate Bessel cylinder.
"""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0133-nonuniform-intrinsic-response")
    x, y, t = s.symbols("x y t", real=True)
    omega, carrier, probe = s.symbols("Omega N M", positive=True)
    zeta, conjugate = x+s.I*y, x-s.I*y
    plus, minus, axis = s.Matrix([1, s.I, 0]), s.Matrix([1, -s.I, 0]), s.Matrix([0, 0, 1])
    speed = 2*omega/carrier

    def gradient(value, axial):
        return s.Matrix([s.diff(value, x), s.diff(value, y), s.I*axial*value])

    def curl(vector, axial):
        return s.Matrix([s.diff(vector[2], y)-s.I*axial*vector[1],
                         s.I*axial*vector[0]-s.diff(vector[2], x),
                         s.diff(vector[1], x)-s.diff(vector[0], y)])

    def project(vector, axial):
        divergence = s.expand(s.diff(vector[0], x)+s.diff(vector[1], y)+s.I*axial*vector[2])
        pressure = s.S.Zero
        if divergence != 0:
            term = divergence
            for order in range(s.Poly(divergence, x, y).total_degree()//2+1):
                pressure -= term/axial**(2*order+2)
                term = s.expand(s.diff(term, x, 2)+s.diff(term, y, 2))
            assert term == 0
        return s.simplify(vector-gradient(pressure, axial))

    def base_generator(vector, axial):
        return project(speed*axis.cross(curl(vector, axial)-carrier*vector), axial)

    vplus, vminus = zeta**2*plus, conjugate**2*minus
    checks.check("finite-amplitude background harmonics are on the same curl shell",
                 s.simplify(curl(vplus, carrier)-carrier*vplus) == s.zeros(3, 1)
                 and s.simplify(curl(vminus, -carrier)-carrier*vminus) == s.zeros(3, 1))
    checks.check("axial drift closes absolute-vorticity Beltrami identity",
                 s.simplify(carrier*speed-2*omega) == 0)
    checks.check("actual same-shell amplitude direction is neutral for the uniform part",
                 base_generator(vplus, carrier) == s.zeros(3, 1))

    # Full actual interaction is (1/2)P[V_± cross(curl-N)w].
    intermediate = conjugate*minus
    shifted = probe-carrier
    first = project(vminus.cross(curl(plus, probe)-carrier*plus)/2, shifted)
    checks.check("first nonuniform Euler scattering is pressure completed",
                 s.simplify(first+2*intermediate) == s.zeros(3, 1))
    checks.check("the opposite first scattering vanishes by actual polarization",
                 s.simplify(vplus.cross(curl(plus, probe)-carrier*plus)) == s.zeros(3, 1))
    return_field = project(vplus.cross(curl(intermediate, shifted)-carrier*intermediate)/2, probe)
    basis = [plus,
             (x*x+y*y)*plus+2*s.I*zeta*axis/probe,
             zeta**2*minus+4*s.I*zeta*axis/probe]
    basis_matrix = s.Matrix.hstack(*basis)
    checks.check("second scattering retains its pressure-induced core velocity and shape",
                 s.simplify(return_field+basis_matrix*s.Matrix([8/probe**2, 2, 1]))
                 == s.zeros(3, 1))
    fundamental = speed*probe-2*omega
    opposite = speed*probe+2*omega
    middle = speed*probe
    matrix = s.Matrix([[-s.I*fundamental, 4*s.I*omega/probe**2, -8*s.I*omega/probe**2],
                       [0, -s.I*fundamental, 0], [0, 0, -s.I*opposite]])
    for index, vector in enumerate(basis):
        checks.check(f"full intrinsic Euler generator on return column {index}",
                     s.simplify(base_generator(vector, probe)-basis_matrix*matrix[:, index])
                     == s.zeros(3, 1))
    checks.check("intermediate physical field has its derived intrinsic frequency",
                 s.simplify(base_generator(intermediate, shifted)+s.I*middle*intermediate)
                 == s.zeros(3, 1))
    observation = s.Matrix([[s.simplify(curl(vector, probe)[0].subs({x: 0, y: 0})/probe)
                             for vector in basis]])
    checks.check("actual core vorticity row includes the pressure/shape cancellations",
                 observation == s.Matrix([[1, -2/probe**2, -4/probe**2]]))
    checks.check("return has zero instantaneous core-vorticity observation despite nonzero velocity",
                 s.simplify(observation*s.Matrix([-8/probe**2, -2, -1])) == s.zeros(1)
                 and return_field.subs({x: 0, y: 0}) != s.zeros(3, 1))

    # Ordered Duhamel coefficient, derived from those Euler matrices.
    laplace = s.Symbol("p")
    initial = s.Matrix([1, 0, 0])/(laplace+s.I*fundamental)
    hidden = -2*initial[0]/(laplace+s.I*middle)
    response = (laplace*s.eye(3)-matrix).inv()*s.Matrix([-8/probe**2, -2, -1])*hidden
    measured = s.simplify((observation*response)[0])
    f = laplace+s.I*fundamental
    g = laplace+s.I*opposite
    checks.check("actual measured second-amplitude transfer retains the resonant return",
                 s.simplify(measured-32*s.I*omega/(probe**2*f**3*g)) == 0)
    correction = s.exp(-s.I*fundamental*t)/probe**2*(
        4*t*t+2*s.I*t/omega-1/(2*omega**2)+s.exp(-4*s.I*omega*t)/(2*omega**2))
    checks.check("physical response correction vanishes through the second time derivative",
                 all(s.simplify(s.diff(correction, t, order).subs(t, 0)) == 0
                     for order in range(3)))
    checks.check("first genuine observed nonuniform correction is the third time derivative",
                 s.simplify(s.diff(correction, t, 3).subs(t, 0)-32*s.I*omega/probe**2) == 0)
    material_transport = s.exp(s.I*(speed*probe-omega)*t)
    checks.check("actual axis material transport and vector rotation recover laboratory tilt",
                 s.simplify(material_transport*s.exp(-s.I*fundamental*t)-s.exp(s.I*omega*t)) == 0)

    # At finite background amplitude the same Kelvin generator prepares a
    # different Eulerian field. Keep that exact order-a initial correction.
    kelvin_generator = -s.I*plus/(2*omega)
    prepared = project(kelvin_generator.cross(carrier*vminus)/2, shifted)
    prepared_coefficient = -s.I*carrier/(omega*shifted)
    checks.check("fixed actual Kelvin generator retains its nonuniform initial pressure field",
                 s.simplify(prepared-prepared_coefficient*intermediate) == s.zeros(3, 1))
    prepared_hidden = (prepared_coefficient-2/(laplace+s.I*fundamental))/(laplace+s.I*middle)
    prepared_return = (laplace*s.eye(3)-matrix).inv()*s.Matrix([-8/probe**2, -2, -1])*prepared_hidden
    prepared_observed = s.simplify((observation*prepared_return)[0])
    expected_prepared = (32*s.I*omega-16*carrier*f/shifted)/(probe**2*f**3*g)
    checks.check("actual Kelvin-prepared physical response differs from fixed Eulerian data",
                 s.simplify(prepared_observed-expected_prepared) == 0)
    checks.check("Kelvin preparation changes the second initial physical-angle derivative",
                 s.simplify(s.limit(laplace**3*prepared_observed, laplace, s.oo)
                            +16*carrier/(shifted*probe**2)) == 0)
    print("Order-a^2 Laplace core-vorticity transfer:", s.factor(measured))
    print("Order-a^2 third initial time derivative:", 32*s.I*omega/probe**2)
    print("Kelvin-prepared order-a^2 second initial derivative:",
          -16*carrier/(shifted*probe**2))
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
