"""Exact full-pressure one-wave response and stationary core continuation."""

import sympy as s

from substrate_framework.euler_fourier import (
    add, cross, curl, divergence, leray, scale, transport, trig,
)
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0151-whole-field-one-wave-second-jet")
    z, phi, c = s.symbols("z phi c", real=True)
    a = s.sqrt(1-z*z)
    normal = s.Matrix([0, 0, 1])
    direction = s.Matrix([a, 0, z])
    pn = s.eye(3)-normal*normal.T
    pk = s.eye(3)-direction*direction.T
    covariance = c*pn
    bare = s.simplify(direction.dot(covariance*direction))
    # Contract the actual initial Kelvin row, not its proposed eigenvalues.
    kelvin_force = pk*(
        bare*s.eye(3)-pn*(bare*s.eye(3)+direction*(direction.T*covariance))
        -covariance*direction*(direction.T*pn)
        -covariance*(direction.dot(pn*direction)))
    velocity_force = pk*(bare*s.eye(3)-2*bare*pn
                        -2*covariance*direction*(direction.T*pn))
    sh = s.Matrix([0, 1, 0])
    sv = s.Matrix([z, 0, -a])
    def clean(matrix):
        return matrix.applyfunc(s.simplify)

    checks.check("the two genuine initial phases give the same observed force",
                 clean((kelvin_force-velocity_force)*pk) == s.zeros(3))
    checks.check("the physical SH coefficient follows from covariance contraction",
                 clean(kelvin_force*sh+c*(1-z*z)*sh) == s.zeros(3, 1))
    checks.check("the full-pressure SV coefficient retains its negative sector",
                 clean(kelvin_force*sv+c*(1-z*z)*(4*z*z-1)*sv) == s.zeros(3, 1))

    # Exact actual Euler/Leray check with an oblique macro direction.
    k = s.Symbol("k", real=True, nonzero=True)
    background = (trig(2), trig(2, kind="sin"), {})
    omega = curl(background)
    checks.check("the circular field is actual constant-pressure stationary Euler",
                 not divergence(background)
                 and all(not component for component in transport(background, background)))

    def linear(field):
        left = transport(background, field)
        right = transport(field, background)
        return leray(tuple(scale(add(left[i], right[i]), -1) for i in range(3)))

    macro = (3*k/5, 0, 4*k/5)
    kk = s.Matrix([s.Rational(3, 5), 0, s.Rational(4, 5)])
    pp = s.eye(3)-kk*kk.T
    nn = pp*normal
    target = s.Rational(1, 2)*(1-kk[2]**2)*pp \
        -(1-2*kk[2]**2)*nn*nn.T

    def coefficient(field, power):
        return s.Matrix([s.simplify(s.diff(s.cancel(component.get(macro, 0)), k, power)
                                   .subs(k, 0)/s.factorial(power)) for component in field])

    for label, polarization in (("SH", sh), ("SV", s.Matrix([s.Rational(4, 5), 0,
                                                             -s.Rational(3, 5)]))):
        eta = tuple({macro: value} if value else {} for value in polarization)
        position = leray(cross(eta, omega))
        common = eta
        pos_force = linear(position)
        velocity_force_field = linear(linear(common))
        checks.check(f"canonical full Euler oblique {label} Kelvin force matches the derived matrix",
                     clean(coefficient(pos_force, 2)+target*polarization) == s.zeros(3, 1))
        checks.check(f"canonical full Euler oblique {label} common-V phase matches independently",
                     clean(coefficient(velocity_force_field, 2)+target*polarization) == s.zeros(3, 1))
        # The first later-time derivative tests the proposed whole-jet
        # reduction; the general proof uses the exact cell equation.
        checks.check(f"canonical next-time {label} mean derivative has no missed second spatial term",
                     coefficient(linear(pos_force), 2) == s.zeros(3, 1)
                     and coefficient(linear(velocity_force_field), 2) == s.zeros(3, 1))

    # Haar average in coordinates kappa=e_z, e=e_x.
    nx = s.sqrt(1-z*z)*s.cos(phi)
    projected_normal = s.Matrix([nx, s.sqrt(1-z*z)*s.sin(phi), 0])
    observed = c*(1-z*z)*s.Matrix([1, 0, 0]) \
        -2*c*(1-2*z*z)*projected_normal*nx

    def sphere(expression):
        return s.simplify(s.integrate(s.integrate(s.expand(expression),
                                                 (phi, 0, 2*s.pi)), (z, -1, 1))/(4*s.pi))

    averaged = s.Matrix([sphere(value) for value in observed])
    second_moment = sphere(observed.dot(observed))
    variance = s.factor(second_moment-averaged.dot(averaged))
    checks.check("whole-field Haar integration yields positive 2/15 of actual energy",
                 averaged == s.Matrix([4*c/15, 0, 0]))
    checks.check("orientation response variance is nonzero and positive",
                 s.simplify(variance/c**2).is_positive is True)
    print("Derived orientation second moment:", second_moment)
    print("Derived orientation variance:", variance)
    checks.check("the earlier perpendicular unstable branch is retained",
                 (c*(1-z*z)*(4*z*z-1)).subs(z, 0) == -c)

    # Actual small-amplitude elliptic tube, not an externally supplied mode.
    x, y, zz, d = s.symbols("x y zz d", real=True)
    coords = [x, y, zz]
    field = s.Matrix([s.cos(zz)+d*s.sin(y), s.sin(zz), d*s.cos(y)])
    actual_curl = s.Matrix([s.diff(field[(i+2) % 3], coords[(i+1) % 3])
                           -s.diff(field[(i+1) % 3], coords[(i+2) % 3]) for i in range(3)])
    pressure = -field.dot(field)/2
    checks.check("small insertion satisfies the exact Beltrami and Euler equations",
                 clean(actual_curl+field) == s.zeros(3, 1)
                 and s.simplify(sum(s.diff(field[i], coords[i]) for i in range(3))) == 0
                 and clean(field.jacobian(coords)*field+s.Matrix([s.diff(pressure, coord)
                                                                 for coord in coords])) == s.zeros(3, 1))
    h = -s.cos(zz)-d*s.sin(y)
    checks.check("the transverse first integral and actual axial speed agree",
                 s.simplify(sum(field[i]*s.diff(h, coords[i]) for i in range(3))) == 0
                 and s.simplify(field[0]+h) == 0)
    transverse = field[1:3, :].jacobian([y, zz]).subs({y: s.pi/2, zz: 0})
    checks.check("the periodic core has a genuine elliptic transverse linearization",
                 transverse == s.Matrix([[0, 1], [-d, 0]])
                 and transverse**2 == -d*s.eye(2))
    angle, action = s.symbols("angle action", real=True)
    quartic = -((s.sqrt(2*action)*s.cos(angle))**4
                +d*(s.sqrt(2*action)*s.sin(angle))**4)/24
    averaged_quartic = s.integrate(quartic, (angle, 0, 2*s.pi))/(2*s.pi)
    checks.check("the Birkhoff coefficient is derived from the actual Taylor Hamiltonian",
                 s.simplify(averaged_quartic+(1+d)*action**2/16) == 0)
    energy = s.Symbol("energy", real=True)
    frequency = s.sqrt(d)-(1+d)*energy/(8*s.sqrt(d))
    twist = s.factor(s.diff(frequency/(1+d-energy), energy).subs(energy, 0))
    checks.check("the return-map twist includes the varying actual axial speed",
                 s.simplify(twist-(8*d-(1+d)**2)/(8*s.sqrt(d)*(1+d)**2)) == 0)
    print("Derived periodic-core twist:", twist)
    print("Scope: actual fixed-time second jet, full retained action, periodic elliptic tube;")
    print("not acoustic-time homogenization or a completed EPS/optical action join.")
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
