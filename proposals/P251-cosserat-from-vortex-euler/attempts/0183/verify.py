"""Exact Cartesian Euler, helical tag and actual operator identities."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0183-stationary-helical-clock")
    x, y, z = s.symbols("x y z", real=True)
    c, cap, a = s.symbols("c C a", positive=True)
    xyz = s.Matrix([x, y, z])
    radius2 = x*x+y*y
    h = s.Matrix([-y, x, c])
    profile = s.Function("f")(radius2)
    u = profile*h
    jac = u.jacobian(xyz)
    curl = s.Matrix([s.diff(u[2], y)-s.diff(u[1], z),
                     s.diff(u[0], z)-s.diff(u[2], x),
                     s.diff(u[1], x)-s.diff(u[0], y)])
    checks.check("every radial profile is solenoidal", s.simplify(s.trace(jac)) == 0)
    checks.check("actual radial acceleration fixes localized Euler pressure",
                 s.simplify(jac*u+profile**2*s.Matrix([x, y, 0])) == s.zeros(3, 1))
    chosen = cap/(c*c+radius2)
    velocity = chosen*h
    pressure = -cap**2/(2*(c*c+radius2))
    checks.check("smooth finite-core profile solves Cartesian Euler",
                 s.simplify(velocity.jacobian(xyz)*velocity
                            +s.Matrix([s.diff(pressure, q) for q in xyz])) == s.zeros(3, 1))
    actual_curl = s.Matrix([s.diff(velocity[2], y)-s.diff(velocity[1], z),
                           s.diff(velocity[0], z)-s.diff(velocity[2], x),
                           s.diff(velocity[1], x)-s.diff(velocity[0], y)])
    checks.check("general radial curl retains the localized-profile derivative",
                 s.simplify(curl[0]-c*s.diff(profile, y)) == 0
                 and s.simplify(curl[1]+c*s.diff(profile, x)) == 0
                 and s.simplify(curl[2]-2*profile-x*s.diff(profile, x)
                                -y*s.diff(profile, y)) == 0)
    checks.check("variable curl factor is derived, not constant lambda",
                 s.simplify(actual_curl-2*c*velocity/(c*c+radius2)) == s.zeros(3, 1)
                 and s.diff(2*c/(c*c+radius2), x) != 0)
    checks.check("force-free pressure equals negative kinetic energy",
                 s.simplify(pressure+velocity.dot(velocity)/2) == 0)
    for m in (1, 2, 3):
        moment = (x+s.I*y)**m*s.exp(-s.I*m*z/c)
        grad = s.Matrix([s.diff(moment, q) for q in xyz])
        checks.check(f"regular helical moment m={m} is materially invariant",
                     s.simplify(h.dot(grad)) == 0)
        checks.check(f"moment m={m} has unit rotation and nonzero axial current",
                     s.simplify(-y*grad[0]+x*grad[1]-s.I*m*moment) == 0
                     and s.simplify(grad[2]+s.I*m*moment/c) == 0
                     and grad[2] != 0)
    theta, phase, r, b = s.symbols("theta phase r b", real=True)
    harmonic = r**2*s.exp(s.I*2*(theta-phase))
    tag = 1+b*r**2*s.cos(2*(theta-phase))/a**2
    angular = s.integrate(s.expand_complex(tag*harmonic), (theta, 0, 2*s.pi))/(2*s.pi)
    checks.check("positive tag has derived nonzero reference material moment",
                 s.simplify(angular-b*r**4/(2*a**2)) == 0)
    vx, vy, vz = s.symbols("v_x v_y v_z", real=True)
    v = s.Matrix([vx, vy, vz])
    rotation = s.Matrix([-vy, vx, 0])
    # Equivariance fixes (h.grad)v=Jv; retain the independent (v.grad)u.
    convection = profile*rotation+jac*v
    radial_derivative = s.diff(profile, x)*vx+s.diff(profile, y)*vy
    checks.check("helical reduction retains both basis rotation terms",
                 s.simplify(convection-2*profile*rotation-radial_derivative*h) == s.zeros(3, 1))
    checks.check("constant helical Killing momentum is a property of the actual profile",
                 s.simplify(velocity.dot(h)-cap) == 0
                 and s.simplify(u.dot(h)-profile*(radius2+c*c)) == 0)
    print("Pressure, variable curl, invariant tag and full helical convection derived.")
    print("Mode action, mechanical spin and Euclidean EPS closure are separate constructions.")
    raise SystemExit(checks.finish())


if __name__ == "__main__":
    main()
