"""Exact full-current/boundary-action map without changing actual U momentum."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0235-current-boundary-map")
    t, x, y, z = s.symbols("t x y z", real=True)
    xyz = (x, y, z)
    q = s.Function("q")(t)
    u = s.Matrix([t*x*y+t**2*z, t*y*z+t**3*x, t*z*x+t**2*y])
    phi = s.Matrix([x**2+t*y*z, y**2+t*x*z, z**2+t*x*y])

    def div_matrix(value):
        return s.Matrix([sum(s.diff(value[i, j], xyz[j]) for j in range(3))
                         for i in range(3)])

    def pair(left, right):
        return sum(left[i, j]*right[i, j] for i in range(3) for j in range(3))

    current = s.Matrix(3, 3, lambda i, j:
                       q*sum(s.LeviCivita(i, j, k)*s.diff(u[k], t) for k in range(3)))
    spin_difference = div_matrix(current)
    gradient = phi.jacobian(xyz)
    delta_l = -spin_difference.dot(phi.diff(t))+pair(current.diff(t), gradient)
    temporal = pair(current, gradient)
    spatial = current.T*phi.diff(t)
    boundary = s.diff(temporal, t)-sum(s.diff(spatial[j], xyz[j]) for j in range(3))
    checks.check("the complete observer difference is an exact spacetime action boundary",
                 s.simplify(delta_l-boundary) == 0)
    checks.check("spin and couple-current improvement preserve the actual angular balance",
                 s.simplify(-spin_difference.diff(t)+div_matrix(current.diff(t))) == s.zeros(3, 1))
    checks.check("changing the spin alone would discard a nonzero physical memory current",
                 s.simplify(spin_difference.diff(t)) != s.zeros(3, 1))

    # Independent jets expose every term of the translation Noether momentum.
    qc, qdot = s.symbols("qc qdot", real=True)
    velocity = s.Matrix(s.symbols("v0:3", real=True))
    acceleration = s.Matrix(s.symbols("a0:3", real=True))
    velocity_gradient = s.Matrix(3, 3, s.symbols("hv0:9", real=True))
    rotation_rate = s.Matrix(s.symbols("w0:3", real=True))
    rotation_gradient = s.Matrix(3, 3, s.symbols("g0:9", real=True))
    rate_gradient = s.Matrix(3, 3, s.symbols("gw0:9", real=True))
    density = sum(s.LeviCivita(i, j, k)*(
        -qc*velocity_gradient[k, j]*rotation_rate[i]
        +(qdot*velocity[k]+qc*acceleration[k])*rotation_gradient[i, j])
        for i in range(3) for j in range(3) for k in range(3))
    charge = []
    without_gradient_return = []
    for k in range(3):
        first = s.diff(density, velocity[k])
        time_return = sum(s.LeviCivita(i, j, k)*(
            qdot*rotation_gradient[i, j]+qc*rate_gradient[i, j])
            for i in range(3) for j in range(3))
        spatial_return = -qc*sum(s.LeviCivita(i, j, k)*rate_gradient[i, j]
                                for i in range(3) for j in range(3))
        charge.append(s.expand(first-time_return-spatial_return))
        without_gradient_return.append(s.expand(first-time_return))
    checks.check("the complete higher-derivative boundary term changes no translation momentum",
                 charge == [0, 0, 0])
    checks.check("omitting the gradient-velocity charge would falsely change physical U",
                 without_gradient_return != [0, 0, 0])
    expected_spin = s.Matrix([-qc*sum(s.LeviCivita(i, j, k)*velocity_gradient[k, j]
                                     for j in range(3) for k in range(3)) for i in range(3)])
    checks.check("the same action gives the derived internal-spin improvement",
                 s.simplify(s.Matrix([s.diff(density, w) for w in rotation_rate])
                            -expected_spin) == s.zeros(3, 1))
    couple = -s.Matrix(3, 3, lambda i, j: s.diff(density, rotation_gradient[i, j]))
    expected_couple = -s.Matrix(3, 3, lambda i, j: sum(s.LeviCivita(i, j, k)*(
        qdot*velocity[k]+qc*acceleration[k]) for k in range(3)))
    checks.check("the same action retains the full time-dependent couple-flux correction",
                 s.simplify(couple-expected_couple) == s.zeros(3))
    j0, ix, s0, trace = s.symbols("j0 IX S0 trace", nonzero=True, real=True)
    c = (j0-ix)/s0
    difference = (trace+c*s0)/3-(ix+c*s0)/3
    checks.check("the observer current difference is independent of the target optical normalization",
                 s.simplify(difference-(trace-ix)/3) == 0
                 and s.diff(s.simplify(difference), j0) == 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
