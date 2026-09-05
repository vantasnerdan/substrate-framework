"""Exact finite-parcel cotangent, connection, Kelvin and reduced-spin identities."""

import sympy as s

from substrate_framework.euler_orbit import hermitian_schur_jet
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0078-finite-parcel-reduction")
    mass, inertia, metric, connection = s.symbols("M I G B", positive=True)
    omega, rate, momentum = s.symbols("Omega zdot pi_z", real=True)
    kinetic = inertia*omega**2/2+connection*omega*rate+metric*rate**2/2
    shape_solution = s.solve(s.diff(kinetic, rate)-momentum, rate)[0]
    routh = s.simplify((kinetic-momentum*rate).subs(rate, shape_solution))
    reduced_inertia = inertia-connection**2/metric
    expected = (reduced_inertia*omega**2/2+connection*momentum*omega/metric
                -momentum**2/(2*metric))
    ledger.check("shape-momentum Routh reduction retains the mechanical connection",
                 s.simplify(routh-expected) == 0)
    completed = (inertia*(omega+connection*rate/inertia)**2/2
                 +(metric-connection**2/inertia)*rate**2/2)
    ledger.check("geometric mechanical connection is the exact kinetic square completion",
                 s.simplify(kinetic-completed) == 0)

    # Six material samples corroborate the general exact continuum integral
    # identity; a genuine 3-D affine shape tangent has a NONZERO connection.
    points = [sign*s.eye(3)[:, axis] for axis in range(3) for sign in (-1, 1)]
    deformation = s.diag(1, 2, 3)
    shape = s.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    centered = [deformation*a for a in points]
    shapes = [shape*a for a in points]
    geometric = sum((h.dot(h)*s.eye(3)-h*h.T for h in centered), s.zeros(3))
    body_connection = sum((h.cross(v) for h, v in zip(centered, shapes, strict=True)),
                          s.zeros(3, 1))
    shape_metric = sum(v.dot(v) for v in shapes)
    ledger.check("explicit centered affine shape remains in the polar/Eckart slice",
                 sum((v*a.T for v, a in zip(shapes, points, strict=True)), s.zeros(3))
                 == 2*shape and s.trace(shape*deformation.inv()) == 0)
    ledger.check("actual geometric inertia and shape connection are both retained",
                 geometric[2, 2] == 10 and body_connection[2] == -2 and shape_metric == 4)
    ledger.check("computed example has positive reduced inertia different from the locked inertia",
                 geometric[2, 2]-body_connection[2]**2/shape_metric == 9)
    translation = s.Matrix(s.symbols("V0:3", real=True))
    body_rate = s.Matrix([0, 0, omega])
    velocities = [translation+body_rate.cross(h)+rate*v
                  for h, v in zip(centered, shapes, strict=True)]
    full_energy = sum(v.dot(v) for v in velocities)/2
    split = (len(points)*translation.dot(translation)/2
             +omega**2*geometric[2, 2]/2+omega*rate*body_connection[2]
             +rate**2*shape_metric/2)
    ledger.check("material mass-centroid kinetic split is exact without discarding affine motion",
                 s.simplify(full_energy-split) == 0)
    spin = sum((h.cross(v) for h, v in zip(centered, velocities, strict=True)),
               s.zeros(3, 1))
    ledger.check("physical spin equals locked inertia times frame rate plus shape momentum",
                 s.simplify(spin-geometric*body_rate-rate*body_connection) == s.zeros(3, 1))

    x, y, z = coords = s.symbols("x y z", real=True)
    e = s.Matrix(s.symbols("e0:3", real=True))
    position = s.Matrix(coords)

    def curl(v):
        return s.Matrix([s.diff(v[2], y)-s.diff(v[1], z),
                         s.diff(v[0], z)-s.diff(v[2], x),
                         s.diff(v[1], x)-s.diff(v[0], y)])

    ledger.check("adding geometric rigid velocity changes Kelvin vorticity by exactly 2Omega",
                 curl(e.cross(position)) == 2*e)
    angle = s.Symbol("theta", real=True)
    rotation = s.Matrix([[s.cos(angle), -s.sin(angle), 0],
                         [s.sin(angle), s.cos(angle), 0], [0, 0, 1]])
    ledger.check("material-label rigid rotation is volume preserving and can be a relabeling",
                 s.simplify(rotation.det()) == 1
                 and s.simplify(rotation.T*rotation) == s.eye(3))

    # Exact radial Helmholtz construction for the smooth finite Euler ball.
    radius, lam = s.symbols("r lambda", positive=True)
    f = s.sin(lam*radius)/(lam**2*radius**3)-s.cos(lam*radius)/(lam*radius**2)
    ode = s.diff(f, radius, 2)+4*s.diff(f, radius)/radius+lam**2*f
    ledger.check("explicit spherical-Bessel radial profile obeys its exact Helmholtz ODE",
                 s.simplify(ode) == 0)
    ledger.check("finite-ball Beltrami example is smooth with nonzero core velocity",
                 s.limit(f, radius, 0) == lam/3
                 and s.limit((2*f+radius*s.diff(f, radius))/lam, radius, 0) == s.Rational(2, 3))
    f0, fprime = s.symbols("f fprime", real=True)
    ez = s.Matrix([0, 0, 1])
    velocity = (f0*ez.cross(position)
                +((2*f0+radius*fprime)*ez-fprime*z*position/radius)/lam)
    radial_residual = s.together(position.dot(velocity)-2*f0*z/lam).as_numer_denom()[0]
    ledger.check("j1 boundary zero gives exact tangent velocity on the material sphere",
                 s.rem(radial_residual, x*x+y*y+z*z-radius**2, z) == 0)
    pressure = s.Symbol("p", real=True)
    ledger.check("pressure normal to the sphere exerts zero physical centroid torque",
                 position.cross(pressure*position/radius) == s.zeros(3, 1))
    omega0 = 2*lam*ez/3
    rotated_variation = s.Matrix([1, 0, 0]).cross(omega0)
    angle_variation = -rotated_variation[1]/omega0[2]
    ledger.check("core vorticity orientation is a physical rotation-equivariant angle",
                 angle_variation == 1)

    b, hessian, p = s.symbols("B h p", real=True, nonzero=True)
    qdot = s.Symbol("qdot", real=True)
    orbit_action = b*p*qdot-hessian*p*p/2
    eliminated = s.simplify(orbit_action.subs(p, b*qdot/hessian))
    ledger.check("same fixed-Kelvin orbit gives I_red=B²/h without adding geometric mass",
                 s.simplify(eliminated-b*b*qdot*qdot/(2*hessian)) == 0)
    ledger.check("physical reduced angular momentum is I_red*qdot",
                 s.simplify((b*p).subs(p, b*qdot/hessian)-s.diff(eliminated, qdot)) == 0)
    shape_hessian = s.Matrix([[4, 1], [1, 3]])
    moment_map = s.Matrix([[1], [2]])
    schur = hermitian_schur_jet((shape_hessian, s.zeros(2), s.zeros(2)),
                               (moment_map, s.zeros(2, 1), s.zeros(2, 1)),
                               (s.zeros(1), s.zeros(1), s.zeros(1)))
    full_reduced = -schur.reduced[0][0]
    ledger.check("several Kelvin reaction directions use a full positive inverse, not a locked metric",
                 full_reduced == (moment_map.T*shape_hessian.inv()*moment_map)[0]
                 and full_reduced > 0)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
