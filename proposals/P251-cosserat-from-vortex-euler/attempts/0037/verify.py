"""Exact full-variation and 3-D scope checks, independent of pytest fixtures."""

import sympy as sp

from substrate_framework.euler_displacement import (
    euler_displacement_perturbation,
    euler_jacobi_density,
    material_derivative,
)
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0037")
    x, y, z, t, eps = sp.symbols("x y z t epsilon", real=True)
    rho = sp.Symbol("rho", positive=True)
    coords = (x, y, z)
    a = sp.Matrix(3, 3, sp.symbols("a:9"))
    determinant = (sp.eye(3)+eps*a).det()
    ledger.check("complete determinant second variation",
                 sp.diff(determinant, eps, 2).subs(eps, 0)
                 == sp.expand(sp.trace(a)**2-sp.trace(a*a)))

    # General planar stream displacement embedded in 3-D licenses div xi=0.
    f = sp.Function("f")(x, y)
    xi = sp.Matrix([-sp.diff(f, y), sp.diff(f, x), 0])
    p = sp.Function("p")(x, y, z)
    gradp = sp.Matrix([sp.diff(p, c) for c in coords])
    flux = p*xi.jacobian(coords)*xi-xi*(xi.dot(gradp))
    divergence = sum(sp.diff(flux[i], coords[i]) for i in range(3))
    lhs = p*sp.trace(xi.jacobian(coords)**2)
    ledger.check("pressure constraint integration by parts with explicit flux",
                 sp.simplify(lhs-(xi.T*sp.hessian(p, coords)*xi)[0]-divergence) == 0)

    # Smooth nonconstant-pressure 3-D Beltrami field, not an EPS topology claim.
    u = sp.Matrix([sp.sin(z)+sp.cos(y), sp.sin(x)+sp.cos(z), sp.sin(y)+sp.cos(x)])
    curl = sp.Matrix([sp.diff(u[2], y)-sp.diff(u[1], z),
                      sp.diff(u[0], z)-sp.diff(u[2], x),
                      sp.diff(u[1], x)-sp.diff(u[0], y)])
    p0 = -rho*u.dot(u)/2
    ledger.check("ABC field is Beltrami", sp.simplify(curl-u) == sp.zeros(3, 1))
    ledger.check("ABC field is on-shell stationary Euler",
                 sp.simplify(rho*u.jacobian(coords)*u +
                             sp.Matrix([sp.diff(p0, c) for c in coords])) == sp.zeros(3, 1))
    # Translation plus constant Galilean velocity is an exact linearized solution.
    shift = sp.Matrix([1+2*t, -1+3*t, 2-t])
    delta_p = -shift.dot(sp.Matrix([sp.diff(p0, c) for c in coords]))
    jacobi = rho*material_derivative(material_derivative(shift, u, coords, t), u, coords, t)
    jacobi += sp.hessian(p0, coords)*shift
    jacobi += sp.Matrix([sp.diff(delta_p, c) for c in coords])
    ledger.check("translation and Galilean motion retain pressure reaction",
                 sp.simplify(jacobi) == sp.zeros(3, 1))
    # A relabeling direction has zero Eulerian perturbation despite finite norm.
    ledger.check("steady streamline relabeling has zero physical velocity tangent",
                 euler_displacement_perturbation(u, u, coords, t) == sp.zeros(3, 1))
    gauge_jacobi = rho*material_derivative(material_derivative(u, u, coords, t), u, coords, t)
    gauge_jacobi += sp.hessian(p0, coords)*u
    ledger.check("relabeling is retained by the complete Jacobi operator",
                 sp.simplify(gauge_jacobi) == sp.zeros(3, 1))

    # A 3-D Coriolis field is not generically a pressure gradient.
    jw = sp.Matrix([-u[1], u[0], 0])
    curl_jw = sp.Matrix([sp.diff(jw[2], y)-sp.diff(jw[1], z),
                         sp.diff(jw[0], z)-sp.diff(jw[2], x),
                         sp.diff(jw[1], x)-sp.diff(jw[0], y)])
    ledger.check("3-D stationarization obstruction is axial variation",
                 sp.simplify(curl_jw+u.diff(z)) == sp.zeros(3, 1))
    ledger.check("ABC Coriolis term cannot be absorbed into pressure",
                 sp.simplify(curl_jw) != sp.zeros(3, 1))

    # Time-kinetic and gyroscopic blocks derived from the complete density.
    q = sp.Function("q")(t)
    r = sp.Function("r")(t)
    profile = sp.Matrix([q*sp.sin(y)+r*sp.cos(y), 0, 0])
    density = euler_jacobi_density(profile, u, p0, rho, coords, t)
    mass = sp.diff(density, sp.diff(q, t), 2)
    ledger.check("material time inertia follows from Euler action",
                 sp.simplify(mass-rho*sp.sin(y)**2) == 0)
    gyro = sp.diff(density, sp.diff(q, t), r)-sp.diff(density, sp.diff(r, t), q)
    ledger.check("omitting transport loses a nonzero gyroscopic block",
                 sp.simplify(gyro+rho*u[1]) == 0 and gyro != 0)
    print("Exact smooth Euler identities; no mode positivity or Cosserat closure inferred.")
    return int(ledger.finish())


if __name__ == "__main__":
    raise SystemExit(main())
