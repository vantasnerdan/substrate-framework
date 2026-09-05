"""Exact regressions of material transport, full spin, and Jacobi signs."""

import sympy as s

from substrate_framework.euler_displacement import (
    euler_displacement_perturbation,
    euler_jacobi_density,
    material_derivative,
)
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0084-material-transport")
    x, y, z, t = s.symbols("x y z t", real=True)
    coords = (x, y, z)
    r = s.Matrix(coords)
    rho = s.symbols("rho", positive=True)

    def curl(v):
        return s.Matrix([s.diff(v[2], y)-s.diff(v[1], z),
                         s.diff(v[0], z)-s.diff(v[2], x),
                         s.diff(v[1], x)-s.diff(v[0], y)])

    u = s.Matrix([-y, x, 0])
    xi = s.Matrix([x*z+t*x, -y*z+t*y, -2*t*z])
    ledger.check("test material field is divergence free",
                 s.trace(xi.jacobian(coords)) == 0)
    v = euler_displacement_perturbation(xi, u, coords, t)
    dt_xi = material_derivative(xi, u, coords, t)
    material_spin = xi.cross(u)+r.cross(dt_xi)
    local_spin = r.cross(xi.diff(t))+2*xi.cross(u)
    transport_flux = (r.cross(xi)).jacobian(coords)*u
    ledger.check("full material spin differs from local formula only by background boundary flux",
                 s.simplify(material_spin-local_spin-transport_flux) == s.zeros(3, 1))
    boundary_integrand = (r.cross(u)).jacobian(coords)*xi
    ledger.check("Eulerian velocity spin plus moving-boundary variation equals material spin",
                 s.simplify(r.cross(v)+boundary_integrand-material_spin) == s.zeros(3, 1))
    ledger.check("omitting the factor two in the local spin is exposed",
                 s.simplify(material_spin-r.cross(xi.diff(t))-xi.cross(u)-transport_flux)
                 != s.zeros(3, 1))

    b = s.Function("b")(x, y, t)
    surface_flow = s.Matrix([x, 0, -z])
    surface_xi = s.Matrix([0, 0, b])
    normal_velocity = euler_displacement_perturbation(
        surface_xi, surface_flow, coords, t)[2]
    ledger.check("surface displacement obeys density transport including surface divergence",
                 s.simplify(normal_velocity-s.diff(b, t)-s.diff(x*b, x)) == 0)
    ledger.check("bare scalar surface advection omits a nonzero geometric contribution",
                 s.simplify(normal_velocity-s.diff(b, t)-x*s.diff(b, x)) == b)
    psi = s.Function("psi")(x, y)
    cutoff = s.Function("cutoff")(z)
    potential = cutoff*s.Matrix([-s.diff(psi, y), s.diff(psi, x), 0])
    collar = curl(potential)
    ledger.check("collar curl is divergence free and has the prescribed normal Laplacian",
                 s.simplify(s.trace(collar.jacobian(coords))) == 0
                 and s.simplify(collar[2]-cutoff*(s.diff(psi, x, 2)+s.diff(psi, y, 2))) == 0)

    kelvin = dt_xi+xi.jacobian(coords).T*u
    grad_pairing = s.Matrix([s.diff(u.dot(xi), c) for c in coords])
    ledger.check("material Kelvin one-form equals velocity minus isovortical term modulo exact form",
                 s.simplify(kelvin-v+xi.cross(curl(u))-grad_pairing) == s.zeros(3, 1))
    ledger.check("arbitrary material restriction does not silently fix its Kelvin momentum",
                 s.simplify(curl(kelvin)) != s.zeros(3, 1))

    pressure = rho*(x*x+y*y)/2
    ledger.check("nonconstant-pressure fixture obeys stationary Euler with the physical sign",
                 s.simplify(rho*u.jacobian(coords)*u
                            +s.Matrix([s.diff(pressure, c) for c in coords])) == s.zeros(3, 1))
    static_xi = xi.subs(t, 0)
    static_density = euler_jacobi_density(static_xi, u, pressure, rho, coords, t)
    k_density = -2*static_density
    f = static_xi.cross(u)
    w = curl(f)
    curvature = rho*(w.dot(static_xi.cross(curl(u)))-w.dot(w))
    flux_scalar = static_xi.dot(u.jacobian(coords)*static_xi)
    flux_div = rho*(s.Matrix([s.diff(flux_scalar, c) for c in coords]).dot(u))
    ledger.check("stationary Jacobi stiffness has the full curl-curvature identity and boundary flux",
                 s.simplify(k_density-curvature+flux_div) == 0)
    ledger.check("changing the curvature helicity sign is exposed",
                 s.simplify(k_density+rho*w.dot(static_xi.cross(curl(u)))
                            +rho*w.dot(w)+flux_div) != 0)

    lam = s.symbols("lambda", positive=True)
    wave = s.Matrix([s.sin(lam*z), s.cos(lam*z), 0])
    ledger.check("sign-exposing background is an exact smooth Beltrami field",
                 s.simplify(curl(wave)-lam*wave) == s.zeros(3, 1)
                 and s.simplify(wave.jacobian(coords)*wave) == s.zeros(3, 1))
    wave_density = euler_jacobi_density(static_xi, wave, s.S.Zero, rho, coords, t)
    advected = static_xi.jacobian(coords)*wave
    ledger.check("constant-pressure material stiffness is minus the full convective norm",
                 s.simplify(-2*wave_density+rho*advected.dot(advected)) == 0
                 and s.simplify(advected.dot(advected)) != 0)
    ledger.check("Beltrami material commutator is curl of displacement cross background",
                 s.simplify(curl(static_xi.cross(wave))
                            -static_xi.jacobian(coords)*wave
                            +wave.jacobian(coords)*static_xi) == s.zeros(3, 1))
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
