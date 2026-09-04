"""Direct stationary tube actions and a physical axial-polarization continuation.

Exact conditional actions on named admissible displacement ensembles. No
claim of dynamical stability follows from a restricted static coefficient.
"""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0038-stationary-tube")
    rho, r, lam = s.symbols("rho r lambda", positive=True)
    v, w = s.Function("V")(r), s.Function("W")(r)
    relations = {s.diff(v, r): lam*w-v/r, s.diff(w, r): -lam*v}
    cylindrical_curl = s.Matrix([0, -s.diff(w, r), s.diff(v, r)+v/r])
    ledger.check("Lundquist radial identities imply the full Beltrami equation",
                 s.simplify(cylindrical_curl.subs(relations)-lam*s.Matrix([0, v, w]))
                 == s.zeros(3, 1))
    pressure = -rho*(v**2+w**2)/2
    radial_pressure = s.simplify(s.diff(pressure, r).subs(relations))
    ledger.check("Beltrami pressure has the exact centripetal radial derivative",
                 s.simplify(radial_pressure-rho*v**2/r) == 0)
    q, qt, qz = s.symbols("q q_t q_z", real=True)
    material_velocity = s.Matrix([-v*q, r*(qt+w*qz), 0])
    pressure_form = r*radial_pressure*q**2
    jacobi_density = s.expand((rho*material_velocity.dot(material_velocity)-pressure_form)/2)
    ledger.check("pressure cancels apparent rigid-spin angle stiffness exactly",
                 s.simplify(jacobi_density-rho*r**2*(qt+w*qz)**2/2) == 0)
    ledger.mutation_sensitive(
        "pressure-Hessian coefficient",
        lambda factor: s.simplify((rho*material_velocity.dot(material_velocity)
                                  -factor*pressure_form)/2
                                 -rho*r**2*(qt+w*qz)**2/2) == 0,
        1, [0, 2],
    )

    # All following spatial coordinates are dimensionless; physical period is 2*pi*ell.
    x, y = s.symbols("x y", real=True)
    ell, a, b = s.symbols("ell a b", positive=True)
    field = s.Matrix([-b*s.sin(y), a*s.sin(x), a*s.cos(x)+b*s.cos(y)])
    def grad(f):
        return s.Matrix([s.diff(f, x), s.diff(f, y), 0])/ell

    def advect(vec):
        return vec.applyfunc(lambda f: (field[0]*s.diff(f, x)
                                        +field[1]*s.diff(f, y))/ell)
    curl = s.Matrix([s.diff(field[2], y), -s.diff(field[2], x),
                     s.diff(field[1], x)-s.diff(field[0], y)])/ell
    pressure = -rho*a*b*s.cos(x)*s.cos(y)
    ledger.check("nonaxisymmetric field is exactly Beltrami on the flat torus",
                 s.simplify(curl-field/ell) == s.zeros(3, 1))
    ledger.check("same pressure solves stationary Euler",
                 s.simplify(rho*advect(field)+grad(pressure)) == s.zeros(3, 1))
    psi = a*s.cos(x)+b*s.cos(y)
    ledger.check("psi level surfaces are actual invariant tube surfaces",
                 s.simplify(field.dot(grad(psi))) == 0)
    core_hessian = s.hessian(psi, (x, y)).subs({x: 0, y: 0})/ell**2
    ledger.check("positive a,b give an elliptic core and closed nearby level curves",
                 core_hessian == s.diag(-a/ell**2, -b/ell**2))

    angle_generator = ell*s.Matrix([-s.sin(y), s.sin(x), 0])
    ledger.check("angle generator is volume preserving",
                 (s.diff(angle_generator[0], x)+s.diff(angle_generator[1], y))/ell == 0)
    angular_jacobian = angle_generator[:2, :].jacobian((x, y))/ell
    ledger.check("angle generator rotates the core cross-section to first order",
                 angular_jacobian.subs({x: 0, y: 0}) == s.Matrix([[0, -1], [1, 0]]))
    generator_advection = advect(angle_generator)
    field_gradient = field.jacobian((x, y)) / ell
    eulerian_angle = s.simplify(generator_advection
                               -field_gradient*angle_generator[:2, :])
    expected_angle = (b-a)*s.Matrix([s.sin(x)*s.cos(y),
                                    -s.sin(y)*s.cos(x), s.sin(x)*s.sin(y)])
    ledger.check("noncircular core angle is a physical Eulerian perturbation",
                 s.simplify(eulerian_angle-expected_angle) == s.zeros(3, 1))
    ledger.check("circular equal-amplitude case is the relabeling limit",
                 eulerian_angle.subs(b, a) == s.zeros(3, 1))

    def mean(expr):
        return s.simplify(s.integrate(s.integrate(s.expand_trig(s.expand(expr)),
                                                (x, -s.pi, s.pi)),
                                     (y, -s.pi, s.pi))/(4*s.pi**2))

    pressure_hessian = s.zeros(3)
    pressure_hessian[:2, :2] = s.hessian(pressure, (x, y))/ell**2
    mass = rho*mean(angle_generator.dot(angle_generator))
    stiffness = mean((angle_generator.T*pressure_hessian*angle_generator)[0]
                     -rho*generator_advection.dot(generator_advection))
    ledger.check("one-action angle mass per volume is rho*ell^2", mass == rho*ell**2)
    ledger.check("one-action angle stiffness is minus rho*(a-b)^2/4",
                 s.simplify(stiffness+rho*(a-b)**2/4) == 0)
    ledger.check("single-angle gyroscopic term vanishes by periodic advection",
                 mean(angle_generator.dot(generator_advection)) == 0)
    for index in (0, 1):
        translation = s.eye(3)[:, index]
        ledger.check(f"translation {index} has mass rho from the same action",
                     rho*mean(translation.dot(translation)) == rho)
        ledger.check(f"translation {index} and angle have no uniform mass cross term",
                     mean(translation.dot(angle_generator)) == 0)
        ledger.check(f"translation {index} and angle have no uniform gyroscopic term",
                     mean(translation.dot(generator_advection)) == 0)
        ledger.check(f"translation {index} and angle have no uniform spring cross term",
                     mean((translation.T*pressure_hessian*angle_generator)[0]) == 0)
    ledger.check("uniform translation spring tensor vanishes",
                 pressure_hessian.applyfunc(mean) == s.zeros(3))

    # A materially different physical core polarization: arbitrary axial velocity.
    theta, t = s.symbols("theta t", real=True)
    azimuthal_m = s.symbols("m", integer=True, positive=True)
    amplitude = s.Function("A")(r)
    omega = v/r
    axial_polarization = amplitude*s.cos(azimuthal_m*(theta-omega*t))
    ledger.check("physical smooth axial polarization obeys exact Euler advection",
                 s.simplify(s.diff(axial_polarization, t)
                            +omega*s.diff(axial_polarization, theta)) == 0)
    phase = s.symbols("phi", real=True)
    # m=1 suffices for the phase-stiffness construction and avoids conditional integrals.
    physical_polarization = amplitude*s.cos(theta-phase)
    energy = s.integrate(rho*((w+physical_polarization)**2-w**2)/2,
                         (theta, -s.pi, s.pi))
    ledger.check("nonzero physical polarization energy is independent of its angle",
                 s.simplify(energy-rho*s.pi*amplitude**2/2) == 0)
    ledger.check("physical phase therefore has exactly zero restoring stiffness",
                 s.simplify(s.diff(energy, phase, 2)) == 0)
    print("Projected Beltrami angle mass/volume:", mass)
    print("Projected Beltrami angle stiffness/volume:", stiffness)
    print("Physical Eulerian angle generator:", eulerian_angle.T)
    print("Axial-polarization energy per r dr dz:", energy)
    print("Scope: exact stationary fields and conditional action restrictions; no stability verdict.")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
