"""Actual rotating-cell transport, physical current and whole-field moments."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0144-physical-mean-response")
    t, omega, moment = s.symbols("t Omega moment", positive=True)
    identity = s.eye(2)
    rotation = s.Matrix([[0, -1], [1, 0]])
    u1, u2, v1, v2 = s.symbols("U1 U2 V1 V2", real=True)
    displacement = s.Matrix([u1, u2])
    velocity = s.Matrix([v1, v2])
    q = s.Matrix([s.Function("q1")(t), s.Function("q2")(t)])
    # For actual v=Omega J r, p=Omega²r²/2 and <r r^T>=moment I,
    # chi_z=q.r. These are direct Euler field moments.
    material_acceleration = omega**2*moment*q
    current_derivative = omega*moment*rotation*q.diff(t)
    transported = q.diff(t)-omega*rotation*q
    ledger.check("physical mean current converts the Jacobi force to the transported axial row",
                 s.simplify(material_acceleration+current_derivative
                            -omega*moment*rotation*transported) == s.zeros(2, 1))
    kelvin_row = omega*rotation*displacement
    physical_force = omega*moment*rotation*kelvin_row
    ledger.check("actual Kelvin displacement phase has positive covariance bending stiffness",
                 physical_force == -omega**2*moment*displacement)
    ledger.check("the local cell tilt solves the actual scalar Jacobi equation",
                 (-omega*rotation)**2*(-displacement) == omega**2*displacement)

    # The actual planar flow is exp(Omega t J). Its TR-paired velocity
    # correlation is derived by the two signed rotations.
    propagator = s.cos(omega*t)*identity+s.sin(omega*t)*rotation
    covariance = omega**2*moment
    correlation = covariance*propagator
    paired = s.simplify((correlation+correlation.subs(omega, -omega))/2)
    ledger.check("whole-field time reversal retains the symmetric velocity correlation",
                 paired == covariance*s.cos(omega*t)*identity)
    tau = s.Symbol("tau", real=True)
    memory_m = s.integrate((t-tau)*s.cos(omega*tau), (tau, 0, t))
    memory_x = s.integrate((t-tau)**2*s.cos(omega*tau)/2, (tau, 0, t))
    ledger.check("common-velocity phase has the stated exact mean-current memory",
                 s.simplify(memory_m-(1-s.cos(omega*t))/omega**2) == 0)
    ledger.check("integrated physical mean has the exact displacement memory",
                 s.simplify(memory_x-(t/omega**2-s.sin(omega*t)/omega**3)) == 0)
    coefficient = covariance*(t*t/2-memory_m)
    ledger.check("an immediate acoustic PDE would miss the actual common-velocity initial jet",
                 s.simplify(s.diff(coefficient, t, 2).subs(t, 0)) == 0
                 and s.diff(covariance*t*t/2, t, 2) == covariance)
    exact_rate = -covariance*(t-s.sin(omega*t)/omega)*velocity
    row = t*omega*rotation*velocity-(propagator-identity)*velocity
    raw_rate = omega*moment*rotation*row
    paired_rate = s.simplify((raw_rate+raw_rate.subs(omega, -omega))/2)
    ledger.check("the full Kelvin-rate calculation independently reproduces the correlation formula",
                 s.simplify(paired_rate-exact_rate) == s.zeros(2, 1))

    gamma, cell_area, core, outer, r, theta = s.symbols("Gamma A a b r theta", positive=True)
    swirl = gamma/(2*s.pi*r)
    cov_x = s.integrate(s.integrate(swirl**2*s.sin(theta)**2*r,
                                    (theta, 0, 2*s.pi)), (r, core, outer))/cell_area
    ledger.check("the positive logarithmic bending coefficient is an actual velocity integral",
                 s.simplify(cov_x-gamma**2*s.log(outer/core)/(4*s.pi*cell_area)) == 0)

    # Fix K along z. The SH unit vector is azimuthal; its projection's
    # azimuthal average is P_K/2, independently of polar angle.
    z, phi = s.symbols("z phi", real=True)
    sh_xx = s.sin(phi)**2
    average_sh = s.integrate(s.integrate(sh_xx, (phi, 0, 2*s.pi)), (z, -1, 1))/(4*s.pi)
    average_axial_sh = s.integrate(s.integrate(z*z*sh_xx, (phi, 0, 2*s.pi)),
                                   (z, -1, 1))/(4*s.pi)
    ledger.check("the actual SH observation has one-half rank weight under the whole-field law",
                 average_sh == s.Rational(1, 2) and average_axial_sh == s.Rational(1, 6))
    rho, cp2, cb2 = s.symbols("rho cp2 cb2", positive=True)
    normalized_mass = rho*4*average_sh
    speed = (cp2*(average_sh-average_axial_sh)+cb2*average_axial_sh)/average_sh
    ledger.check("normalizing raw SH mean produces inertia two rho, not rho",
                 normalized_mass == 2*rho and normalized_mass != rho)
    ledger.check("SH covariance benchmark has its derived orientation weights",
                 s.simplify(speed-(2*cp2+cb2)/3) == 0)
    print("Common-V paired mean coefficient:", s.simplify(coefficient))
    print("Covariance bending integral:", s.simplify(cov_x))
    print("SH-only normalized mass:", normalized_mass)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
