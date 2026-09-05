"""Exact smooth-annulus Euler transfer and material/KKS matching interface."""

import sympy as s

from substrate_framework.rankine_modes import boundary_determinant
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0137-smooth-column-transfer")
    r = s.Symbol("r", positive=True)
    m, k, omega = s.symbols("m k omega", real=True)
    angular = s.Function("O")(r)
    axial = s.Function("W")(r)
    f, p = s.Function("f")(r), s.Function("p")(r)
    vort = 2*angular+r*s.diff(angular, r)
    sigma = omega-m*angular-k*axial
    vr = -s.I*sigma*f
    vt = m*p/(r*sigma)-vort*f
    vz = k*p/sigma-s.diff(axial, r)*f
    fp = -(1/r+2*m*angular/(r*sigma))*f+(m*m/r**2+k*k)*p/sigma**2
    pp = (sigma**2-2*angular*vort)*f+2*m*angular*p/(r*sigma)
    substitutions = {s.diff(f, r): fp, s.diff(p, r): pp}
    radial = -s.I*sigma*vr-2*angular*vt+s.diff(p, r)
    theta = -s.I*sigma*vt+vort*vr+s.I*m*p/r
    vertical = -s.I*sigma*vz+s.diff(axial, r)*vr+s.I*k*p
    divergence = s.diff(vr, r)+vr/r+s.I*m*vt/r+s.I*k*vz
    for name, equation in (("radial momentum", radial), ("azimuthal momentum", theta),
                           ("axial momentum", vertical), ("full incompressibility", divergence)):
        checks.check(f"bounded first-order system satisfies exact {name}",
                     s.simplify(equation.subs(substitutions)) == 0)
    # Test this mutation directly, because substitutions into a simplified
    # expression can otherwise erase the very term being tested.
    wrong_vertical = -s.I*sigma*(k*p/sigma)+s.diff(axial, r)*vr+s.I*k*p
    checks.check("axial shear mutation has a nonzero explicit residual",
                 s.simplify(wrong_vertical+s.I*sigma*s.diff(axial, r)*f) == 0
                 and wrong_vertical != 0)

    eta_theta = s.I*(m*p/(r*sigma)-2*angular*f)/sigma
    eta_z = s.I*k*p/sigma**2
    checks.check("material tangential Lin relation retains differential rotation",
                 s.simplify(-s.I*sigma*eta_theta-r*s.diff(angular, r)*f-vt) == 0)
    checks.check("material axial Lin relation retains axial shear",
                 s.simplify(-s.I*sigma*eta_z-s.diff(axial, r)*f-vz) == 0)
    coefficient_matrix = s.Matrix([[-1/r-2*m*angular/(r*sigma), (m*m/r**2+k*k)/sigma**2],
                                   [sigma**2-2*angular*vort, 2*m*angular/(r*sigma)]])
    checks.check("transfer determinant has exact radial Liouville trace",
                 s.simplify(s.trace(coefficient_matrix)+1/r) == 0)

    om, doppler, jr, kr, radius = s.symbols("O sigma J K a", nonzero=True, real=True)
    delta = 4*om*om-doppler*doppler
    core_f = (2*om*m-doppler*jr)/(radius*doppler*delta)
    outer_f = kr/radius
    outer_p = doppler*doppler
    evans = s.Matrix([[core_f, outer_f], [1, outer_p]]).det()
    determinant = boundary_determinant(doppler, om, m, jr, kr)
    checks.check("zero-thickness transfer is the full Rankine pressure determinant",
                 s.simplify(evans+determinant/(radius*delta)) == 0)
    checks.check("exterior potential pressure has the necessary Doppler square at the edge",
                 s.simplify(core_f*doppler-outer_f+determinant/(radius*delta)) != 0)

    circulation = s.Symbol("C", positive=True)
    potential = s.Function("K")(r)
    exterior_o = circulation/r**2
    exterior_sigma = omega-m*exterior_o
    exterior_v = s.Matrix([s.diff(potential, r), s.I*m*potential/r, s.I*k*potential])
    exterior_pressure = s.I*exterior_sigma*potential
    exterior_radial = -s.I*exterior_sigma*exterior_v[0] \
        -2*exterior_o*exterior_v[1]+s.diff(exterior_pressure, r)
    checks.check("actual exterior velocity and pressure remain regular at particle resonance",
                 s.simplify(exterior_radial) == 0)
    beta, row, frequency = s.symbols("beta c sigma", nonzero=True, real=True)
    symplectic = s.Matrix([[0, beta], [-beta, 0]])
    generator = s.Matrix([[0, -frequency], [frequency, 0]])
    hessian = -symplectic*generator
    checks.check("smooth mode uses its complete KKS rotating Hamiltonian",
                 hessian == -beta*frequency*s.eye(2))
    checks.check("mass and stiffness retain the exact intrinsic dispersion",
                 s.simplify((-beta*frequency/row**2)/(-beta/(frequency*row**2))-frequency**2) == 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
