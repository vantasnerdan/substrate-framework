"""KKS/rotation normalization and exact curvature sign for Euler inertial modes."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0131-actual-rotation-action")
    omega, axial, kinetic, row_square = s.symbols("Omega N h row_square", positive=True)
    radial_square = s.Symbol("K_square", nonnegative=True)
    sigma = s.Symbol("sigma", nonzero=True, real=True)
    azimuth = s.Symbol("m", integer=True)
    amplitude1, amplitude2 = s.symbols("q1 q2", real=True)
    rotation = s.Matrix([[0, -1], [1, 0]])
    beta = -kinetic/sigma
    form = -beta*rotation
    intrinsic_generator = sigma*rotation
    laboratory_frequency = azimuth*omega+sigma
    laboratory_generator = laboratory_frequency*rotation
    physical_rotation_generator = azimuth*rotation
    intrinsic_hessian = -form*intrinsic_generator
    laboratory_hessian = -form*laboratory_generator
    charge_hessian = -form*physical_rotation_generator
    checks.check("rotating kinetic energy fixes the KKS sign from the actual generator",
                 s.simplify(intrinsic_hessian-kinetic*s.eye(2)) == s.zeros(2))
    checks.check("physical SO2 charge includes the full azimuthal generator factor",
                 s.simplify(charge_hessian+azimuth*beta*s.eye(2)) == s.zeros(2))
    checks.check("laboratory action retains the angular momentum contribution",
                 s.simplify(laboratory_hessian-intrinsic_hessian-omega*charge_hessian)
                 == s.zeros(2))
    checks.check("omitting azimuthal generator normalization changes the Hamiltonian",
                 s.simplify(charge_hessian+beta*s.eye(2)) != s.zeros(2))

    # Direct scalar chart: angle=c*q1, angle_rate=-c*omega_lab*q2.
    row = s.sqrt(row_square)
    observation = s.Matrix([[row, 0], [0, -row*laboratory_frequency]])
    wronskian = s.det(observation)
    physical_inertia = s.simplify(beta/wronskian)
    checks.check("scalar inertia comes from the physical chart Wronskian",
                 s.simplify(physical_inertia-kinetic/(row_square*sigma*laboratory_frequency)) == 0)
    for helicity in (-1, 1):
        intrinsic = helicity*2*omega*axial/s.sqrt(axial**2+radial_square)
        laboratory = azimuth*omega+intrinsic
        inertia = physical_inertia.subs(sigma, intrinsic)
        stiffness_jet = s.simplify(inertia*s.diff(laboratory**2, radial_square))
        checks.check(f"helicity {helicity}: positive inertia and positive curvature cannot coexist",
                     s.simplify(stiffness_jet+kinetic/(row_square*(axial**2+radial_square))) == 0)

    # A genuine SO2 pattern phase is angle=psi/m, not psi. The polar action
    # follows from the pulled KKS one-form; it remains only quadratic in field.
    t = s.Symbol("t", real=True)
    radius = s.Function("radius", positive=True)(t)
    angle = s.Function("angle", real=True)(t)
    vector = radius*s.Matrix([s.cos(azimuth*angle), s.sin(azimuth*angle)])
    canonical_term = s.simplify(-(vector.T*form*vector.diff(t))[0]/2)
    checks.check("physical pattern-angle momentum is minus m times wave action",
                 s.simplify(s.diff(canonical_term, s.diff(angle, t))
                            +azimuth*beta*radius**2/2) == 0)
    wave_action = s.Symbol("wave_action", positive=True)
    lagrangian = -azimuth*wave_action*s.diff(angle, t)+laboratory_frequency*wave_action
    checks.check("pattern phase is a cyclic drift, not an optical restoring coordinate",
                 s.diff(lagrangian, angle) == 0
                 and s.diff(lagrangian, s.diff(angle, t), 2) == 0)

    # Exposing full material-spin versus velocity-only mutation.
    r, pressure, derivative = s.symbols("r P Pprime", real=True, nonzero=True)
    denominator = 4*omega**2-sigma**2
    vr = s.I*(sigma*derivative-2*omega*azimuth*pressure/r)/denominator
    vp = (2*omega*derivative-sigma*azimuth*pressure/r)/denominator
    full_spin = r*(2*omega*s.I*vr/sigma+vp)
    checks.check("deleting the moving-tag contribution fails the pressure-spin identity",
                 s.simplify(r*vp-azimuth*pressure/sigma) != 0
                 and s.simplify(full_spin-azimuth*pressure/sigma) == 0)

    # The finite-amplitude pressure completion follows from the defining
    # Beltrami property, not a nonlinear frequency ansatz.
    x, y, z = s.symbols("x y z", real=True)
    phase = axial*z-omega*t
    real_velocity = s.Matrix([(x*x-y*y)*s.cos(phase)-2*x*y*s.sin(phase),
                              -(x*x-y*y)*s.sin(phase)-2*x*y*s.cos(phase), 0])
    coordinates = (x, y, z)
    kinetic_density = s.simplify(real_velocity.dot(real_velocity)/2)
    nonlinear = real_velocity.jacobian(coordinates)*real_velocity
    pressure_gradient = s.Matrix([s.diff(kinetic_density, coordinate)
                                  for coordinate in coordinates])
    checks.check("actual real polynomial mode has only a gradient nonlinear self-interaction",
                 s.simplify(nonlinear-pressure_gradient) == s.zeros(3, 1))
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
