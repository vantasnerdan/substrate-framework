"""Closed point-vortex angle action and explicit continuum correspondence limits."""

import sympy as sp

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0032-closed-hexagon")
    rho, gamma, scale = sp.symbols("rho Gamma S", positive=True)
    imbalance, angle = sp.symbols("x chi", real=True)
    amplitude = 3 * rho * gamma**2 / (4 * sp.pi)
    momentum_scale = 3 * rho * gamma * scale / 4
    denominator = (
        2 + 6 * imbalance**2
        - 2 * (1 - imbalance**2)**sp.Rational(3, 2) * sp.cos(3 * angle)
    )
    reduced_h = amplitude * (-sp.log(1 - imbalance**2) - sp.log(denominator))
    equilibrium = {imbalance: 0, angle: sp.pi / 3}
    ledger.check("stationary radius imbalance", sp.diff(reduced_h, imbalance).subs(equilibrium) == 0)
    ledger.check("stationary relative angle", sp.diff(reduced_h, angle).subs(equilibrium) == 0)
    h_xx = sp.simplify(sp.diff(reduced_h, imbalance, 2).subs(equilibrium))
    stiffness = sp.simplify(sp.diff(reduced_h, angle, 2).subs(equilibrium))
    mixed = sp.simplify(sp.diff(reduced_h, angle, imbalance).subs(equilibrium))
    ledger.check("positive full radial Hessian", sp.simplify(h_xx - amplitude / 2) == 0)
    ledger.check("positive full angle Hessian", sp.simplify(stiffness - 9 * amplitude / 2) == 0)
    ledger.check("mixed Hessian vanishes", mixed == 0)
    inertia = sp.simplify(momentum_scale**2 / h_xx)
    frequency_squared = sp.simplify(stiffness / inertia)
    ledger.check("derived positive relative inertia", sp.simplify(inertia - 3 * sp.pi * rho * scale**2 / 2) == 0)
    ledger.check("derived relative frequency", sp.simplify(frequency_squared - 9 * gamma**2 / (4 * sp.pi**2 * scale**2)) == 0)

    # Roots-of-unity product proves the mutual logarithmic energy exactly.
    z, w, radius_one, radius_two = sp.symbols("z w r1 r2", real=True)
    cube_root = (-1 + sp.sqrt(3) * sp.I) / 2
    product = sp.expand(sp.prod(z - w * cube_root**index for index in range(3)))
    ledger.check("triangle pair product identity", sp.simplify(product - (z**3 - w**3)) == 0)
    relative_factor = radius_one**6 + radius_two**6 - 2 * radius_one**3 * radius_two**3 * sp.cos(3 * angle)
    factor_on_fixed_impulse = relative_factor.subs({radius_one: sp.sqrt(scale * (1 + imbalance) / 2), radius_two: sp.sqrt(scale * (1 - imbalance) / 2)})
    # Squared factor removes real-domain radical rewriting ambiguity |x|<1.
    ledger.check(
        "fixed-impulse radial polynomial terms",
        sp.simplify((factor_on_fixed_impulse + 2 * (scale * (1 + imbalance) / 2)**sp.Rational(3, 2) * (scale * (1 - imbalance) / 2)**sp.Rational(3, 2) * sp.cos(3 * angle)) - scale**3 * (2 + 6 * imbalance**2) / 8) == 0,
    )
    theta_one, theta_two = sp.symbols("theta1 theta2", real=True)
    mutual_h = -amplitude * sp.log(relative_factor.subs(angle, theta_one - theta_two))
    ledger.check("equal and opposite internal torques", sp.simplify(sp.diff(mutual_h, theta_one) + sp.diff(mutual_h, theta_two)) == 0)

    # Independent Cartesian Biot-Savart Jacobian at the regular hexagon.
    radius = sp.symbols("r", positive=True)
    rotation_generator = sp.Matrix([[0, -1], [1, 0]])
    points = [radius * sp.Matrix([sp.cos(index * sp.pi / 3), sp.sin(index * sp.pi / 3)]) for index in range(6)]
    velocities = []
    jacobian = sp.zeros(12)
    for index, point in enumerate(points):
        velocity = sp.zeros(2, 1)
        for other, other_point in enumerate(points):
            if other == index:
                continue
            difference = point - other_point
            squared_distance = sp.simplify(difference.dot(difference))
            velocity += gamma * rotation_generator * difference / (2 * sp.pi * squared_distance)
            block = gamma * rotation_generator * (sp.eye(2) / squared_distance - 2 * difference * difference.T / squared_distance**2) / (2 * sp.pi)
            for row in range(2):
                for column in range(2):
                    jacobian[2 * index + row, 2 * index + column] += block[row, column]
                    jacobian[2 * index + row, 2 * other + column] -= block[row, column]
        velocities.append(sp.simplify(velocity))
    orbital_rate = 5 * gamma / (4 * sp.pi * radius**2)
    ledger.check("direct Biot-Savart relative equilibrium", all(sp.simplify(velocity - orbital_rate * rotation_generator * point) == sp.zeros(2, 1) for point, velocity in zip(points, velocities)))
    for index in range(6):
        for row in range(2):
            for column in range(2):
                jacobian[2 * index + row, 2 * index + column] -= orbital_rate * rotation_generator[row, column]
    radial_mode = sp.Matrix.vstack(*[(-1)**index * point / radius for index, point in enumerate(points)])
    tangent_mode = sp.Matrix.vstack(*[(-1)**index * rotation_generator * point / radius for index, point in enumerate(points)])
    radial_to_tangent = -gamma / (4 * sp.pi * radius**2)
    tangent_to_radial = 9 * gamma / (4 * sp.pi * radius**2)
    ledger.check("Cartesian radial-to-angle mode", sp.simplify(jacobian * radial_mode - radial_to_tangent * tangent_mode) == sp.zeros(12, 1))
    ledger.check("Cartesian angle-to-radial mode", sp.simplify(jacobian * tangent_mode - tangent_to_radial * radial_mode) == sp.zeros(12, 1))
    ledger.check("independent frequency agreement", sp.simplify(-radial_to_tangent * tangent_to_radial - frequency_squared.subs(scale, 2 * radius**2)) == 0)
    ledger.mutation_sensitive(
        "angle stiffness normalization",
        lambda multiplier: sp.simplify(multiplier * stiffness / inertia + radial_to_tangent * tangent_to_radial).subs(scale, 2 * radius**2).simplify() == 0,
        1,
        [2, -1],
    )

    # Restore the common impulse: the actual action retains cross inertia.
    total_momentum = sp.symbols("P", negative=True)
    common_inverse_inertia = sp.diff(-5 * amplitude * sp.log(-total_momentum), total_momentum, 2)
    common_inertia = sp.simplify((1 / common_inverse_inertia).subs(total_momentum, -2 * momentum_scale))
    ledger.check("common inertia from unreduced Hamiltonian", sp.simplify(common_inertia - 2 * inertia / 5) == 0)
    micro_rate, cage_rate = sp.symbols("Phi_dot beta_dot", real=True)
    full_kinetic = common_inertia * (micro_rate + cage_rate)**2 / 8 + inertia * (micro_rate - cage_rate)**2 / 2
    cross_inertia = sp.simplify(sp.diff(full_kinetic, micro_rate, cage_rate))
    ledger.check("load-bearing cross inertia", sp.simplify(cross_inertia + 9 * inertia / 10) == 0)
    ledger.check("fixed-impulse kinetic is relative", sp.expand(inertia * (micro_rate - cage_rate)**2 / 2 - inertia * micro_rate**2 / 2) != 0)

    # Affine rotation of a triangular cage is derived from its three directions.
    h00, h01, h10, h11 = sp.symbols("h00 h01 h10 h11", real=True)
    gradient = sp.Matrix([[h00, h01], [h10, h11]])
    cage_angle = sum(((rotation_generator * points[2 * index] / radius).dot(gradient * points[2 * index] / radius)) for index in range(3)) / 3
    ledger.check("affine triangular orientation map", sp.simplify(cage_angle - (h10 - h01) / 2) == 0)

    # Local-induction helical stretch has explicitly named phase variables.
    tension = sp.symbols("T", positive=True)
    mean_gradient, relative_gradient = sp.symbols("B_prime chi_prime", real=True)
    helical_energy = 3 * tension * scale / 4 * ((mean_gradient + relative_gradient / 2)**2 + (mean_gradient - relative_gradient / 2)**2)
    relative_twist = sp.diff(helical_energy, relative_gradient, 2)
    common_twist = sp.diff(helical_energy, mean_gradient, 2)
    ledger.check("relative twist line-tension coefficient", sp.simplify(relative_twist - 3 * tension * scale / 4) == 0)
    ledger.check("common twist line-tension coefficient", sp.simplify(common_twist - 3 * tension * scale) == 0)
    print("K_relative =", stiffness)
    print("I_relative =", inertia)
    print("frequency_squared =", frequency_squared)
    print("I_common =", common_inertia)
    print("cross_inertia_Phi_beta =", cross_inertia)
    print("C_relative_log_leading =", relative_twist)
    print("Scope: exact closed point-vortex optical angle action; conditional filament gradient.")
    print("Open: finite-core spectral continuation, standard kinetic map, EPS and continuum closure.")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
