"""Exact ellipse-sector angle oscillator; conditional continuum contraction.

Units per unit tube length: H and K are N, p is kg m/s, I is kg m.
This construction does not establish a closed interacting Euler ensemble.
"""

import sympy as sp

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0030-local-angle")
    q, c, omega, strain = sp.symbols("Q C omega e", positive=True)
    theta, rotation = sp.symbols("theta Omega", real=True)
    d = sp.sqrt(q**2 - 4)
    hamiltonian = (
        -c * omega * sp.log(q + 2)
        + c * strain * d * sp.sin(2 * theta)
        - c * rotation * q
    )
    # Canonical angle momentum p=-C Q is minus angular impulse per length.
    angle_rate = -sp.diff(hamiltonian, q) / c
    q_rate = sp.diff(hamiltonian, theta) / c
    expected_angle_rate = rotation + omega / (q + 2) - strain * q / d * sp.sin(2 * theta)
    ledger.check("Kida angle equation", sp.simplify(angle_rate - expected_angle_rate) == 0)
    ledger.check("Kida shape equation", sp.simplify(q_rate - 2 * strain * d * sp.cos(2 * theta)) == 0)
    lam = sp.symbols("lambda", positive=True)
    q_lam, d_lam = lam + 1 / lam, lam - 1 / lam
    ledger.check(
        "area-preserving aspect-ratio transport",
        sp.simplify(2 * strain * d_lam / sp.diff(q_lam, lam) - 2 * strain * lam) == 0,
    )
    ledger.check(
        "Kirchhoff self-rotation normalization",
        sp.simplify(omega / (q_lam + 2) - omega * lam / (lam + 1)**2) == 0,
    )

    equilibrium_rotation = -omega / (q + 2) - strain * q / d
    equilibrium = {theta: -sp.pi / 4, rotation: equilibrium_rotation}
    ledger.check("exact stationary angle", sp.simplify(angle_rate.subs(equilibrium)) == 0)
    ledger.check("exact stationary shape", sp.simplify(q_rate.subs(equilibrium)) == 0)
    stiffness = sp.simplify(sp.diff(hamiltonian, theta, 2).subs(equilibrium))
    inverse_inertia = sp.diff(hamiltonian, q, 2).subs(equilibrium) / c**2
    expected_inverse = (omega / (q + 2)**2 + 4 * strain / d**3) / c
    ledger.check("angle stiffness from full Hessian", sp.simplify(stiffness - 4 * c * strain * d) == 0)
    ledger.check("canonical momentum Hessian", sp.simplify(inverse_inertia - expected_inverse) == 0)
    ledger.check("mixed Hessian vanishes", sp.simplify(sp.diff(hamiltonian, theta, q).subs(equilibrium)) == 0)
    positive_d = sp.symbols("D_positive", positive=True)
    positive_chart = {q: sp.sqrt(positive_d**2 + 4)}
    ledger.check("positive angle stiffness for Q>2", sp.simplify(stiffness.subs(positive_chart)).is_positive)
    ledger.check("positive canonical Hessian for Q>2", sp.simplify(expected_inverse.subs(positive_chart)).is_positive)

    momentum, angle, angle_dot = sp.symbols("delta_p delta_theta delta_theta_dot", real=True)
    inv_i, spring = sp.symbols("inverse_I K", positive=True)
    first_order = momentum * angle_dot - (inv_i * momentum**2 + spring * angle**2) / 2
    momentum_solution = sp.solve(sp.diff(first_order, momentum), momentum)[0]
    reduced = sp.simplify(first_order.subs(momentum, momentum_solution))
    ledger.check("Legendre elimination supplies inertia", sp.simplify(reduced - (angle_dot**2 / inv_i - spring * angle**2) / 2) == 0)
    matrix = sp.Matrix([[0, inv_i], [-spring, 0]])
    ledger.check("canonical oscillator frequency", sp.simplify(matrix.det() - inv_i * spring) == 0)
    ledger.check("strain-free angle stiffness vanishes", sp.simplify(stiffness.subs(strain, 0)) == 0)
    ledger.mutation_sensitive(
        "canonical momentum sign",
        lambda sign: sp.simplify(sign * sp.diff(hamiltonian, q) / c - expected_angle_rate) == 0,
        -1,
        [1],
    )

    # Independent geometric normalization from exact area moments.
    rho, gamma, area, a, b = sp.symbols("rho Gamma area a b", positive=True)
    radius, azimuth = sp.symbols("radius azimuth", real=True)
    x = radius * (a * sp.cos(azimuth) * sp.cos(theta) - b * sp.sin(azimuth) * sp.sin(theta))
    y = radius * (a * sp.cos(azimuth) * sp.sin(theta) + b * sp.sin(azimuth) * sp.cos(theta))
    xy_average = sp.integrate(sp.integrate(sp.expand_trig(x * y) * radius / sp.pi, (azimuth, 0, 2 * sp.pi)), (radius, 0, 1))
    ledger.check("exact ellipse cross moment", sp.simplify(xy_average - (a**2 - b**2) * sp.sin(2 * theta) / 8) == 0)
    c_physical = rho * gamma * area / (8 * sp.pi)
    cross_energy = rho * gamma * strain * xy_average
    ledger.check(
        "physical background interaction coefficient",
        sp.simplify(cross_energy - c_physical * strain * sp.pi * (a**2 - b**2) / area * sp.sin(2 * theta)) == 0,
    )

    # Exact frozen-shape mutual potential: its convergent exterior expansion
    # supplies orientation dependence without assuming a micropolar potential.
    separation, gamma_two, beta = sp.symbols("separation Gamma_two beta", real=True)
    coefficient = rho * gamma * gamma_two * (a**2 - b**2) / (16 * sp.pi * separation**2)
    mutual_quadrupole = coefficient * sp.cos(2 * (theta - beta))
    ledger.check(
        "mutual interaction action and reaction",
        sp.simplify(sp.diff(mutual_quadrupole, theta) + sp.diff(mutual_quadrupole, beta)) == 0,
    )
    ledger.check(
        "mutual angle Hessian at transverse orientation",
        sp.simplify(sp.diff(mutual_quadrupole, theta, 2).subs(theta, beta + sp.pi / 2) - 4 * coefficient) == 0,
    )
    z = radius * (a * sp.cos(azimuth) + sp.I * b * sp.sin(azimuth))
    for order in (1, 2, 3):
        moment = sp.integrate(sp.integrate(sp.expand(z**(2 * order)) * radius / sp.pi, (azimuth, 0, 2 * sp.pi)), (radius, 0, 1))
        ledger.check(
            f"ellipse exterior moment {2 * order}",
            sp.simplify(moment - sp.catalan(order) * ((a**2 - b**2) / 4)**order) == 0,
        )

    # The following is an explicit conditional map, not a closed-Euler claim.
    line_density, inertia = sp.symbols("L_v I", positive=True)
    relative_angle = sp.symbols("relative_angle", real=True)
    averaged_energy = line_density * spring * relative_angle**2 / 6
    alpha = sp.solve(sp.Eq(averaged_energy, 2 * sp.Symbol("alpha") * relative_angle**2), sp.Symbol("alpha"))[0]
    microinertia = line_density * inertia / 3
    ledger.check("conditional sphere-contracted optical gap", sp.simplify(4 * alpha / microinertia - spring / inertia) == 0)

    print("H =", hamiltonian)
    print("Omega_stationary =", equilibrium_rotation)
    print("K_angle =", stiffness)
    print("I_angle =", 1 / expected_inverse)
    print("alpha_if_material_frame_bridge =", alpha)
    print("Scope: exact ellipse-sector oscillator and frozen mutual-energy moments.")
    print("Open: self-consistent ambient Euler action, spatial closure, EPS bridge.")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
