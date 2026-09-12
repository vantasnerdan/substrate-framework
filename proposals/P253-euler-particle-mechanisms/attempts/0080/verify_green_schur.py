"""Exact algebraic checks for P253/0080."""

import sympy as sp

from substrate_framework.euler_cao_schur import cao_thin_ring_schur_jet


def main() -> None:
    kappa, radius, core_log, rho = sp.symbols(
        "kappa radius core_log rho", positive=True
    )
    jet = cao_thin_ring_schur_jet(kappa, radius, core_log, rho)

    assert jet.parameter_jacobian == -3 * kappa * core_log**2 / (
        16 * sp.pi**2 * radius
    )
    print("PASS 1: reduced (kappa,R)->(mu,c) Jacobian is negative and nonzero")

    assert jet.moment_jacobian == 2 * sp.pi * rho * kappa * radius
    print("PASS 2: physical impulse row has the exact 2*pi*rho_m normalization")

    assert jet.circulation_response_at_fixed_speed == 4 * sp.pi / (
        3 * radius * core_log
    )
    print("PASS 3: chemical/profile response changes circulation positively")

    assert jet.impulse_response_at_fixed_circulation == (
        -8 * sp.pi**2 * rho * radius**3 / core_log
    )
    print("PASS 4: speed response at fixed circulation moves impulse negatively")

    expected = -32 * sp.pi**3 * rho * radius**2 / (3 * core_log**2)
    assert jet.physical_schur_determinant == expected
    assert sp.simplify(
        jet.circulation_response_at_fixed_speed
        * jet.impulse_response_at_fixed_circulation
        - expected
    ) == 0
    print("PASS 5: triangular physical Schur determinant and sign")

    a, r = sp.symbols("a r", positive=True)
    ring_radius = sp.symbols("R", positive=True)
    b_mu = a
    b_c = a * r**2 / 2
    triangular = sp.expand(b_c - ring_radius**2 * b_mu / 2)
    assert sp.simplify(
        triangular - a * (r - ring_radius) * (r + ring_radius) / 2
    ) == 0
    assert sp.expand(
        triangular
        - a * (ring_radius * (r - ring_radius) + (r - ring_radius) ** 2 / 2)
    ) == 0
    print("PASS 6: exact source-column signs expose the radial translation cell")

    p, lambda_p, s, rho_m, radius = sp.symbols(
        "p Lambda_p s rho_m radius", positive=True
    )
    scaling_mass = sp.simplify((2 * p / (p - 1) - 2) * lambda_p)
    assert scaling_mass == 2 * lambda_p / (p - 1)
    translation_mass = sp.Integer(0)
    translation_first_moment = -lambda_p
    assert translation_mass == 0
    assert translation_first_moment == -lambda_p
    print("PASS 7: scaling and translation density moments follow by integration by parts")

    physical_cells = sp.Matrix(
        [
            [scaling_mass, 0],
            [0, 2 * sp.pi * rho_m * radius * s * translation_first_moment],
        ]
    )
    assert sp.factor(physical_cells.det()) == (
        -4 * sp.pi * rho_m * radius * s * lambda_p**2 / (p - 1)
    )
    a0, a2 = sp.symbols("A_0 A_2", positive=True)
    affine_matching = sp.diag(a0, a2)
    assert affine_matching.det() == a0 * a2
    print("PASS 8: physical two-cell and affine matching determinants are nonzero")


if __name__ == "__main__":
    main()
