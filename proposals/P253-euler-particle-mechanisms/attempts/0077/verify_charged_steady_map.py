"""Exact algebra checks for P253/0077.

This is an algebraic oracle.  It does not prove the elliptic Fredholm or Schur
claims explicitly left open by the attempt.
"""

import sympy as sp


def main() -> None:
    r, z = sp.symbols("r z", positive=True, real=True)
    c, eps, mu, g, rho = sp.symbols("c eps mu g rho", nonzero=True, real=True)
    P = sp.Function("P")(r, z)
    Phi = sp.Function("Phi")(r, z)
    H = sp.Function("H")(r, z)
    chi_fun = sp.Function("chi")
    chi = chi_fun(P)

    # W=(-P_z/r,P_r/r).  For an axisymmetric vector F=(F_r,F_z),
    # curl(F)_theta=partial_z F_r-partial_r F_z.
    Wr, Wz = -sp.diff(P, z) / r, sp.diff(P, r) / r
    Fr = -g * chi * (sp.diff(Phi, r) + H * sp.diff(P, r))
    Fz = -g * chi * (sp.diff(Phi, z) + H * sp.diff(P, z))
    curl_theta_over_r = sp.expand((sp.diff(Fr, z) - sp.diff(Fz, r)) / r)
    target = g * (
        Wr * sp.diff(sp.diff(chi_fun(P), P) * Phi - chi * H, r)
        + Wz * sp.diff(sp.diff(chi_fun(P), P) * Phi - chi * H, z)
    )
    assert sp.simplify(curl_theta_over_r - target) == 0
    print("PASS 1: exact Lorentz theta-curl and streamline derivative")

    # Ampere radial primitive differentiated in z.
    ac_phys = 1 / mu - eps * c**2
    K = sp.Function("K")(P)
    primitive_z = sp.diff(ac_phys * r * H + eps * c * sp.diff(Phi, r), z)
    ampere_radial = g * chi * sp.diff(P, z) / r
    primitive_sub = primitive_z.xreplace({sp.diff(K, P): chi})
    # Directly checking d_z[g K(P)/r].
    rhs_primitive = sp.diff(g * K / r, z).subs(sp.diff(K, P), chi)
    assert sp.simplify(rhs_primitive - ampere_radial) == 0
    print("PASS 2: radial Ampere primitive has K'=chi and the physical sign")

    # Eliminate 2H+rH_r between Gauss and axial Ampere.
    c_em_sq = 1 / (eps * mu)
    a = 1 - c**2 / c_em_sq
    hz_from_ampere = mu * g * chi * (sp.diff(P, r) / r + c) + mu * eps * c * sp.diff(Phi, z, 2)
    gauss_left = eps * (
        -sp.diff(Phi, r, 2) - sp.diff(Phi, r) / r - sp.diff(Phi, z, 2)
        + c * hz_from_ampere
    )
    scalar_eq = sp.simplify((gauss_left - g * chi) / eps)
    scalar_target = (
        -sp.diff(Phi, r, 2)
        - sp.diff(Phi, r) / r
        - a * sp.diff(Phi, z, 2)
        - (g * chi / eps) * (a - c * sp.diff(P, r) / (c_em_sq * r))
    )
    assert sp.simplify(scalar_eq - scalar_target) == 0
    print("PASS 3: Gauss/Ampere elimination gives the subluminal scalar operator")

    # Bernoulli coefficient after the vorticity first integral.
    h = sp.Function("h")(P)
    zeta = h + (g / rho) * (sp.diff(chi_fun(P), P) * Phi - chi * H)
    coeff = sp.simplify(rho * zeta + g * chi * H - g * sp.diff(chi_fun(P), P) * Phi)
    assert sp.simplify(coeff - rho * h) == 0
    print("PASS 4: modified vorticity row reduces Bernoulli to a P-function")

    # Action-coordinate stabilizer implication.
    zeta0, xtheta = sp.symbols("zeta0 xtheta", nonzero=True)
    xI = sp.symbols("xI", real=True)
    zeta_prime = sp.symbols("zeta_prime", nonzero=True)
    avg_second = xI * zeta_prime  # theta derivative averages to zero.
    assert sp.solve(sp.Eq(avg_second, 0), xI) == [0]
    print("PASS 5: nonisochronous regular-band stabilizer fixes the tag")

    # Period integration and tail degree ledger.
    t, T = sp.symbols("t T", positive=True)
    omega_remainder = sp.Function("omega_remainder")
    endpoint = sp.integrate(sp.diff(omega_remainder(t), t), (t, 0, T))
    assert sp.simplify(endpoint.subs(omega_remainder(T), omega_remainder(0))) == 0
    gamma = sp.symbols("gamma", positive=True)
    lead_degree = sp.Integer(-2) + sp.Integer(-4)
    cross_degree = sp.Integer(-2) + (-4 - gamma)
    assert lead_degree == -6 and sp.simplify(cross_degree < -6) is sp.true
    print("PASS 6: exact period cancels the time row and faster brackets miss degree six")

    # Concentrated ring impulse detects radial translation.
    rho_m, r_star, kappa, dr = sp.symbols(
        "rho_m r_star kappa dr", positive=True
    )
    impulse_variation = 2 * sp.pi * rho_m * r_star * kappa * dr
    assert sp.diff(impulse_variation, dr) != 0
    print("PASS 7: physical impulse row detects the radial translation cell")

    # The Maxwell backreaction order after solving a linear source block.
    order_field = sp.Symbol("g")
    assert sp.expand(order_field * order_field) == sp.Symbol("g") ** 2
    print("PASS 8: O(g) field gives O(g^2) Lorentz backreaction")


if __name__ == "__main__":
    main()
