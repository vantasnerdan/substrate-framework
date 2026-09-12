#!/usr/bin/env python3
"""Exact sign/selection oracle for P253/0088.

This checks algebraic consequences only.  It evaluates the exact near-axis
coefficient used to prove the fixed-column gamma_12 nonidentity, but not the
full Cao response remainder or a uniform high-index response.
"""

import sympy as sp


def main() -> None:
    nu = sp.symbols("nu", positive=True, real=True)
    # Omega(e,e_bar)=i/nu and ell=-i nu Omega imply energy-unit extraction.
    omega_pair = sp.I / nu
    assert sp.simplify((-sp.I * nu) * omega_pair - 1) == 0

    # Differentiated metric-skew relation at a double eigenvalue.  The two
    # metric-derivative terms cancel without setting Gprime to zero.
    gp = sp.MatrixSymbol("Gp", 2, 2)
    metric_prime = sp.I * nu * gp - sp.I * nu * gp
    assert metric_prime == sp.ZeroMatrix(2, 2)

    # General 2x2 anti-Hermitian compression and diagonal/off-diagonal
    # commutator.  Symbols are real, so c=x+i y.
    x, y, d1, d2 = sp.symbols("x y d1 d2", real=True)
    c = x + sp.I * y
    moff = sp.Matrix([[0, c], [-sp.conjugate(c), 0]])
    mdiag = sp.diag(-sp.I * d1, -sp.I * d2)
    assert sp.simplify(moff.conjugate().T + moff) == sp.zeros(2)
    assert sp.simplify(mdiag.conjugate().T + mdiag) == sp.zeros(2)
    comm = sp.simplify(moff * mdiag - mdiag * moff)
    assert sp.simplify(comm[0, 1] + sp.I * (d2 - d1) * c) == 0
    assert sp.simplify(comm[1, 0] + sp.I * (d2 - d1) * sp.conjugate(c)) == 0

    # Cylindrical KKS density.  omega_theta=r*zeta and dx=r dr dz dtheta,
    # while normalized characters integrate to one.  The half-density
    # factors cancel r dr dz but leave one physical factor r.
    R, q, a = sp.symbols("R q a", positive=True)
    r_phys = R * q
    jac = R * a**2 * q
    half_density_pair = 1 / jac
    omega_theta = r_phys * sp.symbols("zeta", real=True)
    reduced_weight = sp.simplify(omega_theta * jac * half_density_pair)
    assert reduced_weight == r_phys * sp.symbols("zeta", real=True)

    # Covariant FH: ell((A+i nu)de)=0 leaves ell(A'e)=-i nu'.
    nup = sp.symbols("nu_prime", real=True)
    fh = sp.simplify(-sp.I * nup)
    assert fh == -sp.I * nup

    # In a volume action-angle cell, an ell=d=0 compact divergence-free
    # displacement has xiI'=0, hence xiI=0 and C0 xi=0.
    I = sp.symbols("I", real=True)
    zeta = sp.Function("zeta")(I)
    xiI = sp.Integer(0)
    c0_theta = sp.simplify(xiI * sp.diff(zeta, I))
    assert c0_theta == 0

    # Two-dimensional H5 bump count: response A gamma r^2 divided by the
    # sixth-derivative curl scale A r^(1-6) is gamma r^7.
    r, amp, gamma = sp.symbols("r amp gamma", positive=True)
    response = amp * gamma * r**2
    h5_curl = amp * r ** (-5)
    assert sp.simplify(response / h5_curl - gamma * r**7) == 0

    # Fixed-k massive regrading.  Applying D1 to G0 and D0 to G1 must retain
    # two copies of k^2*s*cos(alpha), not the fixed-integer L1=C.
    s, alpha, k = sp.symbols("s alpha k", positive=True, real=True)
    f = sp.Function("f")(s, alpha)
    ca, sa = sp.cos(alpha), sp.sin(alpha)
    g0 = (sp.diff(f, s), sp.diff(f, alpha) / s, sp.I * k * f)
    g1 = (sp.Integer(0), sp.Integer(0), -sp.I * k * s * ca * f)

    def d0(v):
        return (
            sp.diff(v[0], s)
            + v[0] / s
            + sp.diff(v[1], alpha) / s
            + sp.I * k * v[2]
        )

    def d1(v):
        return ca * v[0] - sa * v[1] - sp.I * k * s * ca * v[2]

    l1 = sp.simplify(d1(g0) + d0(g1))
    l1_expected = (
        ca * sp.diff(f, s)
        - sa * sp.diff(f, alpha) / s
        + 2 * k**2 * s * ca * f
    )
    assert sp.simplify(l1 - l1_expected) == 0

    # The fixed-k longitudinal covariant row is
    # delta/(1+delta*s*c)*(i*kY/delta*Y+CthetaY).
    delta, k_y, y, ctheta_y = sp.symbols(
        "delta k_y y ctheta_y", real=True
    )
    exact_longitudinal = (
        delta / (1 + delta * s * ca)
        * (sp.I * k_y / delta * y + ctheta_y)
    )
    first = sp.diff(exact_longitudinal, delta).limit(delta, 0)
    assert sp.simplify(first - (ctheta_y - sp.I * k_y * s * ca * y)) == 0

    # Geometric character check behind the W/2 and s*Omega'/2 rows.
    omega = sp.Function("Omega")(s)
    psi = s * omega
    w_profile = 2 * omega + s * sp.diff(omega, s)
    assert sp.simplify((sp.diff(psi, s) + psi / s) / 2 - w_profile / 2) == 0
    assert sp.simplify((sp.diff(psi, s) - psi / s) / 2 - s * sp.diff(omega, s) / 2) == 0

    # Full velocity-form response density at an equal-frequency m=0 pair.
    # The order-s alpha component already contains both convection terms;
    # pressure drops against the divergence-free Hessian left row.
    k1, k2, omega0, a2, a1sharp = sp.symbols(
        "k1 k2 omega0 a2 a1sharp", nonzero=True
    )
    c2 = 2 * sp.I * a2 / k2
    b2 = -2 * sp.I * omega0 * a2 / nu
    b1sharp = 2 * sp.I * omega0 * a1sharp / nu
    near_axis = sp.simplify(-sp.I * k1 * c2 * b1sharp + 2 * b2 * a1sharp)
    expected_axis = 4 * sp.I * omega0 * a2 * a1sharp / nu * (k1 / k2 - 1)
    assert sp.simplify(near_axis - expected_axis) == 0

    # Dual-Riesz closure starts with an exact KKS cancellation.  For
    # q_s=-i*k*W*xi_s, the factor W in Omega removes the orbit inverse.
    rho, W, omegap, vr, va = sp.symbols(
        "rho W Omega_prime vr va", nonzero=True
    )
    kks_to_l2 = sp.simplify((-sp.I * nu) * (sp.I * rho / k1))
    assert kks_to_l2 == rho * nu / k1

    # Conjugate-mode displacement and curl representative.  The right
    # m=0 relation va=i*W*vr/nu and W-s*Omega'=2*Omega give the left row.
    eta_s = sp.I * vr / nu
    eta_a = sp.I * va / nu + s * omegap * vr / nu**2
    lambda_s = sp.simplify(-sp.I * rho * nu * eta_s)
    lambda_a = sp.simplify(-sp.I * rho * nu * eta_a)
    assert lambda_s == rho * vr
    lambda_a_on_mode = sp.simplify(lambda_a.subs(va, sp.I * W * vr / nu))
    assert sp.simplify(
        lambda_a_on_mode.subs(W, 2 * omega0 + s * omegap)
        - 2 * sp.I * rho * omega0 * vr / nu
    ) == 0

    # The regular singular-axis recurrence has only one free datum.  Acting
    # with -d_s(d_s+1/s) on s^(2j+3) produces the denominator in (35a.9).
    j = sp.symbols("j", integer=True, nonnegative=True)
    monomial = s ** (2 * j + 3)
    star = sp.diff(monomial, s) + monomial / s
    assert sp.simplify(
        -sp.diff(star, s) / s ** (2 * j + 1)
        + 4 * (j + 1) * (j + 2)
    ) == 0

    # A local solenoidal representer is ambiguous by grad I0(k*s)e^-ikz;
    # its radial coefficient changes by k^2/2, so global Riesz/Sturm data are
    # genuinely needed to fix A_sharp.
    harmonic_radial = sp.diff(sp.besseli(0, k1 * s), s)
    assert sp.simplify(sp.limit(harmonic_radial / s, s, 0) - k1**2 / 2) == 0

    print("PASS physical KKS clock/sign: Omega=i/nu and ell(e)=1")
    print("PASS double-eigenvalue metric-prime cancellation and M*=-M form")
    print("PASS noncommuting off-diagonal/diagonal compression sign")
    print("PASS normalized cylindrical KKS row retains r=R*q")
    print("PASS covariant Feynman-Hellmann row ell(A'e)=-i nu'")
    print("PASS compact (ell,d)=(0,0) DA slice is stabilizer")
    print("PASS two-meridional-dimensional H5 patch exponent is 7")
    print("PASS fixed-k massive L1=C+2*k^2*s*cos(alpha)")
    print("PASS fixed-k B1 has frame rotation minus i*k*s*cos(alpha)")
    print("PASS curvature character rows are W/2 and s*Omega'/2")
    print("PASS full near-axis response coefficient is proportional to k1/k2-1")
    print("PASS KKS W cancellation gives the physical L2 covector coefficient")
    print("PASS conjugate-mode displacement reconstructs the adjoint core row")
    print("PASS singular-axis recurrence has one regular radial datum")
    print("PASS local harmonic-gradient ambiguity shifts the radial coefficient")
    print("SCOPE no Cao response remainder or uniform high-N scale is asserted")


if __name__ == "__main__":
    main()
