"""Exact exposing algebra for the 0142 analytic Lundquist packet construction.

This checks derived identities, not a discretized spectrum or a tally proxy
for the finite-time theorem. Its analytic remainder is in the two proof files.
"""

from __future__ import annotations

import sympy as sp

from substrate_framework.verification import CheckLedger


def main() -> int:
    ledger = CheckLedger("P251-0142")
    x, r, w = sp.symbols("x r w", positive=True)
    delta, omega, lam, p = sp.symbols("delta Omega lambda p", positive=True)
    m = 2
    laguerre = sp.assoc_laguerre(2, m - 1, x)
    profile = r ** (m - 1) * sp.exp(-r**2 / 2) * laguerre.subs(x, r**2)
    oscillator = (
        -sp.diff(profile, r, 2)
        - sp.diff(profile, r) / r
        + (m - 1) ** 2 * profile / r**2
        + r**2 * profile
    )
    ledger.check("full radial oscillator eigenvalue12", sp.simplify(oscillator - 12 * profile) == 0)
    ledger.check("wrong oscillator eigenvalue rejected", sp.simplify(oscillator - 10 * profile) != 0)
    norm = sp.integrate(x * sp.exp(-x) * laguerre**2, (x, 0, sp.oo))
    ledger.check("positive n2 Laguerre KKS norm", norm == 3)

    ell = (2 / (lam * p**3)) ** sp.Rational(1, 4)
    pressure_scale = sp.simplify(1 / (p**2 * ell**2))
    ledger.check("nonlocal Leray small parameter", pressure_scale == sp.sqrt(lam / p) / sp.sqrt(2))
    doppler_scale = sp.simplify((2 * omega / lam) * lam**2 * p * ell**2 / 4)
    ledger.check("axial Doppler trap scale", doppler_scale == omega * sp.sqrt(lam / p) / sp.sqrt(2))
    numerator_series = sp.series(sp.Rational(1, 2) - sp.besselj(1, x) / x, x, 0, 6).removeO()
    denominator_series = sp.series(1 - sp.besselj(0, x), x, 0, 6).removeO()
    axis_ratio = sp.limit(numerator_series / denominator_series, x, 0)
    ledger.check("no-critical-radius ratio regular at axis", axis_ratio == sp.Rational(1, 4))

    s, az, wz = sp.symbols("s az wz", real=True, nonzero=True)
    nilpotent = sp.Matrix([[0, az, wz], [0, 0, 0], [0, 0, 0]])
    covector = -sp.I * s * sp.eye(3) + nilpotent
    inv = sp.I / s * sp.eye(3) + nilpotent / s**2
    ledger.check("full covector Kelvin inverse", sp.simplify(covector * inv - sp.eye(3)) == sp.zeros(3))
    vr, pres, zeta = sp.symbols("vr P Z")
    xi_r = sp.I * vr / s
    v_theta = m * pres / (r * s) - sp.I * zeta * vr / s
    spin = sp.simplify(r * (v_theta + zeta * xi_r))
    ledger.check("full material-position spin equals mP/s", spin == m * pres / s)
    wrong_spin = sp.simplify(r * (v_theta + (zeta - az) * xi_r))
    ledger.check("omitted angular-shear position term exposed", sp.simplify(wrong_spin - spin) == -sp.I * az * r * vr / s)

    # Both opposite-helicity pressure and gap are derived in the same
    # Cartesian convention J e_plus=-i e_plus.
    derivative = sp.Symbol("Dplus2F")
    opposite_pressure = -sp.I * omega * derivative / p**2
    opposite_velocity = sp.solve(4 * sp.I * omega * sp.Symbol("vminus") + opposite_pressure, sp.Symbol("vminus"))[0]
    ledger.check("opposite-helicity pressure correction", opposite_velocity == derivative / (4 * p**2))

    angle_integral = sp.integrate(x**2 * sp.exp(-w * x / 2) * laguerre, (x, 0, sp.oo))
    normalized_angle = sp.factor(angle_integral / angle_integral.subs(w, 1))
    ledger.check("exact excited transported quadrupole", normalized_angle == (w - 4) * (w - 2) / (3 * w**5))
    angle_mean = sp.simplify(
        sp.integrate(x**3 * sp.exp(-x / 2) * laguerre, (x, 0, sp.oo))
        / angle_integral.subs(w, 1)
    )
    ledger.check("signed physical angle first Doppler moment", angle_mean == sp.Rational(38, 3))
    phase_slope = sp.simplify(sp.diff(normalized_angle, w).subs(w, 1))
    ledger.check("physical phase correction -19/3", phase_slope == -sp.Rational(19, 3))
    sigma = -2 * omega + 6 * sp.sqrt(2) * omega * delta
    gamma = sp.expand(sigma + phase_slope * sp.sqrt(2) * omega * delta)
    leading_square = sp.series(gamma**2, delta, 0, 2).removeO()
    curvature = sp.expand(delta**2 * sp.diff(leading_square, delta, 2) / 4 + 3 * delta * sp.diff(leading_square, delta) / 4)
    ledger.check("positive physical natural-scale curvature", curvature == sp.sqrt(2) * omega**2 * delta)
    intrinsic = sp.series(sigma**2, delta, 0, 2).removeO()
    intrinsic_curvature = sp.expand(delta**2 * sp.diff(intrinsic, delta, 2) / 4 + 3 * delta * sp.diff(intrinsic, delta) / 4)
    ledger.check("intrinsic curvature remains negative", intrinsic_curvature == -18 * sp.sqrt(2) * omega**2 * delta)

    # The polynomial row transformation includes the exponential derivative.
    def carrier_row(poly: sp.Expr) -> sp.Expr:
        return sp.expand(-poly / 2 + sp.Rational(3, 2) * x * (sp.diff(poly, x) - poly / 2))

    matrices = []
    for n in (1, 2):
        ln = sp.assoc_laguerre(n, 1, x)
        spin_poly = sp.expand(ln - 2 * sp.diff(ln, x))
        rows = []
        for first in (spin_poly, x * spin_poly):
            rows += [first, carrier_row(first), carrier_row(carrier_row(first))]
        matrix = sp.Matrix([[sp.expand(row).coeff(x, j) for j in range(6)] for row in rows])
        matrices.append(matrix)
    ledger.check("first excitation still rank deficient", matrices[0].rank() == 5)
    determinant = sp.factor(matrices[1].det())
    ledger.check("second excitation repairs all six rows", determinant == sp.Rational(19683, 8192))
    ledger.check("removing transported-time row destroys rank", matrices[1][:5, :].rank() == 5)
    extended = [sp.Integer(1)] + [x**j * sp.exp(-x / 2) for j in range(6)]
    wronskian = sp.simplify(sp.wronskian(extended, x))
    ledger.check("seven-row physical-marker Wronskian nonzero", wronskian != 0 and sp.simplify(wronskian * sp.exp(3 * x)).is_nonzero)

    # Exact scalar action Euler--Lagrange equation, with all time connections.
    t = sp.Symbol("t", real=True)
    theta = sp.Function("theta")(t)
    c = sp.Function("c")(t)
    gamma_t = sp.Function("gamma")(t)
    beta = sp.Symbol("beta", positive=True)
    connection = sp.diff(c, t) / c
    mass = -beta / (gamma_t * c**2)
    action = mass * ((sp.diff(theta, t) - connection * theta) ** 2 - gamma_t**2 * theta**2) / 2
    el = sp.simplify((sp.diff(sp.diff(action, sp.diff(theta, t)), t) - sp.diff(action, theta)) / mass)
    expected = (
        sp.diff(theta, t, 2)
        + (-sp.diff(gamma_t, t) / gamma_t - 2 * connection) * sp.diff(theta, t)
        + (gamma_t**2 + connection**2 - sp.diff(connection, t) + sp.diff(gamma_t, t) * connection / gamma_t) * theta
    )
    ledger.check("physical scalar action retains both time connections", sp.simplify(el - expected) == 0)
    ledger.check("constant-row negative-frequency mass positive", mass.subs({c: 2, gamma_t: -omega}) == beta / (4 * omega))
    print(f"n2 marker determinant={determinant}; physical p² curvature={curvature}")
    print(f"seven-row Wronskian={wronskian}; all arithmetic exact SymPy")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
