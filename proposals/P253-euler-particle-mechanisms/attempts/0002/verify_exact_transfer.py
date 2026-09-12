"""Exact checks for the P253/0002 source-transfer calculations."""

import sympy as sp


def main() -> None:
    lam, a, rho, theta, r0 = sp.symbols(
        "lambda a rho theta r0", positive=True, real=True
    )
    r, z = sp.symbols("r z", real=True)

    # Hill ball: dx = r dr dtheta dz and 0 <= r <= sqrt(a^2-z^2).
    gamma = sp.integrate(
        lam * r, (r, 0, sp.sqrt(a**2 - z**2)), (z, -a, a), (theta, 0, 2 * sp.pi)
    )
    impulse = sp.Rational(1, 2) * sp.integrate(
        lam * r**3,
        (r, 0, sp.sqrt(a**2 - z**2)),
        (z, -a, a),
        (theta, 0, 2 * sp.pi),
    )
    assert sp.simplify(gamma - 4 * sp.pi * lam * a**3 / 3) == 0
    assert sp.simplify(impulse - 4 * sp.pi * lam * a**5 / 15) == 0

    # Davila et al. k=2, q1=-q2=q. Convention q^perp=(-q2,q1).
    x, y = sp.symbols("x y", real=True)
    radius2 = x**2 + y**2
    hamiltonian = -4 * sp.log(2 * sp.sqrt(radius2)) - 2 * x**2 / r0**2
    dx = 2 * y / radius2
    dy = -2 * x / radius2 - 2 * x / r0**2
    hdot = sp.simplify(sp.diff(hamiltonian, x) * dx + sp.diff(hamiltonian, y) * dy)
    polar = {x: rho * sp.cos(theta), y: rho * sp.sin(theta)}
    drho = sp.trigsimp(((x * dx + y * dy) / sp.sqrt(radius2)).subs(polar))
    dtheta = sp.trigsimp(((x * dy - y * dx) / radius2).subs(polar))
    hamiltonian_polar = sp.simplify(hamiltonian.subs(polar))
    assert hdot == 0
    assert sp.trigsimp(drho + rho * sp.sin(2 * theta) / r0**2) == 0
    assert sp.trigsimp(dtheta + 2 / rho**2 + 2 * sp.cos(theta) ** 2 / r0**2) == 0
    assert sp.simplify(sp.diff(hamiltonian_polar, rho)) == -4 / rho - 4 * rho * sp.cos(theta) ** 2 / r0**2

    # Bilinear, self-adjoint Green operator: exact two-copy energy expansion.
    e11, e22, e12, e21 = sp.symbols("e11 e22 e12 e21", real=True)
    expanded = sp.Rational(1, 2) * (e11 + e12 + e21 + e22)
    expected = e11 / 2 + e22 / 2 + e12
    assert sp.simplify(expanded.subs(e21, e12) - expected) == 0

    print("Hill circulation =", gamma)
    print("Hill impulse =", impulse)
    print("reduced Hamiltonian derivative =", hdot)
    print("rho_dot =", drho)
    print("theta_dot =", dtheta)
    print("dH/drho =", sp.diff(hamiltonian_polar, rho))
    print("two-copy energy = E1 + E2 + <zeta1,G zeta2>")
    print("ALL 7 EXACT CHECKS PASS")


if __name__ == "__main__":
    main()
