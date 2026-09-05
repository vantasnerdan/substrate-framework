"""Exact scaling, solvability, and affine-action receipts for attempt 0036."""

from __future__ import annotations

import sympy as sp

from substrate_framework.verification import CheckLedger


def main() -> int:
    ledger = CheckLedger("P251-0036-smooth-polygon-algebra")
    r, d, circulation, density, total_s = sp.symbols("r d Gamma rho S", positive=True)
    p = sp.symbols("p", positive=True)
    u = sp.Function("U")(r)
    ode_u2 = -sp.diff(u, r) / r - u**p
    ode_u3 = -ode_u2 / r + sp.diff(u, r) / r**2 - p * u ** (p - 1) * sp.diff(u, r)

    def zero(name: str, expression: sp.Expr) -> None:
        ledger.check(name, sp.simplify(expression) == 0)

    scaling_mode = 2 * u / (p - 1) + r * sp.diff(u, r)
    linearized = (
        sp.diff(scaling_mode, r, 2)
        + sp.diff(scaling_mode, r) / r
        + p * u ** (p - 1) * scaling_mode
    )
    zero(
        "Lane-Emden radial scaling solves linearization",
        linearized.subs({sp.diff(u, r, 3): ode_u3, sp.diff(u, r, 2): ode_u2}),
    )
    radial_mass, constant = sp.symbols("m C", positive=True)
    exterior = constant - radial_mass * sp.log(r)
    exterior_mode = 2 * exterior / (p - 1) + r * sp.diff(exterior, r)
    zero(
        "radial scaling has nonzero logarithmic coefficient",
        r * sp.diff(exterior_mode, r) + 2 * radial_mass / (p - 1),
    )

    theta = sp.symbols("theta", real=True)
    dx, dy = d * (1 - sp.cos(theta)), -d * sp.sin(theta)
    zero(
        "each polygon vertex radial gradient",
        sp.trigsimp(dx / (dx**2 + dy**2)) - 1 / (2 * d),
    )
    force = sp.zeros(2, 1)
    for index in range(1, 6):
        angle = sp.pi * index / 3
        displacement = sp.Matrix([d * (1 - sp.cos(angle)), -d * sp.sin(angle)])
        force += (
            circulation * displacement / (2 * sp.pi * displacement.dot(displacement))
        )
    zero(
        "six-core radial force solvability",
        force[0] - 5 * circulation / (4 * sp.pi * d),
    )
    zero("six-core tangential solvability", force[1])
    omega = sp.symbols("Omega", real=True)
    divided_lambda = omega * d - force[0]
    zero("rate border derivative", sp.diff(divided_lambda, omega) - d)
    zero(
        "hexagon rate selected without imposed force",
        sp.solve(divided_lambda, omega)[0] - 5 * circulation / (4 * sp.pi * d**2),
    )

    x, core_impulse = sp.symbols("x Ccore", real=True)
    radii2 = (total_s * (1 + x) / 2, total_s * (1 - x) / 2)
    momenta = [
        -3 * density * (circulation * radius2 + core_impulse) / 2 for radius2 in radii2
    ]
    p0 = 3 * density * circulation * total_s / 4
    zero(
        "intrinsic impulse cancels relative momentum",
        (momenta[0] - momenta[1]) / 2 + p0 * x,
    )
    zero(
        "intrinsic impulse is constant in relative radius",
        sp.diff(momenta[0] + momenta[1], x),
    )

    hxx, hqq = sp.symbols("Hxx Hqq", positive=True)
    momentum, q, qdot = sp.symbols("p_rel q qdot", real=True)
    first_order = momentum * qdot - hxx * momentum**2 / (2 * p0**2) - hqq * q**2 / 2
    eliminated = first_order.subs(
        momentum, sp.solve(sp.diff(first_order, momentum), momentum)[0]
    )
    zero(
        "exact finite-core positive angular inertia",
        eliminated - (p0**2 * qdot**2 / (2 * hxx) - hqq * q**2 / 2),
    )
    point_a = 3 * density * circulation**2 / (4 * sp.pi)
    zero(
        "point limit of finite-core inertia",
        (p0**2 / hxx).subs(hxx, point_a / 2) - 3 * sp.pi * density * total_s**2 / 2,
    )
    zero(
        "point limit of finite-core gap",
        (hxx * hqq / p0**2).subs({hxx: point_a / 2, hqq: 9 * point_a / 2})
        - 9 * circulation**2 / (4 * sp.pi**2 * total_s**2),
    )
    ledger.check(
        "wrong orbital-force sign changes selected rate",
        sp.simplify(
            (omega * d + force[0]).subs(omega, 5 * circulation / (4 * sp.pi * d**2))
        )
        != 0,
    )
    print(
        "Scope: exact algebra; existence and positivity use the accompanying analytic proofs, not a numerical spectrum."
    )
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
