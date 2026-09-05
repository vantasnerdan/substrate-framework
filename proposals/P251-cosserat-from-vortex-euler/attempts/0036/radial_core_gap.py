"""Exact identities supporting radial-core-gap.md; no discretized spectrum."""

from __future__ import annotations

import sympy as sp

from substrate_framework.verification import CheckLedger


def main() -> int:
    ledger = CheckLedger("P251-0036-radial-core-identities")
    r, s, radius, c = sp.symbols("r s R c", positive=True)
    m = sp.symbols("m", integer=True, positive=True)
    w = sp.Function("w")(r)
    mass = sp.Function("M")(r)
    substitutions = {
        sp.diff(mass, r, 2): w + r * sp.diff(w, r),
        sp.diff(mass, r): r * w,
    }

    def zero(name: str, expression: sp.Expr) -> None:
        ledger.check(name, sp.simplify(expression) == 0)

    def laplacian_mode(expression: sp.Expr, order: sp.Expr) -> sp.Expr:
        return (
            sp.diff(expression, r, 2)
            + sp.diff(expression, r) / r
            - order**2 * expression / r**2
        )

    omega = mass / r**2
    zero(
        "angular-velocity derivative",
        sp.diff(omega, r).subs(substitutions) - (w - 2 * omega) / r,
    )
    zero(
        "translation Poisson identity",
        -laplacian_mode(mass / r, 1).subs(substitutions) + sp.diff(w, r),
    )
    zero(
        "Doob normalization integration-by-parts integrand",
        -mass * sp.diff(w, r) - (r * w**2 - sp.diff(mass * w, r)).subs(substitutions),
    )

    kernel_inside = (r / s) ** m / (2 * m)
    kernel_outside = (s / r) ** m / (2 * m)
    zero("Green kernel interior harmonic", laplacian_mode(kernel_inside, m))
    zero("Green kernel exterior harmonic", laplacian_mode(kernel_outside, m))
    zero("Green kernel continuity", (kernel_inside - kernel_outside).subs(r, s))
    zero(
        "Green derivative jump uses radial measure",
        (sp.diff(kernel_outside, r) - sp.diff(kernel_inside, r)).subs(r, s) + 1 / s,
    )
    zero(
        "mode comparison interior",
        m * kernel_inside / kernel_inside.subs(m, 1) - (r / s) ** (m - 1),
    )
    zero(
        "mode comparison exterior",
        m * kernel_outside / kernel_outside.subs(m, 1) - (s / r) ** (m - 1),
    )

    ar, ass, omr, oms, z = sp.symbols("a_r a_s Omega_r Omega_s Z", positive=True)
    h_r = sp.sqrt(omr * ar) * r
    h_s = sp.sqrt(oms * ass) * s
    b_inside = sp.sqrt(ar / omr) * kernel_inside.subs(m, 1) * sp.sqrt(ass / oms)
    p_density = b_inside * h_s / h_r
    pi_density_radial_measure = oms * ass * s**2 / z
    zero(
        "Doob kernel ratio when r<=s",
        p_density / pi_density_radial_measure - z / (2 * omr * oms * s**2),
    )

    bump = c * sp.exp(-1 / (1 - r**2 / radius**2))
    zero("smooth bump central value", bump.subs(r, 0) - c / sp.E)
    zero(
        "smooth bump half-area lower value",
        bump.subs(r, radius / sp.sqrt(2)) - c / sp.E**2,
    )
    z_lower = sp.integrate(r * c**2 / sp.E**4, (r, 0, radius / sp.sqrt(2)))
    omega_max = c / (2 * sp.E)
    zero(
        "explicit minorization bound",
        z_lower / (2 * omega_max**2 * radius**2) - 1 / (2 * sp.E**2),
    )
    zero(
        "bump profile slope is nonpositive",
        -sp.diff(bump, r) / r
        - 2
        * c
        * sp.exp(-1 / (1 - r**2 / radius**2))
        / (radius**2 * (1 - r**2 / radius**2) ** 2),
    )

    omega_boundary = sp.symbols("Omega_boundary", positive=True)
    rankine_frequency = m * omega_boundary * (1 - 1 / m)
    zero("Rankine comparison translation frequency", rankine_frequency.subs(m, 1))
    zero(
        "Rankine comparison Kelvin frequency",
        rankine_frequency - (m - 1) * omega_boundary,
    )
    ledger.check(
        "wrong Green sign destroys translation identity",
        sp.simplify(laplacian_mode(mass / r, 1).subs(substitutions) + sp.diff(w, r))
        != 0,
    )
    ledger.check(
        "missing Green radial measure changes jump",
        sp.simplify(
            (sp.diff(kernel_outside, r) - sp.diff(kernel_inside, r)).subs(r, s) + 1
        )
        != 0,
    )
    print(
        "Scope: algebraic receipts only; positivity and Hilbert-space proof are in radial-core-gap.md."
    )
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
