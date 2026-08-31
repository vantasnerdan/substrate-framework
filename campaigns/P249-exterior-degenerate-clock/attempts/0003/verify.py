"""Complete fixed-J calculus and commuting-route virials for P249."""

from __future__ import annotations

import sympy as sp

from substrate_framework.verification import CheckLedger


def main() -> int:
    ledger = CheckLedger("P249-O3-0003")
    angular_momentum = sp.Symbol("J", nonzero=True, real=True)
    inertia = sp.Symbol("C", positive=True)
    delta_c, delta2_c = sp.symbols("delta_C delta2_C", real=True)
    delta_e0, delta2_e0 = sp.symbols("delta_E0 delta2_E0", real=True)

    rotational = angular_momentum**2 / (4 * inertia)
    omega = angular_momentum / (2 * inertia)
    ledger.check("charge-frequency map exact", sp.diff(rotational, angular_momentum) == omega)
    ledger.check("first inertia derivative exact", sp.diff(rotational, inertia) == -omega**2)
    ledger.check(
        "second inertia derivative exact",
        sp.diff(rotational, inertia, 2) == 2 * omega**2 / inertia,
    )
    first_variation = delta_e0 + sp.diff(rotational, inertia) * delta_c
    second_variation = (
        delta2_e0
        + sp.diff(rotational, inertia) * delta2_c
        + sp.diff(rotational, inertia, 2) * delta_c**2
    )
    ledger.check("complete fixed-J first variation", sp.factor(first_variation - (delta_e0 - omega**2 * delta_c)) == 0)
    ledger.check(
        "complete fixed-J Hessian rank-one term",
        sp.factor(second_variation - (delta2_e0 - omega**2 * delta2_c + 2 * omega**2 * delta_c**2 / inertia)) == 0,
    )
    ledger.check(
        "dropping rank-one term changes the Hessian",
        sp.factor(second_variation - (delta2_e0 - omega**2 * delta2_c)) != 0,
    )

    scale = sp.Symbol("R", positive=True)
    e4, e2, e0, c1, c3 = sp.symbols("E4 E2 E0 C1 C3", nonnegative=True)
    static_scaled = e4 / scale + e2 * scale + e0 * scale**3
    inertia_scaled = c1 * scale + c3 * scale**3
    fixed_scaled = static_scaled + angular_momentum**2 / (4 * inertia_scaled)
    scale_derivative = sp.factor(scale * sp.diff(fixed_scaled, scale)).subs(scale, 1)
    expected_virial = -e4 + e2 + 3 * e0 - omega**2 * (c1 + 3 * c3)
    expected_at_scaled_c = sp.factor(expected_virial.subs(inertia, c1 + c3))
    ledger.check("complete spatial Derrick identity", sp.factor(scale_derivative - expected_at_scaled_c) == 0)

    amplitude = sp.Symbol("a", positive=True)
    p2, p4, i4 = sp.symbols("P2 P4 I4", positive=True)
    radius = sp.Symbol("R_c", positive=True)
    u2 = p2 * amplitude**2 * radius**3
    u4 = p4 * amplitude**4 * radius**3
    e_rot = angular_momentum**2 / (4 * i4 * amplitude**4 * radius)
    commuting_energy = u2 + u4 + e_rot
    radial_virial = sp.factor(radius * sp.diff(commuting_energy, radius))
    amplitude_virial = sp.factor(amplitude * sp.diff(commuting_energy, amplitude))
    ledger.check("commuting radial virial exact", sp.factor(radial_virial - (3 * (u2 + u4) - e_rot)) == 0)
    ledger.check("commuting amplitude virial exact", sp.factor(amplitude_virial - (2 * u2 + 4 * u4 - 4 * e_rot)) == 0)
    eliminated = sp.factor((2 * u2 + 4 * u4 - 4 * 3 * (u2 + u4)))
    ledger.check("virial elimination exact", sp.factor(eliminated + 10 * u2 + 8 * u4) == 0)
    ledger.check("virial elimination is strictly negative", eliminated.is_negative)

    collapse_radius = amplitude**-2
    collapsed = sp.factor(commuting_energy.subs(radius, collapse_radius))
    expected_collapse = p2 / amplitude**4 + p4 / amplitude**2 + angular_momentum**2 / (4 * i4 * amplitude**2)
    ledger.check("branch-preserving collapse formula exact", sp.factor(collapsed - expected_collapse) == 0)
    ledger.check("branch-preserving collapse tends to zero", sp.limit(collapsed, amplitude, sp.oo) == 0)

    d = sp.Symbol("d", nonnegative=True)
    tau = 4 + 2 * d
    ledger.check("timelike branch stays above positive split", sp.simplify(tau - d) > 0)
    ledger.check("timelike branch stays above negative split", sp.simplify(tau + d) > 0)
    ledger.check("timelike branch stays above director eigenvalue", sp.simplify(tau - 1) > 0)

    lock = sp.Symbol("P_lock", positive=True) * amplitude**2 * radius**3
    locked_u2 = u2 + lock
    locked_elimination = sp.factor(-10 * locked_u2 - 8 * u4)
    ledger.check("positive axis lock cannot solve commuting virials", locked_elimination.is_negative)

    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
