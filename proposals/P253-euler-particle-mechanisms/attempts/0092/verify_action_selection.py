"""Independent exact checks for P253/0092.

The script derives the finite algebra used by the attempt.  It is not a
quantization oracle and does not assert particle completion.
"""

from __future__ import annotations

import sympy as sp

from substrate_framework.euler_gauge_action_selection import (
    charge_action_map_jacobian,
    compact_bf_route_ledger,
    compact_phase_moment_charge,
    compact_u1_topology_ledger,
    conditional_self_consistent_action,
    electromagnetic_action_scale,
    gauge_field_rescaling,
    transported_tag_charge,
)


def check(label: str, condition) -> None:
    if condition is not True and condition != sp.true:
        raise AssertionError(f"FAIL {label}: {condition}")
    print(f"PASS {label}")


def main() -> None:
    # Dimension vectors use the order (mass, length, time, charge).
    dim_q = sp.Matrix([0, 0, 0, 1])
    dim_epsilon = sp.Matrix([-1, -3, 2, 2])
    dim_c = sp.Matrix([0, 1, -1, 0])
    dim_action = 2 * dim_q - dim_epsilon - dim_c
    check("electromagnetic action dimensions", dim_action == sp.Matrix([1, 2, -1, 0]))

    epsilon, mu, g_tag, a = sp.symbols(
        "epsilon mu g_tag a", positive=True
    )
    ledger = gauge_field_rescaling(epsilon, mu, g_tag, a)
    check("gauge field speed invariant", ledger.maxwell_speed_squared == 1 / (epsilon * mu))
    check("gauge charge normalization invariant", ledger.charge_action_numerator == g_tag**2 / epsilon)

    lam, c_tag = sp.symbols("lambda C_chi", real=True)
    q_lam = transported_tag_charge(g_tag, c_tag, lam)
    check("transported tag continuous family", sp.diff(q_lam, lam) == g_tag * c_tag)

    dtheta, dlambda, phi_u = sp.symbols("Dtheta Dlambda Phi_u", real=True)
    kappa_m = sp.symbols("kappa_m", nonzero=True, real=True)
    invariant_before = dtheta - kappa_m * phi_u
    invariant_after = (dtheta - kappa_m * dlambda) - kappa_m * (phi_u - dlambda)
    check("compact phase gauge cancellation", sp.expand(invariant_after - invariant_before) == 0)
    phase_number, g_0 = sp.symbols("N_phase g_0", real=True)
    check(
        "signed compact character",
        compact_phase_moment_charge(g_0, phase_number, -1)
        == -g_0 * phase_number,
    )

    topo = compact_u1_topology_ledger()
    check("whole space U1 topology", (topo.pi_3_u1, topo.h2_s3, topo.whole_space_magnetic_flux) == (0, 0, 0))
    check("sphere infinity Chern group", topo.h2_s2 == "Z")

    action, speed = sp.symbols("J c", positive=True)
    jac = charge_action_map_jacobian(g_tag, c_tag, lam, action, epsilon, speed)
    check("fixed charge action direction", jac[:, 1] == sp.Matrix([0, 1, 0]))

    root_plus = conditional_self_consistent_action(kappa_m, epsilon, speed)
    root_minus = conditional_self_consistent_action(-kappa_m, epsilon, speed)
    equation = sp.factor(
        action
        - electromagnetic_action_scale(kappa_m * action, epsilon, speed)
    )
    check(
        "signed charge-per-action coefficient same nonzero root",
        root_plus == root_minus,
    )
    check("imposed nonzero root", sp.simplify(equation.subs(action, root_plus)) == 0)

    bf = compact_bf_route_ledger()
    check(
        "Euler vorticity is not compact BF data",
        not bf.euler_vorticity_is_compact_connection
        and not bf.euler_vorticity_periods_are_integral
        and bf.compact_two_form_is_added_structure
        and not bf.level_is_euler_selected,
    )


if __name__ == "__main__":
    main()
