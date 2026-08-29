"""Exact block-map and action-pullback checks for P248."""

from __future__ import annotations

import sympy as sp

from substrate_framework.optical_gothic import (
    determinant_slaved_optical_metric,
    optical_adm_metric,
    recover_optical_adm,
)
from substrate_framework.verification import CheckLedger


def run() -> int:
    ledger = CheckLedger("P248-METRIC")
    lapse, speed = sp.symbols("N c", positive=True)
    g1, g2, g3 = sp.symbols("g1 g2 g3", positive=True)
    v1, v2, v3 = sp.symbols("v1 v2 v3", real=True)
    spatial = sp.diag(g1, g2, g3)
    flow = speed * sp.Matrix([v1, v2, v3])

    complete = optical_adm_metric(lapse, spatial, flow, speed)
    ledger.check(
        "complete block inverse",
        (complete.covariant * complete.contravariant).applyfunc(sp.simplify)
        == sp.eye(4),
    )
    ledger.check(
        "flow-independent determinant",
        sp.simplify(complete.determinant + lapse**2 * g1 * g2 * g3) == 0,
    )
    ledger.check(
        "Lorentzian determinant sign",
        complete.determinant.is_negative is True,
    )
    ledger.check(
        "nonzero complete-map Jacobian",
        complete.component_jacobian_determinant.is_positive is True,
    )

    # Derive the ten-component point-map Jacobian independently of the API.
    q = sp.symbols("N v1 v2 v3 a b c d e f", real=True)
    N, x, y, z, a, b, c, d, e, f = q
    gamma = sp.Matrix([[a, d, e], [d, b, f], [e, f, c]])
    velocity = sp.Matrix([x, y, z])
    mixed = -gamma * velocity
    outputs = sp.Matrix(
        [
            -N**2 + (velocity.T * gamma * velocity)[0],
            mixed[0],
            mixed[1],
            mixed[2],
            a,
            b,
            c,
            d,
            e,
            f,
        ]
    )
    jacobian = outputs.jacobian(q)
    ledger.check(
        "independent component Jacobian derivation",
        sp.factor(jacobian.det() - 2 * N * gamma.det()) == 0,
    )
    metric_euler = sp.Matrix(sp.symbols("E0:10", real=True))
    optical_euler = jacobian.T * metric_euler
    recovered_euler = jacobian.T.inv() * optical_euler
    ledger.check(
        "invertible variational chain rule",
        (recovered_euler - metric_euler).applyfunc(sp.simplify) == sp.zeros(10, 1),
    )
    ledger.check(
        "dropping the lapse loses one variational equation",
        jacobian[:, 1:].rank() == 9,
    )

    witness = optical_adm_metric(
        2,
        [[2, 1, 0], [1, 3, 1], [0, 1, 4]],
        [3, -1, 2],
        5,
    )
    recovered = recover_optical_adm(witness.covariant, 5)
    ledger.check(
        "complete map exact round trip",
        recovered.reconstructed_metric.covariant == witness.covariant,
    )

    slaved = determinant_slaved_optical_metric(sp.eye(3), [0, 0, 0], 1)
    ledger.check(
        "nine-field image constraint",
        slaved.lapse_constraint_residual == 0
        and slaved.volume_constraint_residual == 0,
    )
    independent_lapse_metric = optical_adm_metric(
        2, sp.eye(3), [0, 0, 0], 1
    ).covariant
    ledger.check(
        "independent-lapse counterexample",
        independent_lapse_metric != slaved.metric.covariant,
        detail=str(independent_lapse_metric),
    )
    expected_gothic = sp.diag(-1, 1, 1, 1)
    ledger.check(
        "undeformed gothic limit",
        slaved.metric.gothic_contravariant == expected_gothic,
    )
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(run())
