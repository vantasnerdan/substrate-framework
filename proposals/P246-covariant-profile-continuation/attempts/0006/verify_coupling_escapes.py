"""P246 attempt 0006: exact xi and baseline compactness boundaries."""

from __future__ import annotations

import contextlib
import io
import json
import traceback
from pathlib import Path

import sympy as sp

from substrate_framework.verification import CheckLedger

HERE = Path(__file__).resolve().parent
G_ACCEPTED = sp.Float("46.80699908016004", 50)
G_CRITICAL = sp.Float("0.0523678366438565", 50)


def execute() -> None:
    ledger = CheckLedger("P246-attempt-0006-coupling-escapes")
    xi, baseline, delta_zero = sp.symbols("xi B Delta_0", real=True)
    inverse_total = baseline + delta_zero * (1 - 6 * xi)
    gravity = 1 / inverse_total
    delta_value = 1 / G_ACCEPTED
    critical_inverse = 1 / G_CRITICAL
    baseline_boundary = sp.N(critical_inverse - delta_value, 30)
    xi_boundary = sp.N((1 - critical_inverse / delta_value) / 6, 30)
    boundary_gravity_baseline = sp.N(1 / (baseline_boundary + delta_value), 40)
    boundary_gravity_xi = sp.N(1 / (delta_value * (1 - 6 * xi_boundary)), 40)
    conformal_inverse = sp.simplify(
        inverse_total.subs({baseline: 0, xi: sp.Rational(1, 6)})
    )

    ledger.check(
        "baseline_boundary_reconstructs_critical_G",
        abs(float(boundary_gravity_baseline - G_CRITICAL)) / float(G_CRITICAL)
        < 1.0e-13,
    )
    ledger.check(
        "xi_boundary_reconstructs_critical_G",
        abs(float(boundary_gravity_xi - G_CRITICAL)) / float(G_CRITICAL) < 1.0e-13,
    )
    ledger.check(
        "baseline_is_bare_dominated",
        baseline_boundary / delta_value > 800,
        f"B/Delta0={float(baseline_boundary / delta_value):.6f}",
    )
    ledger.check(
        "reduced_G_requires_large_negative_xi",
        xi_boundary < -100,
        f"xi_boundary={float(xi_boundary):.12g}",
    )
    ledger.check(
        "conformal_xi_cancels_induced_channel",
        conformal_inverse == 0,
        f"inverse coupling={conformal_inverse}",
    )
    ledger.check(
        "near_conformal_side_strengthens_G",
        sp.simplify(sp.diff(gravity.subs(baseline, 0), xi)).subs(
            {xi: 0, delta_zero: delta_value}
        )
        > 0,
    )

    payload = {
        "accepted_G": float(G_ACCEPTED),
        "profile_critical_G": float(G_CRITICAL),
        "induced_inverse_coupling_delta0": float(delta_value),
        "required_baseline": float(baseline_boundary),
        "baseline_to_induced_ratio": float(baseline_boundary / delta_value),
        "required_xi_boundary": float(xi_boundary),
        "xi_one_sixth_inverse_coupling": str(conformal_inverse),
        "classification": {
            "baseline": ("regularity is reachable only on a bare-dominated branch"),
            "xi": (
                "weaker positive G requires xi approximately -148.8; moving "
                "toward +1/6 strengthens G before the induced channel cancels "
                "exactly at 1/6"
            ),
        },
    }
    (HERE / "coupling-escapes.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    ledger.finish()


def main() -> int:
    capture = io.StringIO()
    failed = False
    with contextlib.redirect_stdout(capture), contextlib.redirect_stderr(capture):
        try:
            execute()
        except Exception:
            failed = True
            traceback.print_exc()
    output = capture.getvalue()
    (HERE / "stdout.txt").write_text(output)
    print(output, end="")
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
