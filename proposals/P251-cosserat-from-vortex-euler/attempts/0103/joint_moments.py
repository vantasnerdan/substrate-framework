"""Exact finite-field extension of0085's frozen compact differential syzygy."""

from pathlib import Path
import runpy

import numpy as np

from substrate_framework.verification import CheckLedger


def main():
    # This immutable prototype builder also prints its historical probe.
    # Only the new augmented ranks below answer the present moment question.
    source = Path(__file__).resolve().parents[1] / "0085" / "jet_probe.py"
    data = runpy.run_path(str(source))
    matrix, angular = data["matrix"], data["angular"]
    indices, prime, rank = data["indices"], data["PRIME"], data["modular_rank"]
    generator = np.zeros_like(angular)
    for axis in range(3):
        unit = [int(i == axis) for i in range(3)]
        for column, alpha in enumerate(indices):
            if sum(alpha) != 1:
                continue
            # G=integral (e cross r).xi. Right-normal adjunction of one
            # derivative gives minus d_alpha(e cross r), independently of u.
            derivative = data["cross"](unit, list(alpha))
            for component in range(3):
                generator[axis, component * len(indices) + column] = -derivative[component]
    generator %= prime
    rank_a = rank(matrix)
    rank_l = rank(np.vstack([matrix, angular]))
    rank_g = rank(np.vstack([matrix, generator]))
    rank_joint = rank(np.vstack([matrix, angular, generator]))
    print("NEW constraint, L, G, joint ranks:", rank_a, rank_l, rank_g, rank_joint)
    ledger = CheckLedger("P251-0103-joint-compact-angular-moments")
    ledger.check("unchanged prototype attains the proved universal constraint bound", rank_a == 235)
    ledger.check("existing induced-velocity spin rows retain their independent rank", rank_l == 238)
    ledger.check("three generator angular rows are independent modulo constraints", rank_g == 238)
    ledger.check("all six actual angular moment rows are jointly independent", rank_joint == 241)
    # Reusing L in place of G would erase the newly required freedom.
    ledger.check("conflating generator and velocity rows loses three directions",
                 rank(np.vstack([matrix, angular, angular])) == 238)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
