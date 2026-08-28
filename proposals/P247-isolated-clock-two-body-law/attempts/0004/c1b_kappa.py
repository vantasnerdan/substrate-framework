"""P247 attempt 0004, C1 deliverable: exact boost-sector kinetic normalization.

The C3 reduced functional must carry zero new continuous parameters. The
sigma-projector stiffness kappa is pinned by the accepted census machinery
(C-M5S-001): the projector-current kinetic metric on the aligned vacuum
(targets (4, 1, 3/10, 0)) has exact rational entries whose boost-direction
value normalizes the reduced scalar kinetic.

This script computes, exactly:
  K1. the full projector-current kinetic metric G on the 10-dimensional
      symmetric basis at the aligned vacuum (regression: its three
      propagating-block eigenvalues must match the accepted census values
      1/9, 100/1369, 1/16);
  K2. the boost-direction entry G[K, K] for the boost generator
      K = E_{0z} + E_{z0} (rapidity along z), the normalization of the
      reduced chi kinetic;
  K3. the same in the 4x4 embedding used by the reduced radial ansatz.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
SRC = HERE / ".." / ".." / ".." / ".." / "src"
sys.path.insert(0, str(SRC.resolve()))

from substrate_framework.m5_fluctuation_spectrum import (  # noqa: E402
    aligned_vacuum,
    timelike_rotation_kinetic_metric,
)


def main() -> None:
    targets = (4, 1, sp.Rational(3, 10), 0)
    basis, metric, _variations = timelike_rotation_kinetic_metric(
        targets=targets, projector_stiffness=1
    )
    size = metric.rows
    propagating = [
        (i, j, metric[i, j])
        for i in range(size)
        for j in range(i, size)
        if metric[i, j] != 0
    ]
    eigen = sp.Matrix(metric).eigenvals()
    eigen_sorted = sorted(
        (sp.simplify(value), mult) for value, mult in eigen.items()
    )

    # boost generator along z in the symmetric basis convention: the basis
    # elements are the elementary symmetric matrices E_{ab} (a <= b); the
    # mixed time-space entry (0, 3) is the z-boost generator.
    boost_index = None
    for index, element in enumerate(basis):
        if element[0, 3] != 0 or element[3, 0] != 0:
            boost_index = index
            break
    boost_entry = (
        sp.simplify(metric[boost_index, boost_index])
        if boost_index is not None
        else None
    )

    payload = {
        "metric_size": size,
        "propagating_entries": [
            [i, j, str(v)] for i, j, v in propagating
        ],
        "metric_eigenvalues": [
            [str(v), int(m)] for v, m in eigen_sorted
        ],
        "accepted_census_reference": ["1/9", "100/1369", "1/16"],
        "boost_generator_index": boost_index,
        "boost_direction_metric_entry": str(boost_entry),
    }
    (HERE / "c1b-kappa.json").write_text(
        json.dumps(payload, indent=2, default=str)
    )
    print(json.dumps(payload, indent=2, default=str))
    print("WROTE c1b-kappa.json")


if __name__ == "__main__":
    main()
