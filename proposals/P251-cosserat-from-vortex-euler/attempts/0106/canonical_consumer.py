"""Thin modular regression on the immutable 0085 fixture, no rational rerank."""

import contextlib
import io
from pathlib import Path
import runpy

import numpy as np

from substrate_framework.euler_compact import compact_isovortical_jet_system
from substrate_framework.verification import CheckLedger


def main():
    source = Path(__file__).resolve().parents[1] / "0085" / "jet_probe.py"
    captured = io.StringIO()
    with contextlib.redirect_stdout(captured):
        fixture = runpy.run_path(str(source))
    print("Immutable modular fixture replay (historical scope warning unchanged):")
    print(captured.getvalue(), end="")
    # These integer representatives are used only modulo PRIME. They are
    # not declared to be the characteristic-zero jets of a new background.
    system = compact_isovortical_jet_system(fixture["jets"], fixture["ORDER"])
    prime = fixture["PRIME"]

    def mod(matrix):
        return np.array(matrix.applyfunc(lambda value: int(value) % prime).tolist(),
                        dtype=np.int64)

    constraints = mod(system.constraints)
    velocity = mod(system.velocity_angular)
    generator = mod(system.generator_angular)
    ledger = CheckLedger("P251-0106-importable-compact-jet")
    ledger.check("operator and output column conventions reproduce the frozen fixture",
                 list(system.operator_indices) == fixture["indices"]
                 and list(system.output_indices) == fixture["outputs"])
    ledger.check("every constraint coefficient reproduces the independently built fixture",
                 np.array_equal(constraints, fixture["matrix"]))
    ledger.check("every induced-velocity angular coefficient reproduces the frozen fixture",
                 np.array_equal(velocity, fixture["angular"]))
    rank_joint = fixture["modular_rank"](np.vstack([constraints, velocity, generator]))
    print(f"canonical six-row fixture rank modulo {prime}: {rank_joint}")
    ledger.check("the same selected fixture retains all six independent angular rows",
                 rank_joint == 241)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
