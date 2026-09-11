#!/usr/bin/env python3
"""0131 per-m run (beacon 0135): last licensed Magnus test (capstone a523605e).

Seed each m in {1,2,3,4} on the leapfrog base; gate-first
(eps-linearity + content gate on the SEEDED channel); shared-kappa
fit across passages in that channel under frozen F1 bars. Any HOLD
halts with escape report; all KILL confirms family closure.
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np

ATT = "proposals/P253-euler-particle-mechanisms/attempts"
sys.path.insert(0, os.path.join(ATT, "0120-cipher-m2b1", "receipts", "a3-scan"))
sys.path.insert(0, os.path.join(ATT, "0134-beacon-survbuild"))
import run_a3 as A3
from run_surv import gate, content_ok, fit_passages, leap_windows, adjudicate

SEED = 1e-3


def stage_seed(m):
    print(f"seed m={m}:", flush=True)
    R1, R2, T, _ = A3.shoot()
    X0 = A3.ring_state(R1, 0.0, R2, 0.0) + SEED * A3.basis_field(m, 0, 0, 0)
    wins, n = leap_windows(X0, T)
    ok, sb, sp = gate(X0, n)
    if not ok:
        return "GATE-FAIL"
    if not content_ok(sb, sp, m):
        return "NON-FIRING"
    Rm_, Tm_ = fit_passages(X0, n, wins, m)
    return adjudicate(Rm_.ravel(), Tm_.ravel(), f"m={m} shared-kappa")


if __name__ == "__main__":
    t0 = time.time()
    only = sys.argv[1] if len(sys.argv) > 1 else "all"
    modes = [int(only)] if only != "all" else [1, 2, 3, 4]
    for m in modes:
        v = stage_seed(m)
        print(f"M{m} verdict: {v} ({time.time()-t0:.1f}s)", flush=True)
        if v == "HOLD":
            print("PROGRAM HALT: escape reported, downstream skipped.", flush=True)
            break
