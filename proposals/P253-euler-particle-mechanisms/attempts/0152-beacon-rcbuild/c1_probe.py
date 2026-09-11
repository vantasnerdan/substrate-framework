#!/usr/bin/env python3
"""C1 static-stiffness probe (beacon 0152, banked per drift repair).

k_dir = -<dV,dx>/<dx,dx> per m-channel on the A3 leapfrog orbit
(centered FD, e = 1e-6). Result: -0.0000e+00 all channels
(gyroscopic zero; |dV| live). A3 code read-only.
"""

from __future__ import annotations

import os
import sys

import numpy as np

ATT = "proposals/P253-euler-particle-mechanisms/attempts"
sys.path.insert(0, os.path.join(ATT, "0120-cipher-m2b1", "receipts", "a3-scan"))
import run_a3 as A3


def main() -> None:
    R1, R2, T, _ = A3.shoot()
    X0 = A3.ring_state(R1, 0.0, R2, 0.0)
    e = 1e-6
    for m in (1, 2, 3):
        for k in (0, 2):
            F = A3.basis_field(m, 0, k, 0)
            nF = np.linalg.norm(F)
            if nF == 0:
                continue
            Fn = F / nF
            dV = (A3.rhs3(X0 + e * Fn) - A3.rhs3(X0 - e * Fn)) / (2 * e)
            P = A3.project(dV)
            tot = np.sqrt(sum(float(np.linalg.norm(v) ** 2) for v in P.values()))
            kdir = -float((dV * Fn).sum())
            print(f"m={m} k={k}: k_dir={kdir:+.4e} |dV|_tot={tot:.3e}",
                  flush=True)


if __name__ == "__main__":
    main()
