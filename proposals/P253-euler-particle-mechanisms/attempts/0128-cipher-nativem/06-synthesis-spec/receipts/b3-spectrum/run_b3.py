"""B3 fluctuation spectrum (FROZEN F-bar in ../../06-synthesis-spec/02-b3-fbar.md):
fixed carrier, many tangle seeds; Fano factor + deciles + seed-convergence.
Reuses the S4 pilot module by explicit path (same repo, banked artifact).
Usage: python3 run_b3.py
"""
import os
import sys
import time

sys.path.insert(0, os.path.join("..", "..", "..", "05-s4", "receipts", "tangle-pilot"))
import numpy as np
import run_tangle as T

S = 256
RB, NT = 3.0, 40


def ensemble(S, Rb=RB, Nt=NT):
    C = T.carrier()
    L = np.zeros(S)
    for s in range(S):
        rng = np.random.default_rng(2000 + s)
        L[s] = T.stats(C, T.tangle(rng, Nt, Rb, T.RT))
    return np.abs(L)


def main():
    t0 = time.time()
    A = ensemble(S)
    F = A.var() / A.mean()
    dec = np.round(np.quantile(A, [0, 0.25, 0.5, 0.75, 0.9, 1.0]), 4)
    B = ensemble(128)
    Fb = B.var() / B.mean()
    conv = abs(F - Fb) / F
    print(f"S={S}: Fano={F:.4f} deciles={list(dec)}")
    print(f"S=128: Fano={Fb:.4f} conv={conv:.3f}")
    hold = 0.5 <= F <= 2.0 and conv < 0.20
    kill = F < 1 / 3
    print(f"B3 verdict: {'HOLD' if hold else ('KILL' if kill else 'UNRESOLVED-gray')} ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
