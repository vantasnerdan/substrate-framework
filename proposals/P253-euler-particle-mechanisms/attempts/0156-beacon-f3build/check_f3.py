#!/usr/bin/env python3
"""F3 feed-consistency check (beacon 0156): frozen design. Banked
numbers only; asserts exit nonzero on clash-confusion, prints verdict.
"""

# Banked inputs (see design.md for provenance)
RESP_LO, RESP_HI = 0.3, 1.0      # 0151 skirt-annulus response
KINK = 0.2                        # 0151 kink band (~0.2x)
A, ELL = 0.15, 1.0                # 0117 core / ring-scale mapping
SLOP = 3.0                        # frozen mapping-weakness price


def main() -> None:
    lo, hi = A / SLOP, ELL * SLOP
    print(f"window [{lo:.3f}, {hi:.3f}] vs response [{RESP_LO}, {RESP_HI}]")
    clash = (RESP_HI < lo) or (RESP_LO > hi)
    print(f"kink {KINK} vs a/3 = {lo:.3f}: "
    print("F3 verdict:", "F3-FIRES (sketch DEAD)" if clash
          else "F3-HOLDS (sketch survives F3; F1/F2 still gated)")


if __name__ == "__main__":
    main()
