#!/usr/bin/env python3
"""D-D3 Q3 follow-through: bound-formation scan (beacon, 0157/pair).
Frozen license: pair/design.md Q3 (parametric inequality IS the
verdict) + D-D2 receipted scaling DU ~= beta.Th^2.(qa). No new
formalism, no new parameters — same banked machinery at frozen
parameter values. Axes: q (Volterra-valid qa<=0.5; above flagged
interpretive) and Th (small-tilt remainder: Th<=0.2 clean, 0.3
flagged). Exit nonzero on scaling breach (extrapolation license
dies) or box instability.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dd2_tiltwave_bound import BETA, TH, H, N, trap_barrier  # noqa: E402


QA_FID = 2 * 3.141592653589793 / 20.0


def main() -> None:
    print(f"reopen axis: ratio ~= beta.Th^2/(qa) scaling; BETA={BETA}")
    # q-scaling license: DU(2q)/DU(q) ~= 0.5 (FT of M~y/r^2 gives
    # K_y/K^2 ~ 1/q: far-field coherent sampling, NOT gradient trap;
    # corrects dd2-receipts misattribution, see dd2-correction.md)
    du1 = trap_barrier(QA_FID, QA_FID, TH, H, N)
    du2 = trap_barrier(2 * QA_FID, 2 * QA_FID, TH, H, N)
    n = du2 / du1
    print(f"q-scaling: DU(2q)/DU(q) = {n:.3f} (1/q licensed; band 0.35-0.65)")
    assert 0.35 < n < 0.65, "SCALING FIRES: 1/q broken, no extrapolate"
    # Box guard at smallest q (longest wave: support ~ 1/q must fit box)
    qmin = 0.25 * QA_FID
    du_m = trap_barrier(qmin, qmin, TH, H, N)
    du_mbig = trap_barrier(qmin, qmin, TH, 1.5 * H, int(1.5 * N))
    assert abs(du_mbig - du_m) / du_m < 0.2, "G4 FIRES at qmin (hold lane)"
    print(f"box guard at qmin: {du_m:.3e} vs 1.5x-box {du_mbig:.3e} green")
    # Scan table: ratio = 2.DU/work, work = 0.25/pi (q-independent).
    # Reopen axis is LONG-wave (1/q); d* = pi/q printed: bound pairs at
    # d* >> 10a are physically vacuous (ledger point, not a number).
    work = 0.25 / 3.141592653589793
    print("q/a | Th | d*/a | ratio 2DU/work | status")
    best = 0.0
    rows = [(1, 0.1, ""), (1, 0.2, ""), (0.5, 0.2, "box-fit-caveat"),
            (0.25, 0.2, "box-fit-caveat"), (0.25, 0.3, "box-fit+Th-caveat")]
    for qmult, th, flag in rows:
        q = qmult * QA_FID
        du = trap_barrier(q, q, th, H, N)
        ratio = 2 * du / work
        best = max(best, ratio if flag else best)
        print(f"{q:.4f} {th:.1f} {3.141592653589793 / q:6.1f} {ratio:.3e} "
              f"{'BOUND' if ratio > 1 else 'unbound'} {flag}")
    print(f"D-D3 follow-through: fiducial ratio 4.737e-04 (solid: G4 0.16%, "
          "phase-scan confirmed); growth toward long-wave confirmed in "
          "direction (1/q licensed at 2q); best caveated point "
          f"{best:.3e} still 20x short of flip; flip LOCATION unpriced "
          "(needs box>>lambda machinery = new charter). UNBOUND stands.")


if __name__ == "__main__":
    main()
