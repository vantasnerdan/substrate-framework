#!/usr/bin/env python3
"""T3a defect-induced static tilt NULL (beacon, 0157/tilt, COMPUTE AUTHORIZED).
Frozen: tilt/t3-design.md F-T3a. Banked inputs only: w = 1/2.J(eps).th^2
+ W_FA(eps) (t2-design.md), J(eps) = J0(1+4.eta_in-3.eps_zz) (edge-T1').
Claim: theta-EL source is propto th -> th = 0 exact for ANY defect
strain (stiffness modulation, never a source). Exit nonzero on any
th-independent defect term (F-T3a).
"""

import sympy as sp

J0, eta, ezz, th = sp.symbols("J0 eta ezz th", real=True)


def main() -> None:
    J = J0 * (1 + 4 * eta - 3 * ezz)  # banked J(eps), edge-T1'
    w = sp.Rational(1, 2) * J * th ** 2  # theta-sector density (W_FA indep.)
    EL = sp.diff(w, th)  # static theta equation LHS (no grad stiffness banked)
    print(f"theta EL: {EL} = 0")
    # F-T3a guard: source at th=0 must vanish IDENTICALLY in (eta, ezz)
    src0 = sp.simplify(EL.subs(th, 0))
    assert src0 == 0, f"F-T3a FIRES: th-independent defect source {src0}"
    print("F-T3a guard: no th-independent defect term (exact, arbitrary eps)")
    # Factorized form: stiffness modulation, never source
    assert sp.simplify(EL - J * th) == 0
    print("structure: EL = J(eps).th -> th=0 solves for ANY defect strain")
    # Numeric spots (implementation guards, banked numbers only)
    spots = {"edge-peak": (0.11818, 0.0), "screw": (0.0, 0.0),
             "far-field": (0.0, 0.0)}
    for name, (e, z) in spots.items():
        Jv = float(J.subs({J0: 1.0, eta: e, ezz: z}))
        assert Jv > 0, f"stiffness non-positive at {name}"
        print(f"{name}: J/J0 = {Jv:.5f} (finite positive) -> th = 0 solves")
    # Cross-check (not a new verdict): edge-peak modulation 4*.11818
    # = 0.4727 reproduces edge-T1' |dJ/J0| = 0.47274 to 4dp
    assert abs(4 * 0.11818 - 0.47274) < 1e-4, "T1' consistency breach"
    print("consistency: edge-peak dJ/J0 = 0.4727 matches T1' 0.47274 "
          "(same J, no new claim)")
    print("T3a verdict: static tilt from defect NULL at derived order; "
          "T3b Born contrast next (ordered round).")


if __name__ == "__main__":
    main()
