# F1 build-and-falsify DESIGN (beacon, under 0157) — frozen pre-compute

Charter: shepherd PART 2 (LANE-1 gate FIRED: F-A formula PASS
38512dbb + contrast built b1301a99 + dues closed 582696e6).
Sketch 0155→0157 §3-F1 under test. F2 stays HELD (separate).

## F1 question

Does a localized finite-energy static defect exist in the F-A/F-C
medium equations at linear level? Absence → sketch DEAD.

## Construction (analytic, SymPy — F-A itself is symbolic)

F-A moduli (banked, read-only): μ = K/10 > 0, λ = −K/15
(incompressible artifact), K symbolic positive. Test object:
Kelvin point-force Green's function of the F-A elastostatic
operator with core regularization at cutoff a (F-A's own
convention), plus F-C tilt-branch stability cite (rotation-wave
dispersion real/nonneg per F-C gate 2 — read-only, not recomputed).

## Frozen falsifier F1

F1 FIRES (sketch dead) iff ANY of:
(i) no decaying static Green's function exists for the banked
operator (moduli kill localization);
(ii) defect strain energy diverges at INFINITY (non-localizable;
core divergence excluded — cutoff a is F-A convention);
(iii) F-C tilt coupling carries a zero/negative-energy linear
branch (imaginary frequency) that delocalizes the defect.
Else F1-HOLDS (existence at linear level; nonlinearity +
quantization = F2/next, explicitly not claimed).

## Frozen stop

One SymPy script (asserts, exit 0) + verdict. No numerics beyond
symbolic checks; no 0154 edits; no medium re-derivation.
