# P250 attempt 0005 — thin-wall bag stability split (analytic; no numerics)

Scope decision (owner, 2026-09-02): the stability statement rests on calculus
over the accepted reduced family; the full Hessian-spectrum numerics is
descoped. This derivation supplies the exact second-variation sign in the
reduced radial mode and the charge-monotonicity criterion across the
certified G2 family.

Primary evidence: the executables in this directory —
`stability_check.py` (18 file-based checks; captured stdout in
`run_stdout.txt`) and `evidence_repair.py` (review-driven evidence repair;
captured stdout in `run_repair_stdout.txt`). This prose is the reading
guide, not the record.

## Exact curvature sign at fixed omega

Within the C-M5W-005 reduced thin-wall family the fixed-omega bag energy as a
function of the wall radius is

    F_omega(R) = 4 pi sigma(omega) R^2 - (4 pi/3) p(omega) R^3 ,

with sigma(omega) the wall tension (C-M5W-001 functional; value layer
C-M5W-006) and p(omega) = -min V_omega > 0 the bulk pressure. Then

    F'_omega(R)  = 8 pi sigma R - 4 pi p R^2 = 4 pi R (2 sigma - p R),
    F''_omega(R) = 8 pi sigma - 8 pi p R = 8 pi (sigma - p R).

Stationarity gives the exact C-M5W-005 selection law R_c = 2 sigma/p, and at
that radius

    F''_omega(R_c) = 8 pi (sigma - 2 sigma) = -8 pi sigma < 0   (sigma > 0).

So at fixed omega the stationary bag is a strict maximum of the reduced
radial energy: Morse index >= 1; the bag is a critical nucleation bubble
(saddle), not a minimum. This closes the radial (volume-vs-surface) part of
the linear-stability frontier named in C-M5W-005's exclusions. The sign is
exact algebra (SymPy-verified in stability_check.py); no discretization
enters.

Numeric cross-check at rung 1 (attempts/0003/bag_results.json): sigma_0 =
0.7292984178707203, p at rung 1 (omega^2 = omega*^2 + 0.001) =
0.001197535401340061 => R_c = 2 sigma_0/p = 1218.00 vs R_read = 1217.47
(chi = 0.999562), and F''(R_c) = -8 pi sigma_0 = -18.3293 < 0. The critical
energy F(R_c) = (16 pi/3) sigma_0^3/p^2 = 4.53198e6 vs E_gk = 4.53903e6
(0.16%), matching the recorded rung-1 agreement.

## Charge-constrained branch

Along the stationary family Q(omega) = omega I(omega) (C-M5W-005); the slice
inertia density f^2 + 4b^2 is nonnegative and nonzero on any nontrivial bag
(recorded interiors iota_B = 2.395..2.404), so the integrated I > 0. Across
the certified C-M5W-007 rungs (delta = 0.001..0.007, i.e. omega increasing
away from omega*):

    R: 1217.47 -> 171.23 strictly decreasing  =>  I(R) strictly decreasing
    Q_true: 2.3351e10 -> 6.5152e7 strictly decreasing (corrected columns)

so dQ/d omega < 0 holds across the certified family (numeric ingredient,
recomputed at the true charge normalization Q = omega I in
stability_check.py after the archived Q_trap column was found to carry an
extra factor of omega).

Scope of the constrained-branch statement (review correction): what is
asserted is criterion satisfaction — the family satisfies dQ/d omega < 0 —
plus the exact envelope identity dE_data/dQ = omega (chain rule on the
accepted Legendre scope dE_w/d(omega) = -Q). A constrained-minimum
conclusion for the fixed-charge ensemble is NOT asserted: the reduced
second-variation calculation that would supply it is open frontier, and no
accepted dependency supplies it.

## Permitted verdict and what is NOT asserted

Permitted: fixed-omega Morse index >= 1 in the reduced radial mode (exact);
criterion satisfaction dQ/d omega < 0 across the certified family at the
true charge normalization; envelope identity (exact).

Not asserted (open frontier, descope by owner direction): the full
second-variation operator spectrum (shape modes, angular-momentum l-lift,
thick-wall regime); any constrained-minimum/Hessian conclusion in the
fixed-charge ensemble; non-spherical or dynamical stability; any full-space
lift beyond the slice.
