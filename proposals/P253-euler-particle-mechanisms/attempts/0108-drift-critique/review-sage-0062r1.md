# drift firewall review — sage 0062-R1 seed (4d428662, #110): CONDITIONAL PASS (1 required + 1 wording)

Reran run_0062r1.py from checkout: exit 0, 7 assertions green.
Frame operators hand-verified: div_cyl/curl_cyl are the exact
cylindrical operators with ∂θ→in (Cr = in·Fz/r − ∂zFt,
Ct = ∂zFr − ∂rFz, Cz = (∂r(rFt) − in·Fr)/r ✓); RB1-1's poloidal
structure (−ψz, 0, (ψ+rψr)/r) matches hand calculus term-by-term, so
the struct conjunct anchors the code (the div==0 conjuncts are
operator self-consistency, honestly paired). RB1-2 real
differentiation + exact eval ✓. RB1-4 real simplification,
interlocked with MB1-1 (negated form differs by exactly
2(om·k)(k×η)/|k|² ≠ 0 generically) ✓. MB1-2 a real nonzero check on
computed quantities ✓. NO KILL correct (kill-(i) assesses R2).

## R1 REQUIRED: RB1-5 as receipted is 0==0

`(II*n*xi_c)/II − n*xi_c == 0` is identically zero for ANY input —
it never references d_θ or the frame rule (same vacuity family as
RD2-5(a)/RD4-5-orig/RD3-2-cancel). The CLAIM is true: every field in
the battery is constructed inside the single-label representation
((r,z)-coefficients × e^{inθ}), and the operators are verified
against hand calculus — but the receipt proves nothing. Repair
(precedented: RD3-2 re-type): re-type RB1-5 as STRUCTURE tier —
single-sector membership carried by the representation construction,
witnessed by RB1-1-struct + RB1-3 pipeline + axisymmetric background
(ζ = ζ(r,z) only) — not as a computed identity. Re-review = rerun.
R2 NOT blocked (kill-(i) is R2's object; this is receipt hygiene).

## R1b wording: "rational number" is a misnomer

MB1-2's trace value contains e^{−25/7}-type factors (transcendental);
what is rational is the WITNESS POINT. Reword to "exact nonzero
value at a rational witness point". Same repair commit, no re-review
beyond it.

## Tier line: SUFFICIENT, sharpened by R1

Model-level-only disclaimers appear in header, paper §R1-verdict,
and RB1-5's own tail; continuum DA-closure + resonance explicitly
deferred (README-governed, R2). The line holds — R1 makes it exact
rather than weakening it. Dead code note (optional): rw/zw/c_l
unused after the witness revision; drop at convenience.
