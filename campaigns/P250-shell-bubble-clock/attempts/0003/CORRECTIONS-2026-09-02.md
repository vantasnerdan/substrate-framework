# CORRECTIONS — 2026-09-02 (appended post-review; original files untouched)

Individual review of claims C-M5W-007/008 found three defects in the
attempt-0003 record. Per append-only attempt discipline this note records
the corrections; the original artifacts above are unchanged.

1. Charge normalization: `bag_bvp.py` (lines 338, 345) forms
   `Q_trap = 4 pi w2v * Integral(inertia)` with `w2v = omega^2`, so the
   archived `Q_trap` column is omega^2 * I — one factor of omega too many
   against the accepted C-M5C-001/C-M5W-005 charge Q = omega I. The true
   charge of every stored rung is exactly `Q_true = Q_trap / sqrt(w2)`;
   corrected columns, the envelope identity recomputed at the true charge
   (max rel 1.84e-4, improved from the archived 2.15e-4), and the strict
   decrease of Q(omega) are regenerated in
   `../0005/bag_results_corrected.json` by `../0005/stability_check.py` and
   `../0005/evidence_repair.py` (captured stdout alongside).

2. Acceptance-conjunct retraction: the production acceptance branch of
   `bag_bvp.py` checks solver status, mesh stability of the radius
   (|dR_mesh| <= 1.5e-3 R), and residuals — NOT per-rung virial identities
   or channel monotonicity. `derivation.md` lines 73-79 and any claim
   wording asserting "exact per-rung virial identity and channel
   monotonicity" as acceptance criteria are retracted. The stored normalized
   virial residuals |virial|/virial_norm range 4.19e-4 (delta = 0.001) to
   7.78e-3 (delta = 0.007) and the `mono` field is a spline residual of
   d(T-V)/dr = -4T/r; both are recorded diagnostics, not gates.

3. Selection-law prose: `derivation.md` lines 70-71 describe chi - 1 as
   linear in delta (~ -1.75 delta). The stored rungs contradict linearity:
   (chi-1)/delta = -0.438, -0.369, -0.484, -0.628, -0.922, -1.290, -1.771
   across delta = 0.001..0.007 (log-log slope 1.718). chi approaches 1
   monotonically as delta -> 0+ consistent with a near-quadratic power law.

Downstream impact: chi values shift by <= 1e-11 relative under the G1
sigma_0 repair (see ../0005/sigma0_results_repaired.json); E_crit and the
dQ/d omega sign are unaffected at recorded precision.
