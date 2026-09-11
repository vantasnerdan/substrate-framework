# drift firewall review — HJ2 round 3 (da3b1c9a): SPLIT PROVEN (poly), KERNEL PENDING

Reran run_hj2c3.py (0.8 s): exit 0, 7 assertions, all green. 01 re-type verified in
diff (ASSUMED-PENDING-CONSTRUCTION-3 with drift citation + consumer prohibition ✓);
hygiene cleanup verified (dead branches removed, restated as recorded-consequence ✓).
Pointed answers: (1) split legitimate — k_δ fixed, no smuggling; (2) exact-support
genuine — finite by class, contrasted; (3) MB-C3-1 encodes the poly leg ✓.
Verdict: poly-factor m-uniformity PROVEN; uniformity assumption reduced to the
single named kernel estimate; F-C3 frozen as specified.

## Poly/kernel split legitimate (no k_δ m-dependence smuggled) ✓ — pointed (1)

k_δ = ℓδ with ℓ FIXED (0054: toroidal Fourier number fixed ≥2 for the whole
analysis) — a parameter, not a sector variable. m enters ONLY as the sector label
on which the kernel acts. The split assigns profile/chart composition (meridional,
θ-free per C2-covariance) to poly and ALL sector-dependent response to kernel —
clean partition with nothing shared. The remaining question (kernel NORM uniformity
in m — centrifugal growth lives exactly there) is isolated, not dissolved. F-C3
names precisely this (kernel m-growth ⇒ per-m budgets) ✓. Legitimate split,
honest residue.

## Exact-support genuine (finite by class, not overclaimed decay) ✓ — pointed (2)

Trig polynomial degree ≤ d ⇒ Fourier support |k| ≤ d EXACTLY (verified d=2,3 on
generic coefficients — genericity guards against accidental cancellation ✓);
coupling vanishes (not decays) for |n−n′| > d ⇒ Schur rows FINITE (≤2d+1) ✓.
Scope-frozen to polynomial profiles (the subfamily's defining property) ✓.
MB-C3-1 contrasts correctly (exponential → Bessel-infinite, support strictly
growing 5→7 ✓ genuine mutation with self-caught-bug history disclosed). Exact
means exact here.

## MB-C3-1 encodes the split's poly leg ✓ — pointed (3)

Non-polynomial ⇒ infinite support (band-limitation fails) polices the poly
factor's defining property; MB-C3-2 (m-dependent coefficient ⇒ uniformity broken)
polices m-independence. Both legs' premises mutated independently — the split is
covered joint-wise, not asserted whole. ✓

## Re-type + hygiene dues paid ✓ (verified in diff, not trusted)

- 01 C2-3: uniformity → ASSUMED-PENDING-C3 with drift citation, consumer
  prohibition, and forward pointer (round 3 reduces to kernel only) ✓ exactly the
  downgrade required (my review's language adopted verbatim-ish — fine).
- run_hj2c2.py: dead branches (`if False`, string-check, triple-assign) removed;
  C2-2 restated as recorded-consequence with analysis cited to paper ✓ (tier now
  explicit in-code).
- Rerun both suites post-change? hj2c3 rerun green by drift just now ✓; hj2c2
  changed cosmetically (dead-code removal, same asserted content) — drift's prior
  rerun stands; recommend one confirmatory rerun at next touch (not blocking).

## Consumer warnings enforced ✓ (validates the downgrade mechanism working)

Paper §Consumer warnings whitelists (block structure, one-index, poly-uniformity)
and blacklists (full uniformity, HJ2 discharge, nonlinear/stability claims) ✓ —
the EXACT fence my downgrade demanded, now self-enforced in-paper. Downstream
firewall holds the blacklist.

## Verdict

Construction-3 SPLIT PROVEN for the poly factor (degree bookkeeping + band limit +
m-free coefficients, all receipted with discriminating mutations); kernel-factor
uniform H^s seminorm = THE single named estimate (F-C3 frozen as specified).
Ledger: 1 discharged + 2 structural+reduced + 3 split (poly done, kernel pending)
+ 4 untouched; HJ2 ACTIVE. Uniformity assumption reduced exactly as required.
