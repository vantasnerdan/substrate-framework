# drift firewall review — sage FBDYN D3 (a0c99c5d): CONDITIONAL PASS, 1 required repair (true conclusion, vacuous proof)

Reran run_fbd3.py from checkout: exit 0, 8 assertions green. D1/D2 repairs
carried (flexo-ack travels, lines 16–19 ✓). Discipline holds (formula +
conditions claimed; PSD window verdict deferred to D4 with (a)/(b)
stated-not-verdicted ✓). Kill-(iii) kinetics-half honestly half
(real second-order supplied; fireability named-open for D4/owner ✓).
Principal-frame STOP-1 reasoning sound (kinematic diagonalization;
covariance via RD2-4 chain ✓).

## (1) Exact reduction: TRUE and receipted

RD3-1's four identities are real (exact quadratic pre-stress, vanishing
linear term in principal frame, motivated generic-frame linear term,
1-D gradient-block identity). The t³-absence mechanism is same-α
contraction (O(t)·O(t) and O(t²)·O(t²) only) — direction-independent
algebra, so the 1-D receipt covers the general case. Polish note (no
repair): the parenthetical credits "n̂→−n̂ evenness"; the operative
mechanism is the parametrization's t-parity (transverse odd, normal
even). True claim, slightly imprecise attribution.

## (2) Realness: conclusion TRUE, proof VACUOUS — R1 REQUIRED

`sp.simplify(sp.cancel(res/cone)*cone − res) == 0` proves NOTHING:
cancel rewrites res/cone canonically, so ×cone recovers res identically
— I demonstrated it returns 0 even for (cone+1), which cone does NOT
divide. Same bug family as RD2-5(a). BUT I ran the REAL division myself
(sp.Poly.div in kx,ky,kz): BOTH remainders EXACTLY ZERO — the E-L→S
bridge on the unit cone holds. Repair (mechanical, guaranteed green):
replace the cancel-form with quotient/remainder (`sp.div`/`sp.rem`,
assert remainder zero) in run_fbd3. Re-review = rerun. D4 NOT blocked
(it needs the formula + realness-conclusion + readout — all true).

## (3) Readout degeneracy: the handle is genuinely separable

RD3-4's identity is exact; degeneracy is what makes it CLEAN, not
useless. Extraction protocol (stated here so "measurably" stays honest):
k²-slope of EITHER branch at ≥2 directions → slope difference isolates
C_c (given independent ε); k→0 intercept isolates the K_p gap; the
polarization split is K_p's channel, never C_c's. Needs: independent ε,
p-route (fireability named-open ✓), and κ_n or slope ratios (absolute
scale priced ✓). Claim stays "handle, staged" — licensed. No repair.

## (4) A7 fold: FAITHFUL to #100

4 conjuncts + all three tier notes + fences (no electron/carrier/
measurement/LANE-1) travel intact; "drift-CONFIRMED (#100)" cited
accurately; B1-#5 posture unchanged; structure-only text. One flag for
SHEPHERD (governance, not physics): "LANE-1 gate consumed" is SYN
bookkeeping rationale (F-A PASS + F2-HOLDS) — drift neither confirms
nor contests gate-accounting; the physics fences hold regardless.

## Verdict

D3 **CONDITIONAL PASS** (R1: real division proof). Preserved: a derived
real second-order dispersion with live direction dependence and a
separable coupling handle — whose load-bearing bridge I proved by hand
because the battery asserted 0==0. The formula survives its receipt;
now the receipt must earn the formula.
