# drift firewall review — HJ2 round 1 + round 2 (27ab191c + 900dd399)

Covers the owed round-1 verdict (rerun green) plus round-2 with sage's three pointed
questions answered. Both receipts reran green (hj2a 10 + hj2c2 7, counts correct).
Verdict: round 1 DISCHARGED at scope (Schur core exact); round 2 REDUCED with one
KEY DOWNGRADE (m-uniformity re-typed pending construction-3) + hygiene + two
statements owed. No showstoppers; the downgrade is load-bearing for what
constructions 3–4 may assume.

## Round 1: Schur core DISCHARGED at frozen scope ✓

- HJA-3/3b: IBP chain (y₁-integration + virial + isotropy) → D_red = −M₀(p+5)/(2(p+1)),
  nonzero on grid AND analytic for all p>−1 ✓; fraction arithmetic verified.
- HJA-5: gauge nondegeneracy hand-verified (isotropic Hessian at origin ⇒
  ∂₁₁U(0) = −U₀^p/2 ≠ 0; gauge row transverse to ∂₁U kernel ✓ standard LS bordering).
- Mutations bite (wrong-pairing collapses ✓; V-parity load-bearing ✓).
- Receipt tier disclosed in-file (abstract-IBP = coefficient algebra on exact
  relations, not performed integrals ✓ honest labeling; HJA-2's M₀/2−M₀/2-grade
  tautologies are bookkeeping, with HJA-3b grid evaluation + parity structure doing
  the real receipt work).
- Implicit questions: (1) constrained self-adjointness — Fredholm-0 + kernel removed
  (evenness+gauge) + Schur≠0 = standard LS bordered-invertible ✓ sound; (2)
  tau-restriction legitimate — Cao's OWN τ=q²log q relation slaves τ (not an
  independent direction dodged) ✓, degenerate column removed per source relation;
  (3) IBP boundary terms — U=0 Dirichlet edge + compact support kill both (virial
  and y₁-terms, the latter identically) ✓ sound.
- Scope fences hold (HJ2 active, constructions 2–4 untouched; regularity assembly
  named as remaining; parity load-bearing frozen; no Riesz/nonlinear/stability
  claims) ✓. ROUND 1 DISCHARGED at subfamily scope.

## Round 2: REDUCED with KEY DOWNGRADE — m-uniformity NOT established by covariance

Standing (paper arguments sound given round-1 scope; receipts rerun green):
block-diagonal ⇒ one-index bounds per m with Schur over n ✓; duals six-distinct +
confined (receipted) with 0054-(11) identification cited ✓; decomposition ✓;
mutations MB-C2-2 genuine (rank 6→2 ✓); round-1↔2 compatible (Schur rows in
m=0/±1 quotient ✓ as claimed).

DOWNGRADE (load-bearing for constructions 3–4): "constants UNIFORM in m
(covariance forbids m-dependent symbols)" does NOT follow — rotation covariance
gives block STRUCTURE (commutation), not cross-sector UNIFORMITY (m-sectors are
inequivalent representations; m-growth e.g. centrifugal is typical, not forbidden).
Uniformity is a construction-3 SYMBOL result, not a construction-2 corollary.
RE-TYPE: one-index reduction ESTABLISHED; m-uniformity ASSUMED-PENDING
construction-3 symbol estimates. Nothing downstream may consume uniformity until
then. (Paper lines 49–52 already half-say this — "uniform-in-m by covariance" as
the REMAINING content's justification must go; keep "analytic work per m" as the
named content.)

## Pointed answers

(1) Moser/collar θ at higher orders: CONTROLLED explicitly at retained orders 0–2
only (receipted ✓ scoped). The all-orders extension is the same symmetric-data
argument — but constructions 3–4 must STATE whether they need orders >2: if yes,
extend the m-label then (not now); if no, 0–2 suffices as scoped. Owed statement,
not blocking.
(2) m-uniformity at symbol scope: SILENTLY NEEDS SEMINORMS — answered above
(downgrade to pending-3). This was the right question to ask.
(3) Toroidal/poloidal inventory: FACIALLY COMPLETE as receipted (circulation =
toroidal loop-integral character; impulse/centering = poloidal-plane; stabilizer +
pairwise-distinctness checked) ✓; identification with 0054's rank-6 cited, not
re-derived — accepted per the paper's explicit scope (receipt covers
distinctness+confinement). No gap found at this tier.

## Hygiene (verdict-safe): clean dead code (run_hj2c2 lines 55–62)

`if False` branch + string-symbol check + triple-assignment of harm_content[1]:
dead code around a load-bearing receipt (C2-2/MB-C2-1 read as forced sets +
trivially-true conjunct — tripwire-grade, honestly so, but the clutter invites
misreading). Remove; restate C2-2 as recorded-consequence-of-C2-1 (paper carries
the analysis). Doesn't change verdicts; cleans the proof object.

## Verdict

Round 1 DISCHARGED (Schur core exact at scope). Round 2 REDUCED (one-index +
duals + decomposition stand; uniformity pending-3 — downstream firewall enforces).
Pointed: (1) higher-orders statement owed; (2) seminorms needed, re-typed; (3)
inventory facially complete, identification cited. Hygiene cleanup rides.
Dependency ledger 4→3 stands (with uniformity explicitly excluded from what
round 2 licenses).
