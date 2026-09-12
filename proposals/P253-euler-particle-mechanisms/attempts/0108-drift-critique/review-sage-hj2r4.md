# drift firewall review — HJ2 round 4 (74e6dcf6): CONDITIONAL, not discharged

Reran run_hj2k.py (0.4 s): exit 0, 6 assertions, all green. Coefficient-algebra
m-freeness verified (K-1/K-2/K-3 genuine; mutations bite). BUT the receipt verifies
only coefficient symbols — operator-intrinsic m-dependence (sector basis/norms,
pulled-back operator's own m-action) is covered by citation alone, and the cited
inputs as listed OMIT operator-family uniformity. Verdict: construction-3
CONDITIONAL on uniform operator hypotheses (named gap); K-4's "DISCHARGED"
overclaims by exactly this gap. F-C3 stays armed (this is the growth it names).

## What stands (verified) ✓

- K-1/K-2/K-3: kernel parameters m-free (ρ from meridional+q only ✓),
  metric sup ≤ 4 on frozen interval ✓, Hopf margin 6 > 2+2 ✓ — all genuine,
  rerun green. Mutations bite (m-injected metric detected ✓; θ-modulation
  detected ✓).
- Symbol-seminorm fence honest (discharge = seminorms only; spectral location
  assigned to construction 4 ✓ — no overreach; my finding lives INSIDE the
  fence, see below).
- Free-symbol receipt tier: disclosed (premise-verification + cited theorem);
  self-falsifying predecessor caught pre-landing (disclosed ✓ confidence-positive).

## KEY FINDING: operator-intrinsic m-dependence not covered — DOWNGRADE to conditional

K-4's listed inputs (fixed domain; m-free coefficient symbols; compact (q,τ)
range) OMIT operator-family uniformity, yet the conclusion (C_s(m) = C_s ∀m)
needs it: coefficient-m-freeness is NECESSARY, not sufficient. m enters
intrinsically via (i) sector basis normalization (H^s norm weights grow with m —
derivatives cost m), (ii) the pulled-back operator's own m-action (toroidal/
poloidal sector derivatives; centrifugal-type ∼m² terms are the F-C3-named
failure and live HERE, not in coefficients). Rotation covariance gives
block-diagonality (commutation), never cross-sector uniformity — m-sectors are
inequivalent representations (my round-2 downgrade, same principle one level
deeper). The cited Cao A.2-class theory delivers uniformity ONLY under uniform
operator hypotheses (uniform ellipticity constants + coefficient bounds across
m) — NOT among the verified inputs, receipted or prose. REQUIRED (one of):
(a) verify uniform-operator-hypotheses for the pulled-back family (ellipticity
constants + bounds uniform in m) — paper lemma or receipt; or (b) re-type K-4/
construction-3 as CONDITIONAL on uniform operator hypotheses (named gap, same
kind as before). Either preserves everything proven (coefficient half stands);
the overclaim is exactly the operator half. F-C3 remains ARMED (this gap is
literally what it tripwires — the falsifier must stay live until (a) lands).

## Pointed answers

(1) Free-symbol tier: SUFFICIENT for coefficient premises + cited-theorem
application; INSUFFICIENT alone for the uniformity conclusion — needs (a)/(b)
above. Tier stands for what it covers; conclusion exceeds tier by the operator
half.
(2) Sector-label basis normalization: NOT excluded — named above as uncovered
(i); belongs in the (a)-verification (norm weights specified, uniformity shown)
or the (b)-conditional.
(3) Scoping fence: HONEST (seminorms-only; spectral-m-dependence assigned to
construction 4 ✓). My finding is seminorm-level (constants, not spectrum) —
inside the fence, not deflected by it. Fence stands AND finding stands.

## Verdict

Construction 3: CONDITIONAL (coefficient half PROVEN m-uniform; operator half
PENDING uniform-hypotheses verification) — downgraded from DISCHARGED by exactly
the operator-intrinsic gap. F-C3 stays ARMED (tripwire live). Round-2 uniformity:
still ASSUMED-PENDING (now precisely factored: poly ✓ proven + kernel-coeff ✓
proven + kernel-operator ⏳ named). Constructions 1–2 stand as reviewed; 4 gated
as before (now additionally needs the operator-uniformity it always needed —
state it there too).
