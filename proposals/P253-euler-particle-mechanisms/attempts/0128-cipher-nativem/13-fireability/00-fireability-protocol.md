# 13 — Fireability protocol for FB-D-waves (structure-only design, no build)

Charter: shepherd, owner-level open. Converts the FB-D-waves falsifier from
UNFIREABLE-named to executable by supplying the missing construction: an
independent-p measurement protocol. Drift reviews. No compute, no receipts,
no verdict change to the FB-D lane — this doc designs the route; a future
chartered build executes it.

Frozen inputs (cited, not re-derived):

- Falsifier + fireability clause: `0160/00-fbdyn-scope.md` §4. Kill-(iii):
  measured director-wave response with ω² NOT ∝ k² at small k, beyond priced
  untracked orders, **with p independently fixed**. Fireability: FIREABLE only
  with a p-route independent of the wave measurement; else UNFIREABLE stands.
- Dispersion: `0160/03-fbdyn-d3.md` —
  ω²_± = [K_n k² + 2C_c k²(k̂ᵀMk̂) + K_p λ_±(T)] / I_n, real second-order,
  k-even, polarization-degenerate coupling readout (RD3-4).
- Stiffness + scaling: `0160/01-fbdyn-d1.md` —
  K_n = (8π/3) K p² M₄ ξ² > 0 on p ∈ (0,1], K_n(0) = 0 identically
  (no director waves in the unpolarized medium; MB-D1-2).
- Window: `0160/04-fbdyn-d4.md` — PSD on {p ∈ (0,1)} × {λ₋(T) ≥ 0,
  k̂ᵀMk̂ ≥ −5}; inertia I_n = ρξ³κ_n priced-not-derived.
- Static moduli (F-B-lite, SYNTHESIS A7-line + 0154/07): μ_⊥ = K(1−p)/10,
  μ_∥ = K(3p+2)/20. Ratio r ≡ μ_∥/μ_⊥ = (3p+2)/(2(1−p)), inverted
  p = 2(r−1)/(2r+3). K-free but NOT model-free (assumes the F-B-lite
  static form + declared dipole statistics).
- Hygiene envelope H (fractional, banked budgets): ξ/λ separation, κ_n
  pricing, multipole/quadrupole footnote, acoustoelastic O(p) budget
  (0154/07 RB10b). MISS inside H = INCONCLUSIVE, never kill.

## 1. Circularity rule (what counts as independent)

Write kill-(iii) with p a free symbol: KILL iff
|O_meas − O_pred(p*)| > H + δO_p, where p* is the inserted value and δO_p
the p-induced band. The route is:

- **CIRCULAR** if p* is a function of (a) the wave data under test, (b) any
  parameter fit to that wave data, or (c) the μ_∥/μ_⊥ ratio (or any static
  modulus derived under the same F-B-lite form whose extension the wave test
  adjudicates — a miss then cannot be localized to the dynamical prediction,
  and p absorbs the discrepancy: fit-then-test on one model).
- **INDEPENDENT** if p*'s calibration chain terminates in standards disjoint
  from the wave channel AND from F-B/FB-D fitted parameters (field strength,
  geometry, counts, an independent physical response with its own budget).
- **Procedural bar (required at build):** p-team seals p* with its error bar
  BEFORE the wave fit opens (blinded insertion). Organizational separation is
  part of independence, not decoration.

The μ-ratio route audit: r → p map is exact K-free algebra, but it rides the
static form; using it as p* makes kill-(iii) a joint static+dynamics test and
lets p soak up dynamical misses. Verdict: CIRCULAR for kill-(iii) purposes
(usable as cross-check only, never as p*).

## 2. Candidate p-routes (observable sequences)

**R0 — p = 0 null test (no p-precision needed).** Preparation: unpoled,
annealed, symmetry-guaranteed p = 0 (no measurement; premise = preparation
symmetry + declared statistics). Prediction: NO director waves, K_n(0) = 0
identically. Any propagating k² response kills (MB-D1-2 fires). Cheapest leg;
run first. Failure mode it cannot touch: it tests existence-at-zero, not the
D3 formula shape — a pass leaves kill-(iii) unfired.

**R1 — Birefringence Δn ∝ p (optical channel).** Sequence: (1) calibrate
Δn vs preparation on the saturation asymptote (p → 1 scale anchor, §3);
(2) read Δn at the wave-test preparation; (3) insert p* ± δp sealed.
Independent physics (optical polarizability vs elastic moduli); shares only
the p definition + dipole density. Tier-1.

**R2 — Dielectric anisotropy Δε/ε̄ ∝ p (electrostatic channel).** Same
sequence shape as R1 with a parallel-plate / cavity readout. Independent
physics; shares only p + dipole density. Tier-1. Cross-locks with R1
(two independent channels, one p* — disagreement above combined bars =
preparation problem, not wave physics).

**R3 — Microstructural order parameter (scattering / microscopy).**
S = ⟨P₂(cosθ)⟩ from SAXS anisotropy or orientation histogram → p-map via
the DECLARED dipole statistics (same declared premise as D1's pair model —
shared premise, disjoint parameters and channel). Tier-1.5: flag the shared
premise at build; it cannot launder a statistics-model error, but it is
independent of every fitted number in the wave formula.

**R4 — Saturation-anchored preparation ladder (p set, not measured).**
Strong-poling asymptote defines the p → 1 scale (self-calibrating: no
absolute standard needed at the top); intermediate preparations ordered by
field/cooling protocol. Tests the STRUCTURAL scaling K_n ∝ p² via
two-state wave-speed ratios (K, C_c, κ_n, ξ cancel at leading order) without
absolute p. Tier-2: fires the scaling leg of kill-(ii)-adjacent content, not
the absolute-formula leg of kill-(iii). The p = 0 foot of the ladder is R0.

**R5 — μ-ratio cross-check (explicitly demoted).** Compute p from r as
consistency readout AFTER adjudication. Agreement = confidence; disagreement
= static-form problem flagged, wave verdict stands as adjudicated with p*.
Never p*.

## 3. Anchors (why the ends are cheap)

- **p = 0:** prediction exactly zero stiffness — qualitative null, no bar.
- **p → 1:** static branch marginally degenerate (sliding mode, priced
  honest) while K_n maximal — maximum wave signal; saturation asymptote
  self-anchors the scale for R1/R2/R4 calibration. Note the asymmetry the
  protocol exploits: statics go soft where dynamics go loud.

## 4. Precision criterion (parametric — numbers at build)

Let O(p) be the tested observable (absolute speed or direction-contrast
Δ = [ω²(x̂) − ω²(ẑ)] / ω̄², uniaxial form banked RD3-3). Let δp be the sealed
half-width from the chosen route. FIREABLE at preparation p iff

|dO/dp · δp| < |O(p) − O_null| − H,

with O_null the kill-(iii) null (non-k² response) and H the hygiene envelope
in the same units. Reading: the p-induced band must sit strictly inside the
kill margin outside hygiene. Because K_n ∝ p², fractional speed error ≈ δp/p
at leading order (pre-stress/gap terms enter at build from the D3 formula).
Consequences frozen here: (a) low-p preparations demand proportionally
tighter absolute δp — prefer mid-to-high p for the formula leg; (b) R0 needs
no δp at all (upper bound p < p_lo with predicted signal below detection
suffices — state p_lo at build from H); (c) if no route meets the inequality
with banked budgets, the protocol returns UNFIREABLE (honest output, §5).

## 5. Recommended sequence + verdict rule

1. R0 null test (p = 0). Kill → lane falsified, stop.
2. R1 + R2 absolute-p at one mid-to-high-p preparation (two channels, one
   sealed p*); R3 if apparatus suits. Check the §4 inequality with banked H.
3. Wave measurement at the sealed preparation; blinded p* insertion;
   kill-(iii) adjudication with the hygiene ledger open.
4. R4 scaling ladder as the follow-up charter (structural p² leg).
5. R5 cross-check after adjudication, never before.

**Verdict rule:** protocol declares FIREABLE iff ≥1 route R1–R3 demonstrates
δp satisfying §4 with banked budgets at a preparation inside the D4 window;
else UNFIREABLE-named stands and the propagation claim stays downgraded per
scope §4. The protocol itself cannot fail — it returns one of the two words
with the inequality ledger attached.

## 6. What this does not do (fences)

No FB-D verdict touched; no number produced; no p measured; no apparatus
chosen (owner/appointee picks R1 vs R2 vs R3 per apparatus at build charter).
κ_n stays priced; flexo-analog acknowledgment (non-uniform-ε terms outside
the uniform formula) travels into the build charter. Drift review requested
on the circularity rule (§1), route tiering (§2), and the inequality (§4).

## 7. Amendments A1–A5 (same-agent review repairs, §§1–6 frozen above)

- **A1 (saturation anchor, closes §1 loophole).** The p → 1 asymptote used
  to calibrate R1/R2/R4 must be identified IN the independent channel (Δn
  or Δε plateau with poling field), never in wave-signal plateau. Locating
  saturation via the wave channel feeds wave data into p* — circular by §1.
- **A2 (R0 zero, reconciles §2 with §4b).** R0 = preparation symmetry PLUS
  an independent-channel zero reading (cheap null on R1/R2 at the same
  preparation). Model p is local; macroscopic unpoling alone does not bound
  residual local correlation, whose weak signal would confound the null.
- **A3 (R3 containment).** R3 is never the sole route: every FIREABLE
  return pairs R3 with R1 or R2 (common-mode statistics error — formula
  and p-map shifting together — is otherwise uncontained).
- **A4 (R4 same-state).** The ladder poles the SAME sample across steps, or
  bounds cross-sample ΔC_c / Δgap inside H at build. Pre-stress and gap
  terms cancel in two-state ratios only under identical mechanical state.
- **A5a (ledger).** Instrument-noise floor joins the budget ledger beside
  H (detection limits enter O_null).
- **A5b (failure legs).** An UNFIREABLE return names the failing leg —
  route precision vs preparation window vs hygiene-dominance — so the
  failure constructs the next attempt.
