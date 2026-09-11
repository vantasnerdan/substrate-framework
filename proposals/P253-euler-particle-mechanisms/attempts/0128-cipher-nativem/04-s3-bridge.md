# S3 F²-bridge paper (FROZEN F2 spec in 03-prefire.md) — verdict: DEAD, residue banked

Domain (F2): Ω_δ = R³ minus closed tubular neighborhoods (radius δ ~ a) of
filament cores; (α,β,φ) smooth on Ω_δ; vorticity δ-supported on filaments
(distributional core terms stated, none dropped); Lin constraint weak (test
variations compactly supported in Ω_δ). δ-gate: conclusions below are FORM
statements (algebraic structure), δ-independent by inspection; any future
quantitative bridge repeats at δ/2 per F2.

## Setup
Clebsch: u = ∇φ + α∇β on Ω_δ. Euler kinetic action S_E = ∫dt∫_Ω ½|u|² (+
pressure/incompressibility constraints, A-sector-neutral — stated, not
carried). Axisymmetric data assumed for concreteness; conclusions are
algebraic (hold generally on Ω_δ).

## Reading 1 (gauged Clebsch): A auxiliary — gauge-dynamics DEAD
Gauge the β-shift symmetry: ∇β → ∇β − κA. Then
|u|² = |∇φ + α∇β − καA|² = |..|²₀ − 2κα(∇φ+α∇β)·A + κ²α²|A|².
A-sector: linear source −2κα(∇φ+α∇β)·A plus MASS term (κ²α²/2)|A|² —
ALGEBRAIC in A (verified by sympy expansion: A² coefficient (k²α²/2) per
component; zero ∂A-dependence anywhere). No F² kinetic term is generated at
any order of the Euler–Clebsch kinetic term: A is auxiliary (constraint),
integrating it out yields current–current interaction, never propagating
gauge dynamics. Falsifier 1 FIRES exactly as the sketch designed: "E,B are
gauge artifacts" — the sketch kills its gauge-dynamics version cleanly.
(No Lin subtlety rescues this: Lin constrains (α,β) variations, not A's
algebraic status.)

## Reading 2 (identification u↔A, ω↔B): Maxwell form absent — Maxwell-dynamics DEAD
Under the R-EM1 kinematic identification: Euler kinetic ∫½|u|² supplies an
"electric-type" (∇φ-ish) term ONLY; there is no −|B|² magnetic term with
independent coefficient (the action has one sign, one coefficient — the
Maxwell F² = E²−B² relative minus is underivable). Helicity ∫u·ω IS the
Abelian Chern–Simons 3-form ∫A∧dA on the Clebsch data (verified structure:
w = ∇α×∇β, h = ∇φ·w = Jacobian det(φ,α,β), α-independent) — FIRST-order
topological, not second-order Maxwell. So the native gauge sector is
Chern–Simons/BF-TYPE, never Maxwell-type, by direct inspection of the action.
Second falsifier FIRES (Lin/propagation): A has no propagation equation from
S_E; any "waves" would need an added kinetic principle (named missing, not
smuggled — cf. 0043's ban, which this confirms from the action side).

## Residue (the honest positive — upgrades S1/S4, converges with canon)
Euler-native EM is TOPOLOGICAL (CS/BF-type, integer-valued charges), not
Maxwellian. This (a) converges with the canon BF ledger (compactness missing
≠ BF wrong — the FORM is BF, the compactness source is the open problem);
(b) converges with S1's integer-linking charge (D1) and S4's relational
linking; (c) reframes the lane: the missing construction is COMPACTNESS
(integral level from Euler data) + PROPAGATION (second-order dynamics),
not coupling form. Maxwell F² needs new physics; BF-type charge does not.
Core-dominance check (F2 pre-registered): distributional core terms carry the
helicity/CS density's singular part — filaments localize the topological
charge (S3→S1 reduction CONFIRMED as convergence, per F2 rule).

## Verdict
S3-as-Maxwell-route: DEAD under both readings (mechanisms named: auxiliary-A;
absent-F²/CS-instead). S3-residue: native EM = topological/BF-type —
BANKED as lane direction. No imports used; no verdicts beyond the bridge.
Next (chartered order): S4 tangle model, last — or lane synthesis if
shepherd calls it.
