---
description: Certify the confined-clock full-band kinetic-normalized spectrum about the committed R=12 family-S root via independent quadrature families, then compute the blocked zero-point mass shift with an honest error budget
author: ox-alpha
created: '2026-08-25T19:30:00+00:00'
updated: '2026-08-26T12:30:00+00:00'
tags:
- substrate-framework
- campaign-proposal
category: proposals
confidence: exploratory
status: active
---

## Question and Positive Deliverable

C-M5S-009 delivered radiative stability but blocked its quantitative half on one named
construction: a certified complete kinetic-normalized spectrum table with
independent-discretization agreement at a declared scale-relative tolerance across the
full band. This campaign constructs exactly that object about the committed R=12
family-S window root (order-16 basis), certifies the fluctuation frequencies whose
quadrature sensitivity was recorded as gate G2 failure in P243 attempt 0008, and
computes the renormalized zero-point mass shift delta-E = (1/2)*sum_i omega_i over the
certified propagating set within the committed model class, carrying an itemized error
budget.

## Base Release and Provenance

Base release v0.165.0, main tip dbf44fa (merge of PR #171, small-ratio-numerics
hardening). Accepted inputs read at source: C-M5S-001..009 statements; C-M5S-006 root
provenance for the R=12 family-S window root with stiffness block scale s_b =
2.982251210281484; C-M5S-009 missing-construction wording; attempt artifacts of
campaigns/P243-clock-sourced-induced-coupling/attempts/0008; committed infrastructure
proposals/P240-m5-kinetic-axis/attempts/0041 cpu_energy.py and solve_radial_1d.py plus
attempts/0042 largeR-roots.json; canonical modules under src/substrate_framework/.

## Invariants, Conventions, and Allowed Imports

The background is FROZEN once and never re-solved; refinement varies only the
discretization of the fluctuation operator about identical field data. Committed
functional conventions are preserved verbatim: envelope factors x^2(1-x^2), (1-x^2),
x^4(1-x^2); potential V(S) = -0.5 tr S^2 - tr S^3 + (tr S^2)^2 + 0.5; measure
2*pi*r^2 dr dmu; Frobenius product; fixed director frames. The tangent-channel family
is an exact kinetic null deflated analytically; the kappa projector-current term stays
excluded with its recorded bias direction; no fitted constants enter any statement;
the Route-2 no-go is respected. Numerics follow the hardened small-ratio-numerics
skill: pinned BLAS threads, compensated reductions feeding quoted quantities,
observed-order checks before extrapolation, per-quantity budgets, zero-mode gauges.
Allowed imports: C-M5S-001/002/006/009 as accepted context and protocols; P240 0041
and 0042 helper APIs and roots as source-audited committed infrastructure; canonical
package modules read-only.

## Candidate Preregistration

Route A is the committed modal machinery itself with autograd-exact derivatives on
Gauss-Legendre nodes. Route B candidates registered: Chebyshev-family quadrature via
rule injection into the verbatim committed machinery (selected once FD ladders proved
derivative-discretization cross-checks floor at their own truncation), with the FD
radial route retained as supplementary coarse consistency evidence.

## Selection Criteria and Blinding

Route selection is structural: the selected second family shares no node set and no
weight rule with route A while keeping every other construction byte-identical, which
isolates the quadrature rule as the only variable between families. Gates were
declared before any route-B value existed. G0 transfer within 1e-6 relative per rung.
G1 symmetry and PSD hygiene at 1e-12 and 1e-11. G2 repaired per-band scale-relative
omega agreement with bands by stiffness decades of s_b and tolerances soft 5 percent,
mid 2 percent, stiff 1 percent of band-local maximum omega. G3 mutations: coefficient
mutation moves the pencil, channel rays pairwise distinct above 1 percent, R10-root
sensitivity in the scale-free form recorded in the attempts table. G4 artifact
adjudication follows C-M5S-006 where applicable. G5 certification margin omega at
least ten sigma with itemized budgets. G6 cross-family entry exactness with its
measured floor recorded rather than assumed.

## Proposed Claim Delta

C-M5S-010 (numeric_evidence): certified full-band kinetic-normalized spectrum table
about the committed R=12 family-S root by cross-family agreement under gates G0-G6,
with per-mode classification and itemized budgets. Dependencies: C-M5S-006.
C-M5S-011 (numeric_evidence): renormalized zero-point mass shift delta-E =
(1/2)*sum_i omega_i over certified propagating modes within the committed order-16
sector about that root, discharging the missing construction that C-M5S-009 names.
Dependencies: C-M5S-010 and C-M5S-009. Identifier collision search executed 2026-08-25
over governance, campaigns, proposals, and memory: all three identifiers free.
C-M5S-009 itself is neither rewritten nor challenged; its stability half stands and
its blocked half is completed downstream.

## Implementation and Oracle Plan

Attempts are append-only under proposals/P244-clock-full-band-spectrum/attempts/.
All verifiers run with PYTHONPATH=src against committed helper modules; stdout is
captured into the attempt directory on first execution.

## Attempts (2026-08-25 session)

Ten append-only attempts separate defect discovery from repair and each
verdict below is backed by captured stdout and JSON artifacts in its directory.

| Attempt | Content | Verdict |
| --- | --- | --- |
| 0001 | Route-A ladder on committed machinery; exposed non-converging M while H saturated at 1.4e-15 | VERIFIER DEFECT FOUND: P243-0008 kinetic_functional reduced as T=4(sum w)(sum rho) instead of per-cell weighting; positive-scalar scaling preserves all scale-free stage-2 verdicts; P240 inertia path verified clean; no accepted claim invalidated |
| 0002 | Corrected per-cell reduction; M saturates at 1.5e-15 like H; mp congruence bug produced a spurious uncertified population | implementation defects recorded; G0/G1 pass all rungs |
| 0003 | mp congruence repaired via Z^T H Z with L Z = I; band labels from per-mode stiffness Rayleigh quotients | all 32 kept modes certified; dE = 72.58859646 plus minus 3.7e-7 |
| 0004 | Scale-free sensitivity gate (R10 max relative shift 0.253 above 0.1); spectrum-table.json artifact | 7 of 7 CHECKS PASS |
| 0005/0006 | FD radial-discretization pencil ladders, second and fourth order stencils | cross-route agreement floors near median 4 percent because pointwise channel oscillation limits grid derivative accuracy regardless of stencil order; FD route demoted to supplementary consistency |
| 0007 | G6 quadrature-exactness by rule injection into the verbatim committed machinery using interior second-kind Chebyshev nodes with solved weights: M entries machine-exact at 3.1e-14; H entries floor at 4.5e-7 from algebraic pole content in the azimuthal channel; pencil 1.04e-7; delta-E replicates to 9.1e-10 | mechanism named; family floor folded into budgets |
| 0008 | Final certified table with cross-family budgets | 4 of 4 CHECKS PASS; all 32 kept modes certified; dE = 72.58859646 plus minus 3.9e-7 RSS and 9.2e-7 linear; the soft band is EMPTY and every certified frequency is O(1) or above with omega_min = 1.0464, independently reinforcing the C-M5S-006 no-fourth-species verdict at frequency level |
| 0009 | PR-172 review repairs (erratum exclusions for the preregistered 0007 pencil-gate miss, reviewer-session identifiers, corrected full-mode validation receipt) plus the preregistered Casimir-style empty-window subtraction: verbatim machinery at the all-zero 48-coefficient background, both quadrature families, pairing-free level sums at N=24/32/40 | BLOCKED WITH MECHANISM: clock side reproduces C-M5S-011 bit-exactly (S32 = 145.17719291997776, diff 0.0) and the mutation gate passes (shift 4.5e+29 relative), but the empty window certifies 0 of 32 kept modes -- its kinetic metric carries a resolution-stable soft direction M[16] = 1.910066e-07 (130x softer than the clock's first propagating kinetic eigenvalue 2.503e-05), so raw pencil frequencies are roundoff junk near 4e41 and no certified empty spectrum exists to subtract. The bare delta-E = 72.58859645998888 stands as the defined zero-point quantity within the committed model class; a renormalized comparison needs a construction that does not rely on kinetic normalization of window-only directions (next candidate, not yet designed) |
| 0010 | De-boxing ladder Delta_E(R) about the window-family roots R10/R12/R14, verbatim machinery with per-rung G0 transfer (worst 7.5e-12 relative), all-mode certification and mutation gates; preregistered bands in the manifest before execution | MIXED per frozen bands: D = 76.369751262 / 72.588596460 (R12 bit-exact C-M5S-011 regression) / 72.536205371; ratios 1.0521 / 0.99928 miss the symmetric plateau window at the inner rung by 0.0021; p = -0.153. Non-monotone shape excludes box decay AND truncated-vacuum growth; per-rung linear half-sum budgets (5.98e-04 / 1.13e-03 / 1.95e-03) sit three orders below the R10 elevation of 3.78 so numerics cannot explain it; omega_min rises outward (1.010 / 1.046 / 1.089), opposite to box quantization. Mechanism: decaying wall-proximity transient atop an apparently flat intrinsic core. Consumer-mass attribution stays formally open; outward extension preregistered as 0011. REVIEW ERRATUM (Reviewer244b): the first draft of this row and of the manifest revision quoted R10's budget as ladder-wide; G0 references for R10/R14 are in-session regression anchors partially corroborated by largeR-roots.json energies to 6.2e-09 / 4.1e-08 relative; the executed script docstring's '4e-13' is 4.1e-10 absolute (7.5e-12 relative); executed artifacts left untouched |
The soft-direction value was subsequently confirmed digit-stable across a
further independent ladder: M[16] = 1.9101e-07 identically at 48x24, 96x48,
144x72 and 192x96 (16x the nodes), while adjacent null-scale eigenvalues
bounce at the 1e-15/1e-16 roundoff scale -- distinguishing converged operator
structure from under-resolution noise.

The scheme-spread frontier item of issue #170 was quantified from accepted canon with
no new claim: R(z) = J_smooth/J_sharp equals one exactly at z=0 with leading
difference EulerGamma*z; the unit-mass value is 2 e K_1(2)/(1 - e E_1(1)) =
1.8837726; R(4) = 15.61; regulator dependence for massive species reopens gradually
and is unbounded.

The PR-172 review also confirmed every headline number digit-exactly and the
full pytest suite (2440 passed); its three governance findings are repaired
in this session's commits.

## Blast-Radius Consistency Audit (2026-08-26)

User directive: examine the full blast radius of delta-E = 72.58859645998888
and resolve framework-wide consistency. Inventory and classification of every
consumer of the clock sector's energy/mass numbers:

| Consumer | Quantity used | Classification under the open attribution question |
| --- | --- | --- |
| C-M5S-004 BVP consumer Phi(R); mass identity M = E_static(mu-avg) = 54.70900884959007 | classical static density, declared convention excluding the fixed-J constraint term | Scope-explicitly-classical claim; unchanged. Whether its sourced mass should gain +Delta_E is exactly the open attribution question -- flagged as frontier-dependent, not an erratum |
| C-M5S-004/007 strong-coupling wall G_total M / R = 213.397 at B = 0 (C-M5S-007 verdict on the C-M5S-004 BVP solution) | same BVP mass with G_total = 46.807 | INVARIANT under attribution: naive full inclusion gives about 496.6, still >> 1 by over two orders; the dichotomy verdict is untouched (the B_min weak-field bound scales linearly with the chosen M, so its numeric values shift under inclusion, but the induced-versus-bare dichotomy does not) |
| C-M5S-005 pairing F = G_total M1 M2 / d^2 | per-clock mass twice | Sign and d-dependence unchanged; magnitude would rescale by (1 + delta)^2 if the shift is included; magnitude flagged attribution-dependent |
| Census thresholds s_b and stiffness band edges (C-M5S-006) | classical Hessian block scales about the root | Already fluctuation-sector objects computed about the classical background; no zero-point input exists at their scope; unchanged |
| Lambda^2 = E_U - E_S = 0.26847204181661866 induced-gravity cutoff | eigenvalue splitting about the root | Unchanged; no zero-point sum enters its definition |
| Issue #170 step-1 figure E = 50.45 stable for R in [~8, ~34] | mu-independent fixed-J branch (N=20, R=8, P240 phase1 ladder) | Provenance clarification only: a distinct root family from the window-family E_cl = 55.10418278043526; recorded here so future readers do not read them as one object |

Robustness headline: every accepted gravity-consumer VERDICT survives both
outcomes of the attribution question -- the strong-coupling wall is 213.4 or
about 496.6, i.e. >> 1 either way. What awaits attribution is only the precise
coupling magnitude and force magnitudes, not any accepted verdict.

## Debt Ledger

Empty inside the proposed scopes. Declared hypotheses (frozen background, excluded
kappa term, order-16 sector scope) and measured mechanisms (algebraic endpoint content
of H densities) are recorded hypotheses or mechanisms, not debt.

## Review and Promotion Plan

One claim-level review per proposed claim by a reviewer without authorship, against
raw attempt artifacts. Then registry entries, release manifest bump, rendered docs,
memory synchronization, scoped validation. The pull request is opened, not self-merged,
per the standing user directive of 2026-08-23.

## Done Gate

Done when both claims pass individual review, enter the registry and a pinned release,
docs and memory agree, and issue #170's frontier item 2 comment carries the certified
result. If review refutes a statement, the objective stays open with the next
registered repair named.
