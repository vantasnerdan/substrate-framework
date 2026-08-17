---
description: Independent review and acceptance of C-WLN-001 (massive einbein density, positive branch elimination, Legendre Hamiltonian)
author: vantasner-review and codex
created: '2026-08-17T20:30:58Z'
updated: '2026-08-17T20:45:00Z'
tags:
- substrate-framework
- claim-review
- P227
- worldline
- massive-einbein
category: decisions
confidence: established
status: active
---
# Review of C-WLN-001

## Claim Under Review

The conditional exact statement is: for a supplied nondegenerate
pseudo-Riemannian metric in mostly-plus convention, a timelike velocity
norm `sigma<0`, and exact positive `e`, `m`, `c0`, the density
`L=sigma/(2e)-e*(m*c0)^2/2` has auxiliary equation
`sigma+e^2*(m*c0)^2=0`; its unique positive solution is
`e=sqrt(-sigma)/(m*c0)`; eliminating it gives exactly
`-m*c0*sqrt(-sigma)`. With canonical covector momentum `p_mu`, the
inverse Legendre map is `dot(x)^mu=e*g^(mu nu)*p_nu` and the Hamiltonian
is the pure constraint `e*(g^(mu nu)*p_mu*p_nu+(m*c0)^2)/2`.

## Sourced Inputs

The review read release `v0.159.0`, campaign P227, the canonical
`massive_einbein_ledger` and `einbein_hamiltonian_ledger` APIs in
`src/substrate_framework/relativistic_particle.py`, the existing tests
in `tests/test_relativistic_particle.py`, and the P227 verifier
(`campaigns/P227-einbein-worldline-harvest/verify.py`). Merged PRs #38
and #39 are provenance only. No accepted scientific claim is proposed
as a premise.

## Independence

The review rederives every term without importing
`massive_einbein_ledger` or `einbein_hamiltonian_ledger`. The full
rederivation is documented below; only `sympy` and the repository's
`pseudo_riemannian` helpers (`coordinate_symbols`, `exact_metric_matrix`,
`metric_inverse`) are used for the metric inversion step. The branch
selection and index placement are verified by direct substitution.

### Independent derivation

1. **Density**: `L = sigma/(2e) - e*(m c0)^2/2` (exact, symbolic).

2. **e Euler equation**: `dL/de = -sigma/(2 e^2) - (m c0)^2/2 = 0`,
   which gives the auxiliary constraint `sigma + e^2*(m c0)^2 = 0`.

3. **Positive branch**: solve the constraint for `e > 0`:
   `e = sqrt(-sigma)/(m c0)`. Substituting back:
   `L(e_pos) = -m*c0*sqrt(-sigma)` (the claim's expected square-root
   action). Diff = 0.

4. **Covector momentum**: with `v^i` the contravariant velocity, the
   canonical conjugate is `p_i = g_ij v^j / e` (from
   `dL/dv^i = (1/e) g_ij v^j`). The inverse Legendre map is then
   `v^i = e g^ij p_j`.

5. **Hamiltonian**: `H = p·v - L = e*(p^T g^-1 p) - L`. Substituting
   `v = e g^-1 p` gives
   `sigma_v = e^2 (p^T g^-1 p)`, so
   `L = e*(p^T g^-1 p)/2 - e*(m c0)^2/2`, and
   `H = e*(p^T g^-1 p + (m c0)^2)/2` (the claim's pure constraint).

6. **Non-diagonal metric verification**: on
   `metric = [[-2, 0, 0], [0, 3, 1], [0, 1, 2]]`, the closed form
   `H = e*(p^T g^-1 p + (m c0)^2)/2` is verified to equal the Legendre
   construction `H = p·v - L(v(p))` exactly (diff = 0 after
   `sp.simplify`).

## Verification Status

Exact symbolic substitution proves the auxiliary constraint,
positive-root elimination, square-root recovery, covector momentum
definition, inverse Legendre map, and pure-constraint Hamiltonian.
The independent derivation closes on the same expressions as the
canonical helpers without using them. The claim earns
`symbolic_verified` on the conditional exact statement.

## Sensitivity and Counterexamples

- **Wrong-sign mass-term**: substituting `e_pos` into
  `L' = sigma/(2e) + e*(m c0)^2/2` gives `L'(e_pos) = 0`, which
  does NOT equal `-m c0 sqrt(-sigma)`. Mutation rejected.
- **Wrong index placement (covariant raising)**: replacing
  `v^i = e g^ij p_j` with `v^i = e g_ij p_j` produces
  `H_wrong != H` (nonzero diff) on the non-diagonal metric.
  Mutation rejected.
- **Negative branch**: substituting `e = -sqrt(-sigma)/(m c0)` gives
  `L(e_neg) = +m c0 sqrt(-sigma)` (sign-flipped). Mutation rejected.
- **Massive eliminator at m=0**: would require `0 * infinity`,
  undefined; the massless branch is the correct specialization, not
  the massive eliminator.

All four mutations are rejected by the canonical expressions. The
positive branch and index placement are load-bearing.

## Framework Compatibility

The claim is a conditional exact classical identity. It does not
select a background metric, quantize the particle, identify a physical
species, or import any field-theory content. Positive `e`, `m`, `c0`
and timelike `sigma<0` are explicit assumptions. The zero-mass
specialization is reserved for C-WLN-002.

## Dependency and Consumer Replay

The direct implementation dependency is
`pseudo_riemannian.exact_metric_matrix` and
`pseudo_riemannian.metric_inverse`. They are conditional infrastructure,
not accepted scientific premises; the independent review also inverts the
displayed exact matrix directly.
Consumers are `massive_einbein_ledger`, `einbein_hamiltonian_ledger`,
the canonical worldline tests, and the P227 verifier. All 14 P227
verifier checks pass on the current main.

## Competing Candidate Audit

P227 preregistered three candidates:
- A: scalar-contraction ledger with generic local-metric and
  Lorentz-algebra APIs (selected).
- B: coordinate-function derivation with `euler_equations` (heavier,
  poorer composability).
- C: document-local formulas (no canonical API).

The independent rederivation above is a candidate-B-style oracle; it
agrees with candidate A exactly, validating the chosen ledger without
its help. No empirical comparator exists.

## Four-Axis Decision

The four governance axes are assigned as follows. Verification is symbolic_verified, review is accepted,
`compatibility=compatible_extension`, `epistemic=active`. No
supersession. Scope ceiling is claim content: classical auxiliary-field
identity, not a quantum-mechanical or field-theoretic statement.

## Cross-References

Source campaign: `campaigns/P227-einbein-worldline-harvest/proposal.yaml`.
Source PR: PR #60 (merged). Review request:
`campaigns/P227-einbein-worldline-harvest/review-requests/C-WLN-001.md`.
Canonical goal: issue #59.
Executable corrective review:
`campaigns/P227-einbein-worldline-harvest/reviews/independent_worldline_lorentz_review.py`.
