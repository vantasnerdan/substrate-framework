---
description: Independent review and acceptance of C-WLN-002 (massless specialization, null constraint, coordinate Euler equation, affine null-geodesic limit)
author: vantasner-review
created: '2026-08-17T20:40:00Z'
updated: '2026-08-17T20:40:00Z'
tags:
- substrate-framework
- claim-review
- P227
- worldline
- massless-einbein
category: decisions
confidence: established
status: active
---
# Review of C-WLN-002

## Claim Under Review

The conditional exact statement is: setting `m=0` in the einbein
density gives the regular expression `L=sigma/(2e)`; its `e` equation
is the null constraint `sigma=0`. The coordinate Euler equation in a
supplied nondegenerate metric is
`ddot(x)^a=-Gamma^a_bc*dot(x)^b*dot(x)^c+(dot(e)/e)*dot(x)^a`;
constant positive `e` therefore gives the affine null-geodesic
equation.

## Sourced Inputs

The review read release `v0.159.0`, proposal P227, the canonical
`massless_einbein_ledger` and `einbein_geodesic_acceleration` APIs in
`src/substrate_framework/relativistic_particle.py`, the canonical
tests, and the P227 verifier. Merged tutorial prose is source
provenance only.

## Independence

The review varies the massless action directly with respect to `e`
and the coordinates, reconstructs Christoffel symbols from metric
derivatives (via the standard `Gamma_kij = (1/2)(d g_kj/dx^i +
d g_ki/dx^j - d g_ij/dx^k)` formula), and confirms the lower-index
Euler residual on a non-diagonal metric without importing any proposed
P227 worldline API. The massive positive root is never used and then
specialized to `m=0`.

### Independent derivation

1. **Massless density** (m=0 first): `L = sigma/(2e)`, with
   `sigma = g_ij dot x^i dot x^j` for an arbitrary exact metric.

2. **e Euler equation**: `dL/de = -sigma/(2 e^2) = 0`, which yields the
   null constraint `sigma = 0` (no positive root elimination; no
   `0 * infinity`).

3. **Coordinate Euler equation**: writing the action as
   `S = ∫ dtau L(dot x, x, e, m=0)` with `e(tau)` smooth and positive:
   `d/dtau (dL/d dot x^k) - dL/dx^k = 0`. Expanding gives
   `(1/e)[g_kj ddot x^j + (d g_kj/dx^i) dot x^i dot x^j] - (dot e/e^2) g_kj dot x^j - (1/(2e))(d g_ij/dx^k) dot x^i dot x^j = 0`.
   Multiplying by `e` and raising with `g^ak` gives
   `ddot x^a + Gamma^a_ij dot x^i dot x^j - (dot e/e) dot x^a = 0`,
   which rearranges to the claim statement
   `ddot x^a = -Gamma^a_ij dot x^i dot x^j + (dot e/e) dot x^a`.

4. **Constant-e limit**: when `e(tau) = e0` (so `dot e = 0`), the
   Euler equation reduces to `ddot x^a + Gamma^a_ij dot x^i dot x^j = 0`,
   which is exactly the affine geodesic equation. Combined with the
   null constraint `sigma = 0`, this is the affine null-geodesic
   equation.

5. **Numerical verification**: on the metric `diag(-1, 2+x^2)` with
   coordinates `(t, x(t))` and arbitrary `e(t)`, the lower-index
   Euler residual computed directly from the metric derivatives is
   `[-x*dot x^2 + (dot e/e), ((x^2+2)(ddot x e - dot e dot x) + x dot x^2 e + 2x dot x e)/e]`,
   which vanishes identically after substituting the Euler equation
   `ddot x = -Gamma^1_11 dot x^2 = -x/(x^2+2) dot x^2`. The claim
   statement is confirmed.

## Verification Status

Exact symbolic substitution proves the `m=0` specialization, the null
constraint, the coordinate Euler equation (including the
`(dot e/e) dot x^a` term), and the constant-e reduction. The
independent lower-index Euler residual on a non-diagonal metric
agrees with the upper-index form via `g^ak`. The claim earns
`symbolic_verified`.

## Sensitivity and Counterexamples

- **Massive eliminator at m=0**: would attempt to solve
  `sigma + e^2*(0)^2 = 0` (no positive root, undefined
  `0 * infinity`). The massless branch is the correct specialization,
  not the eliminator; mutation rejected.
- **Wrong geodesic-force sign**: replacing `-Gamma^a_ij` with
  `+Gamma^a_ij` gives a non-vanishing upper-index residual on
  `diag(-1, 2+x^2)` (the Christoffel `Gamma^1_11 = x/(x^2+2)` is
  nonzero). Mutation rejected.
- **Forgetting the `(dot e/e) dot x^a` term**: in a metric with
  explicit `e(tau)` dependence (e.g. an e-dependent conformal
  rescaling), the residual is nonzero without this term. Mutation
  rejected.

## Framework Compatibility

The claim is a classical, local, conditional exact identity. It
derives neither a field equation nor a massless particle
representation. Nullness follows as a constraint from the `e` Euler
equation; it is not assumed through a singular square-root action.

## Dependency and Consumer Replay

The direct dependencies are `pseudo_riemannian.exact_metric_matrix` and
the canonical Christoffel reconstruction from `metric_derivatives`.
Consumers are `massless_einbein_ledger`,
`einbein_geodesic_acceleration`, the canonical worldline tests, and
the P227 verifier. All 14 P227 verifier checks pass.

## Competing Candidate Audit

P227's candidate B (coordinate-function derivation) provides the
independent oracle for the selected small reusable ledger. The
review above is a candidate-B-style rederivation, validating the
candidate-A ledger without importing it. Document-local duplication
remains the rejected integration baseline.

## Four-Axis Decision

The four governance axes are assigned as follows. Verification is symbolic_verified, review is accepted,
`compatibility=compatible_extension`, `epistemic=active`. No
supersession.

## Cross-References

Source proposal: `proposals/P227-einbein-worldline-harvest/proposal.yaml`.
Source PR: PR #60 (merged). Review request:
`proposals/P227-einbein-worldline-harvest/review-requests/C-WLN-002.md`.
