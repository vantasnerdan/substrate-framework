---
description: Constructive independent review of C-M5S-012
author: ReviewCMM5S012
created: '2026-08-26T20:15:00+00:00'
updated: '2026-08-26T20:15:00+00:00'
tags:
- substrate-framework
- claim-review
- C-M5S-012
- nonlinear-gravity
category: decisions
confidence: established
status: active
---

## Claim and Positive Role

The reviewed claim gives the nonlinear spherical fate of the accepted certified confined-clock source rather than extending its invalid weak-field potential. In `c=1` mostly-plus areal gauge it states `m'=4*pi*r^2*rho`, `f=1-2*G*m/r`, and `C=2*G*m/r`; a sign-bracketed `f=0` is a necessary obstruction to a globally static horizonless completion of that same source. Applied to C-M5S-004 at the accepted purely induced coupling, it locates the first crossing, maximum compactness, critical coupling and exterior Schwarzschild radius, then establishes that the exact curvature/potential homothetic family also remains trapped at its global compactness minimum.

## Frozen Transaction

The review froze base release `v0.166.0` at `3e3a613`, additive claim `C-M5S-012`, the exact transaction at `local://p245-review-transaction.md`, implementation `src/substrate_framework/nonlinear_clock_gravity.py`, additive root exports, tests, attempts 0001–0003, proposal and memory contract. Accepted dependency propositions used are C-M5S-003's coupling, C-M5S-004's typed density/mass/radius convention, and C-M5S-007's nonlinear escape branch. C-M5S-004/005/007 and issue #170 are affected consumers. The independent reviewer ran no tests, linters or validators; development receipts available at freeze were 29 targeted API/convention tests and 61 targeted consumer tests passing. This record was materialized by Main from the reviewer's `ACCEPT` verdict because reviewer policy prevented direct file mutation; the authoritative independent transcript is `history://ReviewCMM5S012`.

## Strongest Supported Positive Statement

The strongest supported statement is the proposal without narrowing. The exact Hamiltonian constraint supplies the necessary obstruction, the frozen-source primary and independent routes establish the crossing by large margins, and the homothetic scaling plus convex global minimum rules out pure dilation. The statement correctly stops before arbitrary covariant M5 profile relaxation, local pressure/stress, matter Euler equations, black-hole interior dynamics or the optional zero-point mass ledger.

## Evidence Map

The evidence groups have distinct roles and no attachment is credited beyond its proposition.

| Evidence | Proposition established | Role | Bridge to claim | Limit |
| --- | --- | --- | --- | --- |
| `nonlinear_clock_gravity.py`, package exports and tests | Exact mass constraint encoding, finite-volume integration, compactness, crossing, homothetic scaling and analytic soluble limits | exact proof / regression | Implements the mathematical objects used by both routes | Encodes the Einstein constraint; it is not a fresh tensor derivation |
| attempt 0002 | Frozen-source crossing, refinement, critical coupling, exterior radius and homothetic minimum at frozen gates | applicability | Applies exact objects to C-M5S-004's accepted source | Frozen source and homothetic family only |
| attempt 0003 | Independent analytic NumPy source derivatives plus PCHIP/DOP853/Brent corroboration | corroborating subclaim | Removes autograd, radial finite-volume and root-finder common modes | Shares the accepted source definition/root |
| attempt 0001 | `1e-12` absolute identity abort and scale-relative repair pedigree | provenance | Shows the only repaired gate was below measured float64 accumulation scale | No scientific verdict artifact |

## Oracle Audit

The strongest practical oracle is the paired 0002/0003 numerical construction under the exact mass constraint. The reviewer confirmed the mostly-plus areal-gauge signs and factors, and that lapse or pressure equations cannot repair the same failed Hamiltonian constraint after `f` crosses zero. Primary refinement gives crossing drift `1.43e-3`, maximum-compactness drift `2.80e-6`, and trapped signal/error `3.57e5`; independent differences are `8.00e-4` in crossing, `2.11e-6` in maximum compactness and `1.86e-6` in homothetic minimum, with mass-ODE residual `7.85e-6`. Zero density, negative `G`, and `10^-6 G` mutations exercise vacuum, sign and coupling failure channels. The `R^-1` curvature and `R^3` potential ledgers make each pointwise compactness `a/R^2+bR^2`, convex/coercive in log scale; the maximum is convex, so the interior optimizer minimum is global within the continuum radial objective up to the separately refined radial discretization.

## Findings

The review found no blocker and retained three boundary notes as nonblocking scope controls.

| Finding | Direct evidence | Classification | What would extend it |
| --- | --- | --- | --- |
| No blocking finding | Reviewer verdict: `overall_correctness=correct`, confidence `0.97`, `ACCEPT` | Transaction accepted | New contrary evidence would reopen the exact affected claim |
| Package exact surface encodes rather than freshly derives the Einstein tensor | Module plus C-STG-002 convention precedent | Nonblocking; registry must say exact constraint encoding/regression | Independent tensor derivation could add corroboration but is not needed for this accepted convention |
| C-M5S-011 local compactness is unused | Claim exclusions and open de-boxed attribution | Nonblocking frontier | A typed local radiative stress/density map would support a separately stronger mass ledger |
| Arbitrary-profile nonlinear completion remains open | Claim scope and issue #170 | Nonblocking adjacent frontier | Derive curved M5 local stress/Euler equations and continue a deformed profile as a new object |

## Compatibility and Consumers

The claim is native to the accepted mostly-plus, areal-radius, positive-density and C-M5S-004 source conventions. No scalar proxy, fitted mass or equation of state enters. C-M5S-004's weak-field applicability warning is sharpened into an exact nonlinear constraint failure for the same source; C-M5S-005's far-field law remains formally correct but inapplicable to this trapped source; C-M5S-007's nonlinear escape branch is closed for the accepted frozen and homothetic source, while bare-dominated and arbitrary-profile alternatives are unchanged. Issue #170 gains a resolved certified-source nonlinear fate and retains arbitrary-profile work honestly open. Targeted consumer replay passed 61 tests; GitNexus detected low risk and no affected execution process for the additive API/export change.

## Four-Axis Decision

The axes remain independent and the numeric application is not inflated into exact continuum verification.

- Verification: `numeric_evidence`, with exact structural constraint/regression and independent numerical corroboration.
- Review: `accepted`.
- Compatibility: `native`.
- Epistemic: `active` within the frozen-source and homothetic scope.
- Relationship: additive leaf; no accepted claim is superseded.
- Strongest accepted statement: the exact spherical constraint plus refinement-stable trapped-surface and homothetic no-dilation result stated in the frozen transaction.

## Promotion Transaction

Promotion requires moving the adjudicated P245 record into `campaigns/`, adding the exact reviewed statement and C-M5S-003/004/007 edges to `governance/claims.yaml`, pinning `v0.167.0`, rendering generated docs, synchronizing accepted claim/release memory, updating issue #170/#174, and running the impact-selected validation once at the final tree. C-M5S-010/011 remain scope context rather than dependencies.

## Correction Check

No correction was requested, so no correction check is needed. Materializing this record from the reviewer's immutable `ACCEPT` output changes no scientific statement or dependency edge.

## Result and Frontier

`C-M5S-012` is accepted for promotion as drafted. The accepted certified source is not a regular horizonless nonlinear clock at the purely induced coupling, and pure dilation cannot rescue it; the mechanism is a Misner–Sharp compactness crossing by factors near `10^3`, not a failed nonlinear solver. Arbitrary covariant M5 profile relaxation and a local radiative stress map remain adjacent future objects and do not weaken this result.

## Cross-References

See `local://p245-review-transaction.md`, `history://ReviewCMM5S012`, proposal `proposals/P245-nonlinear-self-gravity/proposal.yaml`, attempts 0001–0003, dependencies C-M5S-003/004/007, consumers C-M5S-004/005/007, canonical issue #174, parent frontier #170, and the final promotion receipt recorded with release `v0.167.0`.
