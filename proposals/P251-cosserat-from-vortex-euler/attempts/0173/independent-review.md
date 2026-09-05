---
description: Independent review of the actual prepared Euler mean and action transfer C-CST-013
author: construction_review
created: '2026-09-05'
updated: '2026-09-05'
tags:
- substrate-framework
- claim-review
category: decisions
confidence: established
status: completed
---

# C-CST-013 independent scientific review

## Claim and positive role

The proposed claim is the actual fixed-time, response-specific second
directional Bloch jet and its uniform-period transfer to a selected
finite stationary constant-curl insertion retaining the local optical
object. I accept the statement frozen in0173/README with the explicit
mean-zero background, transverse preparation, norms and ordering in
0170/prepared-mean-jet.md. The useful new result is not simply existence
of a Taylor coefficient: it eliminates the apparent inverse microscopic
gap from the OBSERVED coefficient, permitting a finite same-field
positive physical acoustic action while retaining prescribed local
optical/EPS margins.

This is one substantive claim review. I supplied no derivation or
implementation input to0170 or its new response estimate before taking
the review role. I authored earlier unchanged Euler-Fourier utility
work used by its algebraic verifier and earlier, separately declared
campaign inputs including0145's twist attachment. Those unchanged
inputs are not represented as independently re-proved here. I did not
author0151,0153,0154 or the new0170 derivation. I have not commissioned
another reviewer, rerun a scientific check, or edited those sources.

## Frozen transaction

The transaction is based on v0.177.0 and checkpoint dfe495c. The central
0173 registration passed265 claims before the review. C-CST-013 is a
new proposed claim; C-CST-011 and the existing release remain unchanged.
The reviewed raw artifacts have these SHA256 pins:

- 0170/prepared-mean-jet.md: 6796d42c37da53843310199a21c45e00a98ea69e4fb000b400cd5bc8b1455cea
- 0170/verify.py: e688d1d9ede8aa8aec49b12d94e895ce39ef90d317ca08fe710b1dda09096022
- 0170/receipt.md: 88215fb55269901426a0de98edf047dcd7fa9a0a69fd0bb60a30c43f64c24bca
- 0153/same-field.md: 0adbb25c02c81102fffad4cac3886038c54a8511a17ded009fd1e8613325d14d
- 0154/coherent-physical-action.md: 09281417184a987583874e3c59cd99ad2326d22a5c230ebf150821d1b92df8d3
- 0151/one-wave-response.md: 5a012db4150cd5028e4bf262533358fea513a7e5504e179846982f13b43ed9bb

The scientific source and verifier hashes agree with the native receipt.
The source, actual verifier and final19/19 process0 stdout were read.
The importable extraction and its directly affected tests are the
remaining promotion implementation transaction, not a missing analytic
step in this result. No broad replay was needed for this read-only audit.

## Strongest supported positive statement

For a smooth stationary mean-zero periodic Euler field, fixed transverse
Kelvin-displacement D and independently specified common-velocity V,
the complete first and second fixed-time directional spatial jets of
the physical mean X,m are determined by0170 equations(3),(8). The
forced first-cell problem is the actual full Euler/Lin problem; its
second unknown cell is absent from this physical averaged row by an
exact divergence/pressure identity.

For the particular bounded constant-curl quadrature construction in0153,
normalized H2 convergence to the circular one-wave background implies
uniform convergence of these response coefficients on each fixed time
window, with constants independent of the growing quadrature period.
After averaging the actual coherently prepared phase actions over the
whole-field law, a selected finite field and sufficiently small nonzero
macro carrier have strictly positive physical-chart mass and acoustic
second-jet stiffness. The same selected field retains the separately
fixed C-CST-011 optical and material-tube margins at their actual scope.
Every physical/canonical momentum and time-dependent action connection
in the source remains part of this statement.

No correction or narrowing of that positive statement is required.

## Evidence map and the exposing oracle

The analytic calculation is the strongest practical oracle. Its
algebraic anchors and imported local construction have distinct roles.

| Evidence | Contribution | Role and limit |
| --- | --- | --- |
| 0170 sections1--3 | Actual preparation, harmonic convention, first-cell forcing, elimination of second cell, physical mean current | Exact response proof for general smooth stationary mean-zero periodic Euler |
| 0170 section4 with0153 normalized norm bounds | Uniform-period H1 forced Euler/Lin comparison and stress coefficient estimate | Analytic transfer proof; not full Bloch-operator C2 convergence |
| 0170 section5 with0151/0154 | Ordered finite-field selection, whole-action averaging and Wronskian chart | Same-action positive physical second-jet result on fixed time |
| 0170/verify.py and captured19/19 stdout | Independent pressure signs, exact ABC example, exposed current difference, nonautonomous Wronskian identities | Exact algebraic corroboration, not a simulated PDE or numerical spectral proof |
| C-CST-011 and0153 local transfer | Retained actual finite packet, material tag, full pressure and tube margins | Declared local optical compatibility input, not a new joint continuum theorem |

I checked the following load-bearing edges independently from the
equations, rather than accepting the pass tally as their proof.

1. On the directional transverse bundle, the slow harmonic is P_kappa.
   The exact identity P_K grad f=-ik P_K kappa f includes it. Applying
   this to grad(U0.grad p) differentiates an actual gradient source,
   not a naked projector whose derivative really has a 1/q singularity.
   Expanding A_K² U0 gives i k (A a)U0 and A a=-kappa.grad p. These
   two terms supply both forcing pressure columns in equation(3).

2. Kelvin preparation is not interchangeable with arbitrary mean
   velocity data. Expanding its actual material operator gives
   chi_t(0)=-P[aD+kappa(u.D)], with zero microscopic mean. Common V
   has its independent circulation preparation. This is compatible with
   eta0=D+tV and does not assert that both columns lie on one fixed leaf.

3. For w=chi_t+Achi-(Du)chi, direct differentiation gives
   w_t+P[A w+(Du)w]=P[chi_tt+2A chi_t+(A²+Hess p)chi].
   The identity A(Du)+(Du)²=-Hess p fixes the sign. Thus the source's
   second-order equation really factors into full linear Euler followed
   by Lin transport, with order-zero pressure and no assumed bounded
   second-order propagator. Divergence and microscopic means are preserved.

4. The second divergence equation is div zeta=kappa.chi. Periodic
   integration gives <Hess p zeta>=-<grad p(kappa.chi)>. Also
   <a Achi>=-<(A a)chi>=<(kappa.grad p)chi>. In the complete averaged
   A_K² term, <A(a chi)>=0, leaving the displayed second pressure row
   with its correct sign. This removes zeta without solving it. The
   harmonic pressure multiplier is killed by P_kappa, not declared zero.

5. Exact Lin reconstruction and integration by parts give
   m=U_t+ik<a eta-u(kappa.eta)>. Differentiating its second coefficient
   changes 2<a chi_t> into <a chi_t+u(kappa.chi_t)> in R_u. Both
   pressure terms survive. The ABC verifier's nonzero pressure and
   current rows expose omission of these terms; a constant-pressure
   elementary-wave regression alone would not do so.

6. In normalized L2/H1, the Euler/Lin energy estimate uses transport
   integration by parts, fixed C3 bounds and P commuting with spatial
   derivatives with norm one. Subtracting the smooth reference columns
   places h.grad w0 and (Dh)w0 in H1 at cost ||h||H2, without a
   period-dependent Sobolev embedding. The Lin equation bounds chi_t
   in L2 using chi in H1. Pressure products are controlled by normalized
   L2 derivatives and fixed C2 bounds, not by assuming h tends to zero
   in L-infinity. Pairing the resulting columns with the actual stress
   row proves the stated uniform coefficient estimate. This argument
   genuinely handles arbitrarily low insertion self-beats.

7. The finite field is chosen after its local optical margins and
   global total-variation bound, but before the small macro carrier.
   The exact identity ||h_N||H2²=(1+lambda²)² sum|a_j|² uses the
   Fourier H2 convention explicitly fixed in0153. Quadrature can make
   it small without weakening the local insertion. The analytic
   k-neighborhood and higher remainder may depend on this selected
   finite period; the result does not exchange these limits.

8. The initial constrained momentum is
   pi_D=-i rho P_K[K(u.D)] and pi_V=rho V. Zero microscopic mean
   gives the exact normalized pairing rho J. If using a real periodic
   encoding, choose the final small carrier below half the selected
   minimum nonzero microscopic frequency to avoid its twice-carrier
   resonance; this is already available within the source's final
   sufficiently-small-k selection. No added mean or rigid inertia is
   needed. Whole-field action averaging uses shared INITIAL phase
   coefficients and each realization's own evolved reaction, before
   any physical chart is eliminated.

9. For actual physical rows f,g, the Wronskian formulas in(12) follow
   by pulling back the conserved phase symplectic form and its moving
   connection. They yield M=rho/W and
   M X_tt+M_t X_t+K X=0. At the one-wave second jet W=1+O(k³);
   a coefficient perturbation of size epsilon leaves an actual
   O(epsilon k²) mass/connection term. It is retained in(13), as is
   rho(1-1/W)m. Small epsilon then small k preserve both positive
   mass and stiffness uniformly on the selected fixed window.

The tests inspect these defining objects rather than literal pass
booleans. The nonautonomous scalar example checks the chart identity,
not a supplied oscillator model for Euler. No small-ratio numerical
enclosure or new solver run is required: all signs are separated from
their explicit analytic errors by the ordered construction.

## Findings and compatibility

There is no load-bearing scientific blocker in this claim boundary.
In particular the period-uniform estimate is response-specific as
advertised, and the action is a pullback of actual coherently prepared
histories, not an off-shell local law inferred from a force coefficient.

The claim is a compatible extension of incompressible constant-density
Euler with full pressure, physical mean, specified circulation classes
and the declared whole-field law. It introduces no fitted sound speed,
new spin identity, or change to C-CST-011. A stationary background does
not make every conditioned packet realization stationary or its moving
coefficients autonomous; the cited C-CST-011 scope and0170 already
retain that distinction. The finite positive optical density is not
claimed uniform as the quadrature volume tends to infinity.

Generic coupled spin/current closure and acoustic-time homogenization
are active parent achievements, not debt hidden inside this statement.
The useful next consumer is the actual same-background optical/acoustic
phase-current join with the newly licensed mean second jet. This review
does not impose its unfinished conclusions as new gates on C-CST-013.

## Four-axis decision and promotion transaction

The scientific decision is acceptance of the stated new theorem.

- Verification: symbolic_verified, with the analytic PDE estimate as
  load-bearing proof and exact algebra as exposing corroboration.
- Review: accepted, one independent substantive pass.
- Compatibility: compatible_extension.
- Epistemic: active upon registry/release promotion.
- Relationship: additive response/action theorem using unchanged
  local optical construction, not supersession or full parent closure.

Promotion still requires the root's reusable extraction and direct unit
tests, registry/adjudication/release materialization, generated records
and impact-bounded validation. Those implementation steps should reuse
the current native scientific receipts while their sources stay fixed.
No full-suite rerun or second scientific review is requested.

## Correction check and result

Correction check: not needed; no source correction was requested.
The strongest positive scope in0173/README is supported by the raw
proof and evidence. This accepts C-CST-013 scientifically, while leaving
the overall PR199/P251 objective active exactly as the parent contract
requires.

Signed: construction_review, independent non-author of0170,
2026-09-05. Reviewed only this frozen transaction and its declared uses.
