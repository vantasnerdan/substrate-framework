# 0147-sage-s3exact — SCOPE (frozen before any derivation)

Charter: shepherd S3-EXACT identities paper (synthesis order (b), 0143 §3).
Noether/action exact identities for the charge observable, on the F²-bridge
computation's domain and F2 core treatment. D-08-exempt route. PAPER ONLY;
drift firewalls; no builds, no numerics beyond symbolic algebra receipts.

## Inputs consumed at source (never recomputed as if new)

- cipher 04-s3-bridge.md: F²-bridge computation, FROZEN F2 spec reference,
  both Maxwell readings DEAD (auxiliary-A; absent-F²/CS-instead), R-a/R-e
  SymPy banking (A-sector algebra; helicity = Jacobian, α∇β term drops).
- cipher 03-prefire.md §F2: domain rule (Ω_δ tubular excision, distributional
  core terms printed-not-dropped, weak Lin, δ-gate).
- atlas 0139-atlas-s3gaps: G1 pressure-as-A₀ CLOSED-negative (gauge-fixing
  not Gauss); G2 core terms printed (H_sing = ΣΓ²Slk + 2ΣLkΓΓ, CWF class);
  Kelvin-gauge coupling = named reopen price (heavy).
- 06-synthesis-spec: charge = framing/self-linking integer L; P2 phase-charge
  lock; missing-5 = dynamical charge-observable (open).
- drift ledger D-08 boundary clause: EXACT analytic identities stay
  equalities — this paper operates inside that license; no inequality
  predicate is claimed anywhere in it.
- drift assumption-hunt F1/F2/F3: the S3 kill-assumptions this paper must
  respect when consuming the kill (F1 printed by G2, F2 closed negative by
  G1, F3 strong — no axisymmetry assumption anywhere below).

## Domain and conventions (per prefire F2, unchanged)

Ω_δ(t) = R³ minus closed tubular neighborhoods (radius δ ~ a) of filament
cores; (α,β,φ) smooth on Ω_δ; vorticity distributional on filaments (Γ_n
circulations); Lin constraint weak (test variations compactly supported in
Ω_δ); core-adjacent variations excluded (named). All conclusions are
FORM/topological statements, δ-independent by construction (prefire δ-gate:
any quantitative bridge would repeat at δ/2 — none is claimed here).

Native action (standard Lin/Clebsch form for incompressible Euler, sign
convention fixed here once):

  S = ∫dt ∫_Ω [ ½|u|² + α(∂_tβ + u·∇β) − p(∇·u) ],
  u transverse after pressure elimination; variation w.r.t. u gives
  u = ∇p − α∇β (Clebsch ansatz with φ ≡ p, cipher's ∇φ + α∇β up to sign).

No imports: no EM constants, no Maxwell laws, no comparator enters any
derivation. Comparisons to Maxwell/BF structures are FORM statements only.

## Verdict vocabulary (pre-registered)

- ESTABLISHED: identity stated as an equality with exact domain + quantifiers,
  derivation in-paper, machine receipt green, must-FAIL mutation verified to
  fail.
- NOT-DERIVABLE: absence claim inside the bare action's identity algebra,
  carried by the ledger's completeness argument (below), never by fatigue.
- OUT-OF-SCOPE: named, priced elsewhere (Kelvin-gauge coupling: 0139;
  compactness: SYN hole #1; dynamical coupling construction: missing-5).
Paper-level kill condition: any identity whose mutation does not fail, or
any ESTABLISHED that reduces to an inequality/fit, is withdrawn by the
author before firewall review.

## Identity ledger to establish (targets; derivations may sharpen them)

- I-CS: w ≡ ∇α×∇β is an identically conserved current (∇·w = 0); ω = w
  exactly (curl-grad); u·ω = det J(φ,α,β) exactly (α∇β term drops).
- I-Helicity: dH/dt = −∮_∂ [P ω·n + ((ω×u)×u)·n] dS with P = p + |u|²/2;
  compact support / flow-transported tube walls ⇒ dH/dt = 0 exactly.
- I-Noether2: the gauged action has vanishing second functional derivative
  w.r.t. every ∂_iA_j: Hessian tensor ≡ 0 ⇒ A-equation is order-0 ⇒ A is
  auxiliary (constraint) — independent route to bridge Reading 1; any
  propagating completion must ADD a (∂A)-kinetic term (named, not imported).
- I-Sing: on Ω_δ with frozen-in transport (Kelvin: Γ_n const), H splits as
  H_bulk + ΣΓ²Slk + 2ΣLkΓΓ with every term an exact dynamical invariant:
  d/dt(Lk_ij) = 0, d/dt(Slk_n) = 0 for closed transport without reconnection
  (reconnection outside the F2 representation — prefire boundary).
- I-Decouple (the classification): the conserved-current algebra of the bare
  action (energy, linear/angular momentum, relabeling/Casimir charges) is
  DECOUPLED from the topological charge: no identity couples Lk to
  translation/relative motion. Charge-motion coupling is NOT-DERIVABLE from
  the bare action; every such coupling requires extending the action or the
  constraint set (the exact content of missing-5).

## Firewall pre-audit (drift checklist, self-applied)

- Input/output: imported = Clebsch/Lin formalism (standard mathematics),
  CWF/linking-number topology (cited), Kelvin's theorem (classical, derived
  in-paper as receipt-checkable algebra). Derived = everything in the ledger.
  No fitted constant appears; nothing numeric is called exact.
- Forbidden transfers pinned: no axisymmetry assumed (all algebra general);
  I-Helicity is a dynamical identity (no steady-map transfer); Ω_δ is
  Euclidean, no periodic-box transfer; no L1-diagnostic or closure-as-range
  steps exist (identities are equalities, not range claims).
- Load-bearing steps: I-Decouple's completeness = enumeration of the bare
  action's symmetry group (translations, rotations, Galilean boosts,
  relabeling SDiff of the Lagrangian labels, β-shift gauge) — the
  enumeration is the claim and is itself auditable; it licenses absence
  INSIDE the bare action only, never inside extensions.
- Oracles: each identity ships one must-FAIL mutation in the receipt; the
  mutation for I-Noether2 is exactly the "added (∂A)-kinetic term" that any
  propagating completion must introduce — the certificate detects the
  required extension by design.

## Non-goals (OUT-OF-SCOPE, stated to prevent scope creep)

No Kelvin-gauge coupling construction; no compactness source; no dynamical
charge-observable build; no reconnection model; no electron identification;
no promotion; no claims beyond the ledger's scope.
