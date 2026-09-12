# 00 — HJ2 augmented-map parameter-row Schur core (Route A, Obligation A)

Attempt: P253/0159 (sage; 0054 contract continuation — Route A first
step). Owner directive charter: resolve HJ2 at its stated scope (one
polynomial-profile Cao subfamily, routes HJ2-A/B/C); either the upgrade
or a named-mechanism kill banks. THIS ROUND BANKS: the parameter-row
Schur core of missing-construction 1 (the exact augmented Cao carrier
map) is DISCHARGED at the even/polynomial subfamily scope. HJ2 itself
REMAINS ACTIVE: constructions 2–4 are untouched. Receipts:
receipts/run_hj2a.py + run_hj2a.log (10 assertions — 8 identity + 2
mutations — exit 0; vacuity-checked per standing rule; one real bug
caught pre-commit: the draft's own parity table mislabeled d1U and
carried an `or True` clause — both removed before landing).

## Setting (all conventions = 0054 derivation.md)

Fixed-core equation (6): −Δw + q/(1+qy₁)·∂₁w = (1+qy₁)²(w₊)^p; limiting
cells (7)–(8): L_U V = 2y₁U₊^p − ∂₁U, with L_U = −Δ − pU₊^{p−1};
kernel ker L_U = span{∂₁U, ∂₂U} (source Lemma 3.8); even symmetry in y₂;
center gauge ∂₁w(0) = 0; circulation normalization −2πU′(1) = M₀ :=
∫U^{p+1} (15). Frozen range: integer p ≥ 6.

## The result

**Theorem (subfamily scope).** For the even, center-gauged
polynomial-profile Cao subfamily at p ≥ 6, the augmented map's
parameter-row Schur determinant is

    D_red = −M₀·(p+5)/(2(p+1)) ≠ 0,

and the constrained linearization is invertible; hence Lyapunov–Schmidt
produces an actual constrained steady branch w(q), τ = q²log q, with the
candidate jet (16) as its expansion — the formal polyhomogeneous
candidate of 0054 §2.2 becomes a CONSTRUCTED branch at this scope.

Receipt chain (each step an exact identity, all receipted):
- **HJA-1 virial:** ∫|∇U|² = ∫U^{p+1} = M₀ (IBP; boundary term zero
  because U vanishes on the free boundary).
- **HJA-2 isotropy:** ∫(∂₁U)² = M₀/2.
- **HJA-3 q-row:** ⟨∂₁U, L_UV⟩ = 2∫y₁U^p∂₁U − ∫(∂₁U)²
  = −2M₀/(p+1) − M₀/2 = **−M₀(p+5)/(2(p+1))** (y₁-IBP through
  ∂₁(U^{p+1}) = (p+1)U^p∂₁U; boundary vanishes on the support edge).
- **HJA-4 parity triangularity:** the τ-column and all cross terms
  vanish — every τ/cross integrand is y₂-ODD (∂₂U odd; U^p, y₁, V, Z,
  Z_log sources all y₂-even) against the y₂-even domain/data. The
  parity table is receipted entry-wise.
- **HJA-5 gauge nondegeneracy:** radial ODE at the origin gives
  ∂₁₁U(0) = −U₀^p/2 < 0 ≠ 0 — the center-gauge row kills the ∂₁U
  kernel direction without residue.
- **HJA-6 constrained invertibility:** L_U = −Δ − (bounded, compactly
  supported potential) is a compact perturbation of the Dirichlet
  Laplacian ⇒ Fredholm index 0; self-adjoint ⇒ coker = ker;
  ker = span{∂₁U, ∂₂U} with ∂₂U removed by evenness and ∂₁U removed by
  the gauge row (HJA-5) ⇒ constrained kernel and cokernel are {0}.
- **HJA-7 reduction:** by HJA-4 the Schur matrix is triangular with a
  zero τ-column; the Cao restriction τ = q²log q removes the degenerate
  column rather than requiring it — the REDUCED determinant is D_red =
  M₁₁ ≠ 0 on the frozen range (rational-grid verified p ∈ {6,7,8,9,12};
  strictly negative for all p > −1).
- **Mutations (both bite):** MB-HJA-1 — pairing the q-row with ∂₂U
  collapses D to 0 (the pairing is load-bearing); MB-HJA-2 — assuming V
  y₂-odd breaks the vanishing table (parity['V'] = even is load-bearing).

## Honest scope and what this does NOT license

- **HJ2 is NOT discharged.** The graph-theorem heart — common DA graph
  domain (construction 2), two-index symbol bounds (3), all-sector
  X*→D* resolvent (4) — is untouched. 0052's full-sector Grushin and
  Riesz contour remain gated by those three.
- Construction 1 is discharged **at this scope only**: the Schur core +
  constrained invertibility. The full C^{3,α} augmented map still owes
  the standard elliptic regularity assembly (source Lemma 3.8/3.9-class
  a-priori estimates + Theorem 1.6/C.2 applicability), stated as
  remaining assembly, not new mechanism.
- **Parity is load-bearing and scope-frozen:** the triangularity dies
  outside the even subfamily (MB-HJA-2). The odd/even two-mode
  continuation is a different, unbuilt route.
- No Riesz transfer claim, no nonlinear branch claim, no stability/
  quantization/electron/neutrino claim. 0054's does_not_license list
  travels except that "discharge of HJ2" is now PARTIAL: construction 1
  core done, constructions 2–4 remain.

## Route ledger update (0054 verdict vocabulary)

- HJ2-A Obligation A (chart jet (12) construction): the augmented-map
  parameter row + constrained invertibility DISCHARGED at subfamily
  scope (this round); uniform inverse/Jacobian BOUNDS remain named.
- HJ2-A Obligations B/C, Route B, Route C: UNCHANGED (blocked with the
  named constructions of 0054 result.yaml).
- Continuation ladder: next decisive unit = construction 2 (common DA
  graph domain) or Route B form-side construction; neither is started
  here.
