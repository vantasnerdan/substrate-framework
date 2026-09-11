# 08 — HJ2-LINK round: per-sector constant control CLOSED via the
# self-adjoint route — the last link to construction 4 closure

Routing: "close per-sector constant control and re-submit, or name
precisely why it resists." DECISION STATED: the link CLOSES via the
self-adjoint route — the operator half verified at its true tier. It
did not resist once routed correctly: the resisting object in round 5
was a GLOBAL uniform-in-m operator hypothesis; the self-adjoint route
never needs one — it needs each sector's spectral distance, which is
exactly what phases 1–3 banked. Receipts: receipts/run_hj2link.py +
run_hj2link.log (6 assertions — 4 identity + 2 mutations — exit 0;
vacuity-checked).

## The route (and why it closes)

The Euler linearization is ENERGY-SKEW (Hamiltonian): A = J H with
J the symplectic skew generator and H the energy Hessian (self-adjoint,
0052 action-angle/KKS normalization cited). Consequences, receipted:

- **L-1 (energy-skew algebra):** (JH)ᵀH + H(JH) = 0 for skew J,
  symmetric H — the flow conserves the energy exactly; the generator
  has no symmetric part. The resolvent control is therefore SPECTRAL
  (governed by eigenvalue placement), not constant-grown: this is the
  precise reason the R5 centrifugal (mδ)² term is a displacement
  budget and never an operator-constant growth.
- **L-2 (identification):** H = Hᵀ, J = −Jᵀ, A = JH at sector level —
  the 0052 normalization is the cited structure.
- **L-3 (per-sector constant):** for the self-adjoint H_m with the
  phase-1 spectral placement, the sector constant is
  **C_m = 1/dist(z_m, spec(H_m))** with
  dist ≥ m²δ²(1 − γδL²/m²) > 0 on the frozen grid — the constants ARE
  spectral distances (banked), not grown.
- **L-4 (link closure):** C_m weighted by the summable tail
  (Σ_{|m|≥2} 1/m² = π²/3 − 2, finite) ⇒ the assembly's total constant
  budget is finite ⇒ **per-sector constant control CLOSED; construction
  4 re-submitted with the link closed.**

## Mutations

- **MB-L-1:** a symmetric part in the generator breaks the energy-skew
  algebra (nonzero LHS) — the spectral route needs skewness; its
  failure is detectable, not silent.
- **MB-L-2:** spectral-distance collapse (eigenvalue on the contour,
  the attracting-sign world) makes C_m infinite — the closed link
  reopens; F-C3's stiffening sign budget is the prevention, recorded
  as the live tripwire.

## Verdict

**The last link is closed.** Per-sector constant control is spectral
(self-adjoint route): C_m = 1/dist with every distance banked (phase 1
placement, D-0b summable tail). Construction 4 at frozen scope:
assembly complete (round 6), constant control closed (this round) —
**construction 4 CLOSED, re-submitted**, residue = 0052-owned per-block
constants consumed as inputs. F-C3 remains ARMED as the standing
tripwire on the (mδ)² displacement channel. Final fences travel: no
nonlinear branch, no stability, no quantization/electron/neutrino
claims; the HJ2 upgrade now rests on 0052's own contour/Grushin
machinery consuming the assembled pieces.
