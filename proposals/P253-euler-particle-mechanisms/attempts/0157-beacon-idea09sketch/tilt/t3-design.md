# T3 readout-staging mini-freeze — FROZEN PRE-COMPUTE (beacon, 0157/tilt)

Status: FROZEN 2026-09-12. Staged in tilt/design.md T3; pair lane
STOP fired (D-D3 ledger). Freeze commit precedes any compute
commit. Fences inherit (FB-C1 NEVER fires, edge not ALIVE, MA-4,
static tilt only, no measurement/carrier claim).

## Object

Defect signature in the tilt channel: how a probe tilt wave reads
a defect through banked J(eps) + gate-2/3 machinery. Staging only:
a signature is NOT a measurement (FB-C1 stays unfired).

## Frozen structure (predictions, not results)

- T3a (guard-level): defect-induced STATIC tilt is NULL at derived
  order (MA-4 consequence): theta-equation source propto J(eps).th
  vanishes at th = 0 — stiffness modulation, never a source.
  Receipt is analytic (no defect term surviving th=0).
- T3b (physics): probe plane tilt wave + J-heterogeneity -> Born
  scattering, cross-section per length, parametric in (beta, Th, q).
  Pre-registered CONTRAST: screw INVISIBLE at derived order
  (T1a: J(eps_screw) = J0 EXACT — no shear channel) while edge
  VISIBLE (eta pattern fires J). If edge is also null the readout
  channel is dead — lane closes, disclose.
- T3c (ledger): what the tilt trilogy banks (T1 source + T2
  back-reaction + T3 signature) and what a real readout would still
  need (named gaps, no construction).
- F-T3a (guard): any th-independent defect source in the
  theta-equation -> implementation bug (MA-4 travels).
- F-T3b (guard): screw Born amplitude != 0 -> bug (T1a exactness
  travels); nonzero screw scattering is a surprise, hold lane.

- F-T3c (physics): edge cross-section identically 0 -> channel
  DEAD (verdict, not failure).
- STOP: bank T3a/T3b/T3c; then STOP — detector/probe-design and
  carrier physics need own charters. No new formalism: Born +
  banked J(eps)/gate-2/gate-3 + Volterra only.
