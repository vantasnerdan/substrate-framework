# 0137-atlas-xdomain — cross-domain survey: native charge-like coupling

Charter: shepherd IDEATION (closures-only → generate). Scope: how other
programs derive charge-like coupling from neutral substrates; concrete
mechanisms with importable mathematics + bridge sketch to Euler each. NO
imports proposed — native math only. Confidence: standard-physics recall
(unverified at source here — verify before building); campaign walls cited
from attempts/ reviews. Brief per physics-discovery: tension = neutral Euler
must yield quantized selection + back-reaction, but Magnus-filament closed
(0134/0135), punctures divergent (0065), imports declined (R-EM2).

## X1 — Symmetry reduction curvature (Marsden–Weinstein; STRONGEST)
- Mechanism: reducing a symmetric Hamiltonian system by a symmetry group
  leaves curvature (magnetic) terms in the reduced Poisson bracket —
  gauge-like forces from pure geometry, no charge postulated. Cotangent-lifted
  action → reduced space with curvature 2-form acting as effective magnetism
  on slow variables.
- Importable math: reduction theorem + reduced bracket with curvature term;
  Hannay–Berry connection for cyclic slow variables (classical Hannay angles
  in fluids literature).
- Bridge sketch: fast filament oscillations (A3 monodromy modes, banked)
  as fibers over slow vortex-position base; adiabatic circuit of the slow
  variables picks up Hannay holonomy = effective Lorentz-like deflection.
  Compute holonomy 1-form from banked A3 mode shapes; test deflection sign
  vs S1-3D wrong-way record (discriminating check).
- Wall-check: needs clean scale separation (fast/slow) — 0124 floor (~1e-2)
  may smear adiabaticity; needs the orbit chart 0071 named missing. But uses
  ONLY banked modes + reduction machinery. Rank 1: native, quantitative,
  falsifiable by sign.

## X2 — BF observables as linking charge (already S3-residue; SHARPEN)
- Mechanism: in BF theory, charge IS linking number — Wilson lines/surfaces
  measure topology; no local source needed, hence no puncture divergence.
- Importable math: abelian BF action `∫B∧dA`; observables = linking of
  surfaces/lines; first-order (no F²) — matches S3 residue exactly.
- Bridge sketch: Euler helicity + vortex-tube linking as the charge
  variable (integer by topology). DYNAMICS AS CONSERVATION TEST (drift
  reword — frozen-in vorticity forbids reconnection jumps in ideal Euler):
  linking-class conservation becomes the PREDICTION (charge conserved
  exactly); any charge-changing event needs non-ideal physics, which is
  itself the discriminating test. Next check: measure linking-class drift
  in banked member states (0117 npz) — zero drift CONFIRMS, nonzero drift
  bounds non-ideality.
- Wall-check: S3 bridge died as Maxwell-replacement — this sketch does NOT
  revive S3 as dynamics; BF is kinematics + conservation law. Rank 2:
  kinematics native-exact, prediction sharp.
## X3 — Ponderomotive separation (plasma; HONEST NEGATIVE with residue)
- Mechanism: neutral high-frequency wave fields separate charges by
  mass/mobility (ponderomotive force); effective separation from neutral drive.
- Bridge sketch: Euler acoustic/strain field as the pump; tag-density
  contrast as the separated quantity. KILLED as charge mechanism: single
  Euler fluid has no second species — nothing to separate. This is the
  S2-acoustic kill from another angle (type error, cf. 0130).
- Residue: two-timing/averaging machinery itself is reusable for X1's
  adiabatic step. Rank: dead as mechanism, alive as method.

## X4 — Gauged-Ginzburg–Landau phase stiffness (CONDITIONAL)
- Mechanism: GL vortex lines carry quantized circulation; gauging the U(1)
  phase yields flux quantization + Maxwell-like equations from phase
  stiffness. Charge = winding number; robust by topology.
- Importable math: GL free energy, winding quantization, London limit.
- Bridge sketch: needs a complex order parameter Euler lacks natively.
  ONLY native candidate: micropolar orientation field (C-CST-018's sector)
  as the phase — orientation defects as charge carriers with winding
  selection. Check: does the accepted micropolar action admit stable
  line defects with quantized winding? Pure analysis, no compute.
- Wall-check: S3 warns against auxiliary fields (auxiliary-A death);
  orientation must be ACCEPTED-sector (C-CST-018 is), not invented.
  Rank 3: one assumption (defect stability), high payoff (quantization).

## X5 — Filament Hasimoto → NLS U(1) (RECOMBINATION)
- Mechanism: Hasimoto map takes filament curvature/torsion to NLS; NLS
  power (L² norm) is a conserved U(1) charge-like quantity with its own
  continuity equation — charge-like conservation from pure geometry.
- Importable math: Hasimoto transform, NLS conserved densities.
- Bridge sketch: banked A3 monodromy modes already live near filament
  dynamics; map the m1–6 signatures to NLS-side densities; ask whether
  the NLS power current couples back to filament motion (self-induced
  velocity correction) = native back-reaction in NLS variables.
- Wall-check: self-stretching may break the map's validity (vortex-stretching
  vs NLS integrability); the 0124 floor may dominate the correction.
  Rank 4: elegant, needs validity audit first.

## X6 — Kinetic-moment charge separation (EXCLUDED, recorded)
- Mechanism: Vlasov→fluid moment hierarchy can carry charge separation via
  pressure-tensor anisotropy from neutral kinetics.
- Exclusion: needs a kinetic substrate beneath Euler — violates
  Euler-native scope (new physics, not new math). Recorded so the lane
  shows it was considered, not missed.

## Ranking + next checks
1. X1 reduction curvature — compute Hannay 1-form from A3 modes; sign test
   vs S1-3D record. Needs 0071 orbit chart (budget it).
2. X2 BF-linking kinematics — measure linking-class drift in 0117 npz;
   reconnection dynamics = missing-5 neighborhood.
3. X4 micropolar defects — winding-stability analysis in C-CST-018 sector.
4. X5 Hasimoto/NLS — validity audit (stretching vs integrability).
- Dead as mechanisms, banked as methods: X3 averaging → X1; X6 excluded.
- All sketches native-math-only; no imports proposed; verify recalled
  physics at source before any build.
