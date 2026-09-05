# 0057 — one-action finite-coherence material/orbit assembly

Parent P251 / issue #198; owner `/root/orientation_construction`; this
directory only. Base release v0.171.0. Parent objective unchanged.

Frozen child obligation: derive a single declared finite-coherence affine
ensemble from the material Euler action whose centroid kinetic sector and
centered EPS orbit internal sector are counted exactly once. Retain all
pressure/shared-face constraints, exterior kinetic motion and impulse,
and the physical-to-canonical current map. Coordinate with 0055's actual
compact coordinate constructor, including cages outside the original
invariant tube. No whole-space boost is assigned a finite energy.

Candidates: (A) a full finite material partition, with exact centroid
decomposition before internal reduction and explicit ambient bookkeeping;
(B) a finite material observation region plus a relative exterior orbit,
with its interface action and impulse retained. Selection criteria are
one Euler action, no duplicated mean energy, genuine independently
testable slow-affine kinematics, exact physical/current observables, and
no assumed target constitutive energy. A local RVE ensemble is explicitly
distinguished from an arbitrary global gluing theorem.

Oracle: exact material kinetic and symplectic decompositions, boundary
variation and partition identities, plus finite-dimensional symbolic
checks. No numerical solver, spectrum, empirical comparator or fitted mass.
Status: frozen before substantive assembly work.

Append-only representation expansion: a single decaying EPS realization
has a finite defect action but does not itself provide a normalized
positive-density fluid ensemble. Candidate C uses the source-constructed
stationary Gaussian Beltrami measure of Enciso--Peralta-Salas--Romaniega,
arXiv:2006.15033. It has finite pointwise variance, full local Beltrami
support, and a positive density of knotted invariant structures. Reconstruct
the compact EPS core/cage coordinates on positive-probability good patches
and derive the material/orbit action per unit TOTAL fluid volume. This is
a genuine new representation; the old decaying-field coefficients are not
transferred unchanged. The parent approved and activated this route.

The new common-rotation density and physical-angle reconstruction are being
checked independently by `/root/construction_review` in 0059. This child
owns the normalized material ensemble, compact moment reconstruction,
stationary full Hessian, source audit and ambient bookkeeping.

Selected continuation: 0059 removes the unrequested independent global
common rotor and uses the physical compact relative angle
`q=Phi-curl U/2`, projected off all eleven incompressible-affine moments.
This preserves the original U,Phi target and retains its complete
gradient kinetic cross. The earlier six-response route remains explicit
branch evidence, not an assumed source of body inertia.

Completed records:

- `stationary-assembly.md`: normalized one-fluid material action, genuine
  good-patch ensemble, moment construction, full stationary Leray/Hessian
  density, operator-level reaction and all-fluid current bookkeeping.
- `affine-spectral-energy.md`: exact same-field affine energy,
  `mu_affine=2 rho U_*²/5`, the full negative reaction Schur correction and
  finite positive bound, plus distinct coherence/background-removal limits.
- `slow-locality.md`: compact smooth double-divergence construction from
  the eleven moments and exact C² slow jets retaining all Leray tails,
  without an added mixing or infrared regularization premise.
- `source-audit.md`: checked primary theorem numbers and SHA receipts;
  independent source review is in 0060.

`verify.py` passes 26/26 exact checks and calls the canonical
`hermitian_schur_jet` API for the complete noncommuting reduction. The first
20-check execution remains in `stdout.txt`; the 23-check canonical-API
extension remains in `replay-stdout.txt`; the final analytic cross/limit
extension is in `final-stdout.txt`, each with its matching stderr file.
Ruff and scoped diff checks pass. No numerical eigenvalue or fitted
constitutive coefficient is used.

Route verdicts: direct isolated-defect normalization is blocked by its
missing normalized positive-density volume measure and is superseded here;
the stationary material/compact
assembly and exact affine/slow-jet constructions are established as stated,
with individual construction review pending. This child does not declare
the parent's kinetic normal form, whole claim promotion, or PR complete.
