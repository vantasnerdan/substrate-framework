# S3 Dávila–del Pino–Musso–Wei verified (beacon 0110)

Source: *Leapfrogging vortex rings for the 3-dimensional incompressible Euler
equations*, arXiv:2207.03263v4 [math.AP] 11 Nov 2023. Fetched abs + full PDF
this session (V2; 3089-line extraction, artifact://35).

## Exact consumed scope

- Theorem 1: for every collisionless solution `q(τ)` of the reduced Hamiltonian
  system (1.14) on `[0,T]` (`inf|q_i−q_j| > 0`), there is a SMOOTH exact Euler
  solution `W_ε` with `k ≥ 2` rings of core size `ε`, centers at mutual distance
  `|log ε|^−1/2`, tracking the reduced law with stated `O(...)` remainder (1.17).
- Reduced law (1.14)–(1.15): Hamiltonian point-filament system, periodic orbits
  on closed level curves (k=2 antisymmetric case).
- Building blocks: Kaufmann–Scully profile `U(y) = 8/(1+|y|²)²` + non-radial
  corrections (anisotropy NOT a perturbation); inner-outer gluing; centers
  `P_j(τ) = (r_0,0) + |log ε|^−1/2 q_j(τ) + corrections`.
- Hypotheses: axisymmetric NO-SWIRL class; all rings SIMILAR POSITIVE
  circulation (dipole/opposite-sign excluded — singularity-scenario related);
  finite window `[0,T]` in scaled time `τ = |log ε|·t`.
- Source's own limits: phenomenon NOT expected to persist as `t → ∞`
  (numerics: few crossings, then mixing); k=2 antisymmetric case NOT reducible
  to travelling/rotating waves; time-periodic case called "very difficult".

## Perturbation class / norm

Construction (ε-families tracking a reduced orbit), not a stability theorem:
no perturbation neighborhood, no norm, no return. Core-shape control is the
method's point (vs Marchioro-type results whose supports may expand).

## S2-carrier transfer verdict

- Compatible regime: O(1) ring radius, thin core, no-swirl — same family as
  the S2 Cao carrier; S4 explicitly adopts this regime (§1.2 of S4 paper).
- What transfers NOW: nothing as a theorem. The S2 carrier is a specific
  variational steady family; S3 builds its own ε-family from compact-polynomial
  type blocks. Join needs: (a) S3-type construction with Cao profiles, or
  (b) proof the S2 member sits in an S3 ε-family — both open.
- Preserved distinction (frozen-issue requirement): exact finite-window
  solutions vs reduced law (1.14) — the law is an input selector here, not a
  derived interaction observable of S2.
