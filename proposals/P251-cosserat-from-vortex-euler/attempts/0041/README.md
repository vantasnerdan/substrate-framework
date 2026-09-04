# Attempt 0041: full Biot--Savart finite-core twist gradient

## Frozen delegated contract

Parent objective: the exact smooth-Euler micropolar continuum with the
original slowly varying affine coarse-graining premises, finite cores,
translational action, and stationary EPS compatibility. This child owns only
this new directory. Attempt 0036 is frozen for independent review; no shared
API or earlier attempt is changed here.

Positive deliverable: derive exact finite-core integral expressions for the
two-triangle axial twist-gradient matrix from the complete three-dimensional
Biot--Savart isovortical energy and prove its positive leading logarithmic
part with a controlled bounded remainder. Include the second-order vorticity
variation, mutual-core terms, and an explicitly typed infrared regularization.
The related common-bend sector will be derived if the same representation
supplies it. An unrestricted all-wavenumber equality is not part of the
original slow-varying affine objective and is not added here.

Inputs are constant-density incompressible Euler, the smooth separated
six-core transverse profile from 0036, the exact volume-preserving twist
pushforward, and the three-dimensional Newton kernel. No local-induction
coefficient, empirical value, or fitted target is an input. Candidate
representations are axial Fourier/K0 kernels and a rotation-invariant fixed
outer separation window; select by exact second-variation closure and a
uniform thin-core remainder. The strongest oracle is independent symbolic
variation plus analytic kernel estimates. No numerical spectrum, soft
eigenvalue, or energy-difference computation is yet designed.

## Active record

The complete Biot--Savart variation gives an exact positive projector Gram
matrix for the two triangle twists. Its finite local coefficient has a
rotation-invariant real-space integral and a controlled leading logarithm;
all mutual terms and the second vorticity variation are retained. Common
bend is derived from the same representation, with its distinct physical
infrared logarithm explicitly stated.

The full energy also exposed a necessary companion: the conjugate radius
imbalance has gradient energy. Its elimination modifies optical dispersion
at the same derivative order. The exact finite-core matrix rule and the
affine-cage collective-field correction are derived, rather than applying
the free-angle factor to the constrained field.

Artifacts:

- `full-biot-savart-twist.md`: exact functional, second variation, positive
  projector, finite-core integral, uniform leading-log bound, common bend.
- `radial-momentum-gradient.md`: conjugate-radius gradient, momentum
  elimination, exact finite-core affine-cage matrix rule, and separately
  the derived thin-core asymptotic coefficient.
- `full_twist_energy.py`: 21/21 first-run exact checks, captured in
  `stdout.txt` and `stderr.txt`.
- `radial_gradient_elimination.py`: 8/8 first-run exact checks, captured
  in `radial-stdout.txt` and `radial-stderr.txt`.
- `cage_gradient_matrix.py`: 7/7 first-run exact checks, captured in
  `cage-stdout.txt` and `cage-stderr.txt`.

All three scripts pass Ruff. These are exact algebraic receipts; the analytic
Gram positivity and logarithmic remainder are proved in the documents.
An exposition audit distinguished the omitted second-variation term (wrong
static common-rotation energy) from an omitted longitudinal resolvent change
(doubled twist logarithm). The last mutation-check label was corrected to
name the latter mechanism; the original first output remains preserved and
`replay-stdout.txt`/`replay-stderr.txt` capture the corrected-label replay.
No numerical energy difference, soft-eigenvalue estimate, or fitted cutoff
was used. The twist local coefficient is infrared finite by the triangle
centroid constraints. The parent stationary-EPS and remaining continuum
constructions remain active and are not declared complete here.
