# Validation and failure-derived oracle receipt

## Boundary

Only `P253/0019`'s bounded reviewer-requested normalization wording and
`P253/0025` were edited.  No source module, central registry, proposal
manifest, governance file, memory, skill, other attempt, or commit was
changed.  Central activation existed with exit zero before 0025 source-body
work.  The new full paper remains in `/tmp/primary-source-cache/P253-0025`;
its URL, version, SHA-256, and exact theorem/equation locations are in
`access-inventory.md`.

## Analytic oracle

The strongest oracle is the agreement of four exact layers:

1. Baldi's source action chart yields the deformation and cotangent shear, and
   its published coefficients prove that the derivative-frequency ratio is
   nonconstant;
2. direct Kelvin--Leray calculus retains the coefficient-two global pressure
   row, preserves transversality and polarization triple area, and derives
   `c_dot=Lc+(vorticity dot k)A`;
3. the same exact amplitude system supplies the bounded pressure-normal basis
   and the periodic resonant `SL(2,R)` monodromy; and
4. a fixed-time profile argument with normalization by the exact Biot--Savart
   norm transfers any genuinely unbounded ray-gain sequence under the
   quantifiers `forall j exists delta_j exists N_j`.

The result deliberately does not claim a sign for the resonant monodromy
discriminant.  Closeness to a nontrivial Jordan matrix fixes a transient but
does not fix whether the exact repeated return is elliptic, parabolic, or
hyperbolic.

## Exposing verifier history

`verify_resonant_skew_product.py` is an exact symbolic verifier; no sampled
multiplier, eigenvalue tolerance, or soft numerical quantity occurs.  Its
first two runs failed on SymPy structural equality rather than unequal
expressions.  Those receipts are preserved.  The third run exposed a genuine
analytic error in the attempted route: `c=k cross A` is not a Cauchy vector
when `vorticity dot k` is nonzero.  That receipt is also preserved, the claimed
growth route was removed, and the verifier was strengthened to assert both
the corrected identity and nonzero shortcut residual.

The first execution of the repaired predicate set is preserved in
`first-passing-run.command.txt`, `first-passing-run.stdout.txt`, and
`first-passing-run.exit`.  A final predicate was then added to expose the
claimed bounded complement directly; its first execution is preserved in the
three `strengthened-run` receipts.  The final verifier checks exactly:

- the source-series derivative ratio and its leading decreasing slope;
- the action-angle cotangent shear and derivative resonance;
- canonical symplecticity of the six-dimensional cotangent lift;
- transversality preservation by the full pressure coefficient and failure
  when that row is omitted;
- the exact vorticity correction to the false Cauchy shortcut;
- conservation of the polarization triple area;
- the exact neutral pressure-gradient amplitude `A=u`; and
- exact evolution of its normalized complement when
  `grad|u|^2` is pressure-normal; and
- both terms in the coadjoint-accessibility triple product.

No small-ratio numerical observable exists, so the small-ratio-numerics skill
does not bind this execution.

## Quantifier and global-domain audit

The PDE-transfer theorem does not reuse a fixed-time error at an increasing
time without control.  For each already fixed finite circuit count it first
shrinks the action band, then chooses frequency.  Constants may depend
arbitrarily on the circuit count and band width.  This proves the conditional
operator-norm implication but not a single uniform-in-growing-time expansion.

The profile comparison retains the whole-space Leray multiplier, order-minus-
one Biot--Savart normalization, adjacent commutators, separated smoothing
kernels, flat support boundary, noncompact velocity tail, and Euclidean
symmetry quotient.  It introduces no annulus wall.

## Claim boundary

The attempt constructs the exact repeated returned-fibre object and closes
the continuum-transfer quantifiers, but it does not decide the object's
discriminant.  Its strongest unconditional physical statement is the bounded
pressure-normal principal sector.  The remaining dependency is a concrete
ordered-integral sign/enclosure for the resonant monodromy, followed by nearby
profile control if elliptic.  This is a route gap, not a physical-impossibility
or quantum-mechanics claim.
