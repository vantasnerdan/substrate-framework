# Validation and claim boundary

The analytic proof in `derivation.md` supplies the outgoing boundary-value
sign, exact power prefactor, full shell gradient, coarea weight, shifted-
anisotropic Helmholtz LAP reduction, switched continuity identity, late free-
field spectral norm, the full-real pushforward/autocorrelation plateau bound,
and the conditional reciprocal-amplitude conversion.

`verify_shell_flux.py` independently rederives nine algebraic consequences:

1. the negative source-work delta term and positive phasor-half power;
2. the star-shaped positive-frequency shell;
3. `|partial_r D|=2*c_EM*|omega|`;
4. the exact sphere coarea weight;
5. the full-gradient subluminal lower bound;
6. cancellation of the envelope continuity defect by `a_prime*K`;
7. the directly integrated Gaussian curl-current power;
8. the factor two between amplitude width and fractional energy decay; and
9. the negative-frequency/negative-speed conjugate-shell relation without
   double counting.

The exact receipt exits `0` with nine `PASS` lines and empty stderr.  The
focused public-API receipt reports six passing tests and empty stderr.  Tests
also expose zero frequency, luminal speed, superluminal speed, nonpositive
permittivity, and `|n_z|>1` inputs.  These are exact symbolic/regression checks,
not production radiation numerics.

The oracle does not prove the weighted outgoing resolvent or time-domain
limit by execution.  Those are analytic arguments with explicit hypotheses in
Sections 3--4.  The `T*P+O(1)` result requires the scalar resonant first
correlation moment, endpoint cross-correlation norms, finite endpoint energy,
the conjugate cross first moment, and fixed smooth ramps stated in
(23)--(23h).  The exact late spectral norm remains valid without those norms.
The sharp switch remains a weak-source limit with endpoint energy to be
computed in its chosen regularization.  Equation (25) proves only the finite
moving-volume identity; the infinite-world-tube/late-field equivalence stays
open.  The one-block width formula remains conditional on a reciprocal
amplitude equation and a nonzero physical current trace.

No full-suite run was launched for this additive conditional API.  The earlier
2784-pass suite began before the 0087 and 0090 files existed and is not cited as
evidence for either attempt.
