# Primary-source transfer for the one-sided trace route

The sources below were opened only after the central 0263 schema receipt
passed.  Each is used at its exact scope; none is cited as a substitute for
the new barrier calculation or the invariant-extension argument.

## Degenerate Feynman--Kac representation

[Feehan, Gong and Song, *Feynman--Kac Formulas for Solutions to Degenerate
Elliptic and Parabolic Boundary-Value and Obstacle Problems with Dirichlet
Boundary Conditions*, arXiv:1509.03864](https://arxiv.org/abs/1509.03864)
treat Markov generators whose symbol degenerates as the `2 alpha` power of
distance, including the quadratically degenerate case `alpha=1`.  Their
elliptic results distinguish natural/entrance boundaries from regular/exit
boundaries and give uniqueness of the stochastic representation with a
partial Dirichlet condition on the reachable boundary.

This licenses the representation logic used in 0263 after the actual SDE,
normal drift, boundary classification, exit moment, and killing sign are
checked.  Those checks are equations (3)--(9) of `trace-estimate.md`.  The
paper does not supply the sharp `|m|^{-2/3}` trace multiplier, a smooth exit
density uniform at this characteristic edge, or Hanzawa-parameter tame
bounds.

## Maximum principle for sums of squares plus drift

[Bony, *Annales de l'Institut Fourier* 19 (1969), 277--304,
DOI 10.5802/aif.319](https://doi.org/10.5802/aif.319) studies operators
`sum X_j^2+Y+c`, proving strong maximum-principle and related uniqueness
results under the corresponding Lie-algebra accessibility hypotheses.  The
0261 fields and first drift commutators span at the edge.  Bony therefore
supports the comparison step once the sign of the killing term and the
inward accessible set are fixed.

It does not select the physical trace, quantify its angular regularity, or
give a moving-boundary Poisson operator.  In the radial proof those jobs are
done by the exact entrance diffusion and the explicit barriers (14) and (18).

## Hypoelliptic smoothness and derivative bounds

[Kusuoka and Stroock, *Applications of the Malliavin Calculus III*, J. Fac.
Sci. Univ. Tokyo Sect. IA Math. 34 (1987), 391--442,
DOI 10.15083/00039484](https://repository.dl.itc.u-tokyo.ac.jp/records/39493)
develop smooth-density and derivative estimates for diffusions satisfying
uniform Hörmander conditions.  Together with the already transferred
Hörmander and Bramanti--Zhu local estimates in 0261, this supports interior
hypoelliptic smoothing for the smoothly extended coefficient family.

It does not directly estimate a killed exit kernel at a characteristic
boundary.  Attempt 0263 does not make that unsupported transfer.  Instead it
extends the **diffusion and its stopped resolvent** across the invariant edge,
then applies interior hypoellipticity to the constructed two-sided solution.
The collar shrink keeps the estimated region away from the true exit
boundary.  This is why only a low norm of the remote exit datum occurs.

## Net license

The primary literature licenses partial-boundary stochastic uniqueness,
entrance/exit logic, maximum-principle comparison, and interior hypoelliptic
regularity.  The exact radial all-mode estimate and the invariant-extension
bridge are the new analytic work in 0263.  Their composition closes the
nonradial smooth tame entrance estimate without claiming that any cited paper
itself contains the moving-free-boundary theorem.
