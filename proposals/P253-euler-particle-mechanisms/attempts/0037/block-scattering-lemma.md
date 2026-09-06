# The exact abstract block estimate and its Euler entries

The preceding calculations leave an operator estimate, rather than an
unspecified stability argument.  This section records the estimate in the
form needed by the exact Euler projection and separates the entries already
proved from the two entries still requiring construction.

## 1. A two-channel smooth-perturbation lemma

Let `H=H_P direct-sum H_Q`.  Suppose `A_P` and `A_Q` generate uniformly
bounded groups after the translation kernel of `A_P` is modulated.  Let
`B_P:H_P->Y_P` and `B_Q:H_Q->Y_Q` be Kato-supersmooth in both time
directions.  In particular,

    integral_R ||B_a exp(t A_a)x||^2 dt
       <=K_a^2 ||x||^2,  a=P,Q.                         (1)

Let `A_0=A_P direct-sum A_Q` be closed on the product graph domain. Assume
either that `A_0+V` with the following factorization is already known to be a
closed generator, or define its resolvent by the Kato formula (3) and prove
the resolvent identities and growth bounds needed to recover a generator.
Let the off-diagonal perturbation have the factorization

    V_PQ=B_P^* C_PQ B_Q,
    V_QP=B_Q^* C_QP B_P,                               (2)

where adjoints are taken in the positive channel metrics. Diagonal
Kato-smooth remainders may be included by enlarging `B_a` and the matrix
`C`. Both forward and backward inhomogeneous supersmooth estimates are
required; the homogeneous bound (1) alone does not control the Duhamel
series. Under these domain and two-sided estimates, the sandwiched resolvent
identity on either open half-plane is

    B(lambda-A_0-V)^-1 B^*
      =[1-B(lambda-A_0)^-1B^* C]^-1
        B(lambda-A_0)^-1B^*,                           (3)

with `A_0=A_P direct-sum A_Q`.  Here supersmoothness includes the uniform
sandwiched-resolvent bound

    sup_(Re lambda !=0)
      ||B(lambda-A_0)^-1B^*|| <=K^2,                  (4)

up to the fixed Fourier normalization.  Therefore

    K^2 ||C||<1                                       (5)

makes the inverse in (3) a uniform Neumann inverse. The same Neumann series
in the Duhamel equation for `B exp(t(A_0+V))`, using the assumed generator
realization and the two-sided inhomogeneous estimates, gives bounded incoming
and outgoing wave operators and a uniform group norm on the modulated
continuous subspace. Equation (3) then excludes spectrum in either open
half-plane covered by that realization. This is the conditional
resolvent/Duhamel supersmooth-perturbation lemma; no finite-dimensional
spectral truncation occurs.

The point of (5) is its product structure.  A source of size `mu` supported
over axial length `L_mu` has `K^2||C||=O(mu L_mu)`, not `O(mu)` multiplied by
an infinite time.  Since `mu L_mu->0`, it is perturbative on the fast channel.

## 2. The exact Euler choices

On `Q`, use the positive column metric from `0030`.  The strict noncritical
group-speed gap and the exact Plancherel calculation in
`scattering-block-reduction.md` give (1) for a cutoff to the solitary region,
with `K_Q^2<=C L_mu/delta_2`.  The principal transport perturbation is put
into the one-way transport before factorization; its characteristic speed
stays below `-delta_2/2` for small `mu`.  The remaining full-pressure and
coefficient-gradient terms are order zero in the energy graph and have size
`O(mu)`.  Thus their `Q-Q` Birman--Schwinger norm is `O(mu L_mu)`.

On `P`, the limiting generator is

    A_P,0=+sigma partial_y L_*,
    L_*=-partial_yy+1-3 sech^2(y/2).                  (6)

Here `y=X-sigma T` is the right-moving profile coordinate. Linearizing
`A_T+sigma partial_X(A_XX+beta A^2)=0` there gives the displayed plus sign.

After fixing the KdV mass/momentum row and quotienting `A_*'`, the spectrum
and constraint slopes of `L_*` give the constrained energy boundedness
recorded in `threshold-reduction.md`. They do not by themselves prove the
two-sided Kato supersmooth or local-smoothing resolvent estimate for the
nonselfadjoint generator `partial_y L_*`. That scalar resolvent estimate on
the stated constraint/translation domain is an additional conditional
hypothesis of this lemma and an open step before the Euler transfer.

For an actual Euler critical coordinate, use the regular-label density
`delta M_D` from `mixed-casimir-flux.md`, together with its flux companion.
The conservation law puts an output `partial_z` in that density row, but this
pair is not yet the adjoint oscillator spectral pair.  A bounded change of
coordinates must retain both propagation branches and expose the free
off-diagonal entries before any smoothing factor or `O(mu L_mu)` round-trip
bound is claimed.  The opposite `P->Q` row has the proposed localized scale
`O(mu)`; its actual graph estimate remains part of the same construction.

`ray-exclusion.md` proves zero high-frequency bicharacteristic-amplitude
exponent for the one-pass perturbation. That principal-ray statement alone
does not exclude high-frequency point spectrum or a Birman--Schwinger pole
for the full nonlocal Euler operator. The actual Euler factorization,
high-frequency operator exclusion, adjoint-oscillator conversion, and the
proposed `O(mu L_mu)` cancellation therefore remain part of the construction,
along with the fixed low/intermediate graph window and exact `k=0` exterior
row.

## 3. Two remaining estimates, with verdict consequences

The all-time axisymmetric linear theorem follows from (3)--(6) once the
following two estimates are proved for the actual `0034` operator:

1. **Critical graph convergence.**  The bounded regular-label density/flux
   pair is converted to the physical adjoint oscillator pair on a fixed
   low-frequency window.  The resulting projection is uniformly bounded in
   the finite-excess Hodge graph and its diagonal generator converges to (6)
   in the sandwiched-resolvent norm, with the translation and fixed mixed-
   Casimir/momentum rows retained.
2. **Full-pressure factorization.**  After the conservation-law derivative is
   carried through that bounded spectral change, every `P/Q` remainder,
   including the derivative of the Bessel Dirichlet-to-Neumann map, factors
   as in (2) with total Birman--Schwinger norm `o(1)`.  The proposed scale is
   `O(mu L_mu)`, but it is a conclusion to prove rather than a consequence of
   the local conservation law alone.

If both hold, (5) proves bounded axisymmetric linear scattering modulo
translation.  A nonzero limiting defect in estimate 1 becomes the corrected
critical amplitude operator; a nonzero order-one defect in estimate 2 is an
explicit Evans/Birman--Schwinger coupling and must be tested for a pole.  In
either case the calculation has a route verdict.

Even the positive linear conclusion would not supply the requested nonlinear
neighborhood.  Axisymmetric Euler with swirl has no general global-regularity
theorem that can be imported here.  The nonlinear continuation must either
derive a small-data global estimate in this one-way/critical decomposition or
move to a carrier class with an existing nonlinear variational theorem.  The
activated `0040` Cao-ring route pursues the latter without treating a positive
periodic metric as a necessary condition on all candidates.
