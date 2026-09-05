# Explicit slow material continuation of the Galilean pair

This append-only candidate uses the actual material representation. Its
gradient form is computed here; the positive coadjoint H of 0085 is not
substituted for it. It complements the full orbit/reconstruction repair
in 0095 and the stationary normalization in 0098.

## Exact linear volume-preserving map

Let U(t,x) be any solenoidal slowly varying displacement. The actual
near-identity volume-preserving map is the spatial flow of U at each time,
to first order `Phi=id+epsilon U`. Its material realization is
`g=Phi composed with g0`, where g0 is the stationary Euler flow. Thus

    delta u = U_t + (u0.grad)U - (U.grad)u0,
    delta chi_a = -U.grad chi_a.

These include the microscopic return/advection terms, not just U_t.
For one Fourier macro amplitude `U=a(t) exp(i k.x)`, `k.a=0`, the formula
is exactly

    delta u = exp(i k.x) [adot + i(k.u0)a - (a.grad)u0].

The displacement is divergence free without a vector-potential gauge;
the Eulerian velocity is divergence free because it is its time derivative
plus the bracket of two divergence-free fields. Every invariant tube tag
and the ambient tag are moved by this SAME field.

At k=0, U=X0+Vt gives the exact translation/boost tangent described in
`material-cotangent-bridge.md`. Constant V is a curl-free exact one-form
on R3. At nonzero k the boost datum has vorticity `i k cross adot`; it
cannot silently be treated as an independently variable harmonic velocity
on one fixed zero-mean vorticity leaf.

## The complete mean material quadratic action

The actual incompressible material Jacobi action about a stationary Euler
state is

    L2[eta] = (rho/2) integral |eta_t+(u0.grad)eta|^2
                 - (1/2) integral eta.Hess(p0).eta,
    div eta=0.

Pressure is the constraint multiplier of the original global material
action. The second term is retained before averaging. For the stationary
isotropic law in 0098, write

    E[u0]=0,        E[u0_i u0_j]=sigma^2 delta_ij,
    E[partial_i partial_j p0]=0.

The last equality follows from translation stationarity and the finite
derivative moments in the bounded-event/finite-variance construction;
it is not a pointwise statement about an EPS tube. Since U is a prescribed
slow macroscopic test field independent of the translated microscopic
sample, direct substitution and expectation give exactly

    E[L2[U]] = integral [rho |U_t|^2/2
                         + rho sigma^2 |grad U|^2/2].

Thus the genuine zeroth mass is rho. The bare material gradient term has
the displayed POSITIVE Lagrangian sign, i.e. NEGATIVE elastic stiffness.
It is not the positive coadjoint affine shear of 0043/0057. Positive
material STF cages in 0094 are a concrete way to repair that computed
gradient form without changing representation. This exact calculation
prevents promotion of a false same-action shear coefficient.

The formula already retains every order-k and order-k^2 term of this
material lift: the order-k mean gyroscopic cross vanishes by E[u0]=0,
the k^2 coefficient is the displayed covariance, and the pressure term
vanishes only after its actual stationary average. There is no asymptotic
tail or fitted box constant in this quadratic-in-U calculation.

## Adding actual compact material directions

For `eta=U+sum z_a Xi_a`, with compact divergence-free Xi_a and independent
cell amplitudes z_a, substitute the SAME expression into L2. The complete
mass, gyroscopic and stiffness blocks are respectively the polarizations
of

    rho integral |eta_t|^2,
    rho integral eta_t.(u0.grad)eta,
    rho integral |(u0.grad)eta|^2 - integral eta.Hess(p0).eta.

All cross blocks are retained before any Routh elimination. For a
spatially uniform translation a and a compact Xi,

    integral Xi=0,
    integral (u0.grad)Xi=0,
    integral a.Hess(p0).Xi=0.

The first identity uses compact divergence freedom, the second compactness
and div u0=0, and the third integration by parts against div Xi=0. Hence
the constant translation has no internal mass, gyroscopic or stiffness
cross in the full material action. This is a SAME-action derivation of
the leading decoupled rho block. At the next affine jet the mass cross is

    rho integral (h r).Xi
       = rho beta.integral r cross Xi,

with its STF part zero. Notice that this is the moment of Xi, not the
moment of the coadjoint induced velocity `Xi cross omega0`. This is the
specific reason that the moment-normalized material construction 0091
and the compact orbit construction 0085 cannot exchange their j values.

## Kelvin data and reconstruction, explicitly retained

For this actual material map, the first variation of the pulled-back
velocity one-form is, modulo exact forms,

    delta(g^* u_flat) = g0^* [delta u_flat - (eta cross omega0)_flat].

If its entire Kelvin datum is fixed, the bracketed one-form must be exact
with the required periods. For arbitrary initial circulation perturbation
it is instead the appropriately advected initial one-form. Conservation
of that datum is an equation of the unrestricted material Euler action;
it is not automatically an equation of a finite restricted trial action.
The explicit Fourier lift above allows its residual to be calculated
without treating a small-k circulation change as a new rotor mass.

Accordingly this construction establishes an exact material Cauchy--Born
pullback, actual centroid/tag map, leading rho mass, and all its mean
macro gradient terms. Completing its fixed-Kelvin reduction requires the
corresponding complementary variables or the explicitly declared closure
reactions. The original conditional affine scope permits deriving this
restricted action; it does not make two different restricted actions
identical. The full physical spin/current map is retained in 0091 for the
material route and in 0095 for the reconstructed orbit route.
