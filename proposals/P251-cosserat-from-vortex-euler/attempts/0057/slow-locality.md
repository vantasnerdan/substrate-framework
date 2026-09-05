# Exact second slow jet without a decorrelation premise

This representation uses the eleven incompressible-affine moment
constraints constructed in 0059: three translations and all eight
tracefree linear maps. Let F=xi cross omega be the force of a compact
selected column, in local patch coordinates y. The constraints imply

    integral F_i=0,
    integral y_j F_i=c delta_ij.

Indeed, the affine pairing is minus the contraction of the latter matrix
with an arbitrary tracefree matrix; its orthogonal complement is precisely
the scalar identity. This uses all eight affine constraints, not rotations
alone.

Choose a fixed smooth compact scalar chi0 with integral one and set
`f=-c chi0`, `G=F-grad f`. Then G has zero zeroth and ALL first moments.
This is an exact representation of the same force modulo a gradient:
`P G=P F`, `curl G=curl F`, and its helicity quadratic integral is unchanged
by integration by parts. It is not a subtraction from the physical action
or a different microscopic displacement.

## Explicit compact smooth double-divergence primitive

Take another fixed smooth unit-integral cutoff chi supported inside the
convex coherence ball. For each vector component define

    B_ijl(x)=integral_0^1 (1-t) integral integral
       chi(z) (y-z)_j (y-z)_l G_i(y)
       delta(x-(1-t)z-t y) dy dz dt.

Fourier Taylor's integral remainder about the center z is

    exp(-ip.y)=exp(-ip.z)[1-ip.(y-z)
      -p_j p_l (y-z)_j(y-z)_l integral_0^1(1-t)
                                      exp(-it p.(y-z)) dt].

Pairing with G kills its first two terms. Averaging over z therefore proves
`G_i=partial_j partial_l B_ijl` exactly, including its sign.

The average over Taylor centers is important. A single-center primitive
can be singular; this one is smooth and compact. Its support is contained
in the convex hull of the supports of G and chi, hence in the same finite
coherence ball. For `0<=t<=1/2`, integrate out z and differentiate chi,
whose argument denominator `1-t` is bounded below. For `1/2<=t<=1`,
integrate out y and differentiate G, whose denominator t is bounded below.
All polynomial factors are bounded on the ball. For each finite m this
gives a finite explicit bound of the form

    ||B||_C^m <= C_m(R) [||G||_L1 ||chi||_C^m
                                      +||chi||_L1 ||G||_C^m].

For example, in the fixed reference length units, expanding derivatives
of the quadratic factors gives the sufficient constant
`C_m(R)=2^(m+10) (m+1)^2 (1+R)^2`. The precise optimum is irrelevant;
the displayed finite bound follows directly from the two half-interval
representations. Good-patch C^m bounds therefore give uniform primitive
bounds. No compact-divergence inverse is imported or presumed.

## Modulated patches and the complete Leray quadratic form

For an accepted patch center X_a give its amplitude the Bloch factor
`exp(ik.X_a)`, constant on that patch. The exact global force, modulo the
compact gradient above, is

    G_total=partial_j partial_l sum_a exp(ik.X_a) B_a,ijl(x-X_a)
           =exp(ik.x) (partial_j+ik_j)(partial_l+ik_l) B_k,ijl(x),
    B_k(x)=sum_a exp(-ik.(x-X_a)) B_a(x-X_a).

Here B_k is stationary under joint field/grid translations. Since all
local coordinates have |y|<=R, its dependence on k is analytic as a
bounded operator from square-integrable patch amplitudes into every
controlled stationary Sobolev space. Coefficients are not differentiated
as though exp(ik.X_a) were exp(ik.x) inside a patch.

Let `L(p)B=P(p) p_j p_l B_ijl`. The FULL kinetic quadratic symbol is
`A(p)=L(p)^*L(p)`. It is smooth away from zero and homogeneous of degree
four. Defining it as zero at p=0 makes it C^3 there; in particular,

    ||partial_p^alpha A(p)|| <=C_alpha |p|^(4-|alpha|), |alpha|<=3.

Thus its second derivatives are bounded by C|p|². The spectral integral
for the stationary kinetic Gram is finite and twice differentiable using
the uniform Sobolev bounds of B_k. There is no infrared singular second
derivative to regularize, even if the stationary spectral measure has an
atom or long-range correlations. The same bounds yield operator-norm C^2
dependence for the patch-amplitude momentum block (one extra controlled
derivative bounds the high-frequency tail uniformly).

The helicity and KKS parts are local support integrals with smooth bounded
phase factors, so their second slow jets also exist. Their entire finite
matrices/operators, including momentum-gradient and odd first-gradient
terms, must be differentiated before Schur elimination and ensemble
pairing. Coercivity makes the inverse momentum operator C^2 as well, by
the usual exact identities `dP^-1=-P^-1(dP)P^-1` and its second derivative.

This establishes the required second-order slow-action locality for the
specified affine ensemble without assuming statistical mixing or deleting
Leray tails. It is not an assertion of an all-wave-number local PDE.
