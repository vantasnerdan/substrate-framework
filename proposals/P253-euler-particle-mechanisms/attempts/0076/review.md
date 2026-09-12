# P253/0076 independent joined review of 0072 and 0075

## Review identity, integrity, and boundary

This is the preregistered non-author/non-implementer review of the final
author-stage `0072` and `0075` artifacts. The two units were adjudicated
separately. The corrected review contract had SHA-256
`0c5abac3bf36a2d6faab55591a806541f680c955c9c34b1078c6728d32462457`;
central activation reported `WORKFLOW VALID` and exit `0` before a target body
was opened. All frozen claim-bearing, verifier, API, and test hashes matched at
opening.

The earlier README SHA
`789021ca36f17f86a0bc3cca55d0c750447e8efed9bf040bc50ec7ec0e0f1d76`
is preactivation-only, transparently superseded provenance. It pinned `0075`
before the exact stationary-residual continuation. No scientific evidence was
consumed under that stale inventory.

Active `0074` remained unopened and supplied no definition, lemma, test, or
conclusion. This review makes no P5, electron, neutrino, scalar-charge,
particle, force-law, action-selection, or all-time-persistence conclusion.

## Source applicability

The only new primary mathematical source used by the targets is Roman
Shvydkoy, *Homogeneous solutions to the 3D Euler system*,
arXiv:1510.03378v1, cached with SHA-256
`71b277e2b78c2e8a1d14994c7267458795ee630bd555e0cb856fcc51402f5706`.
Direct inspection confirms that its equations (3a)--(3c) give, for
`V=r^(-alpha)(v+f n)`,

    (2-alpha)f+div_S v=0,
    v.grad_S f=|v|^2+alpha f^2+2 alpha p,
    (1-alpha)f v+nabla_v v=-grad_S p.

At `alpha=2` these are exactly the sphere equations used by the attempts.
Proposition 5.2 classifies the axisymmetric `alpha=2` case only, and the paper
explicitly leaves the general case open. The targets respect this boundary:
their nonaxisymmetric fixed-frame result comes from a direct residual, not an
imported classification.

Standard whole-space Sobolev Euler local well-posedness is used by `0075`
only for the ordinary `H^(s+1)` solution and common lifespan. The weighted
tail persistence is proved in the attempt rather than attributed to an
external theorem.

## Unit A — final 0072

### Harmonic, flux, and smooth-completion calculus

For

    u=r^(-2)[f(n)n+T(n)]+R,

the radial divergence cancels and the leading divergence is
`r^(-3) div_S T`. Hence `T` is a toroidal vector harmonic and the radial
coefficient is unrestricted except for the global flux row. Smooth
source-free completion requires `integral_S2 f=0`: the excluded invariant
tail `q n/r^2=q x/r^3` has divergence `4*pi*q*delta_0`. Conversely, for each
zero-mean radial harmonic, the annular field

    -(chi'(r)/r) grad_S Delta_S^(-1) f

cancels the cutoff divergence exactly; toroidal harmonics can be cut off
directly. This proves necessity and sufficiency at the stated smooth
finite-energy datum scope. The radial and toroidal parity assignments and the
absence of an `SO(3)`-trivial transverse source-free line are correct.

### Fourier transforms and translated cross energy

With `fhat(k)=integral exp(-i k.x)f(x)dx`, the homogeneous scalar coefficient
quoted in `0072` has the standard gamma-function normalization. Differentiating
the degree-three scalar transform gives

    FT[r^-2 Y_l n]=i*c_(l,3)|k|^-1 grad_S Y_l,

while applying `x cross grad` to the degree-two scalar transform gives

    FT[r^-2 n cross grad_S Y_l]
      =c_(l,2)|k|^-1 khat cross grad_S Y_l.

Thus `A_1=pi^2` and `B_1=-4*pi*i`. Parseval with the translated second field
`u_2(x-d)` supplies `exp(+i k.d)`. Inverting the transverse `|k|^-2` tensor
gives `(I+dhat tensor dhat)/(8*pi*d)`. These facts reproduce all four entries
of the displayed `l=1` multiplicity block, including the `-2*pi*rho_0/d`
mixed chiral sign. The factor `rho_0`, and the absence of an extra `1/2` in
the cross term from `rho_0/2 integral |u_1+u_2|^2`, are correct.

The weighted remainder is strong enough to make the nonhomogeneous Fourier
part bounded at zero; crossing it with a `|k|^-1` tail is one power faster.
The general angular product therefore produces the asserted leading `d^-1`
kernel with the orientation tensor retained.

For `u_a=(a cross x)/(1+r^2)^(3/2)`, direct curl calculus gives the stated
vorticity, zero helicity density, and norm
`||u_a||_2^2=pi^2|a|^2/2`. Its cross energy is

    (2*pi*rho_0/d)[a.b+(a.dhat)(b.dhat)]+O(d^-2).

It is positive and decreasing for aligned like amplitudes but differs by a
factor two between parallel and transverse separation. Calling it a
repulsive-sign effective potential is therefore properly conditional on a
mechanical separation coordinate; it is not an Euler force theorem.

### Fixed-frame Gaussian channel and covariance boundary

The symbol

    F(n)=[P_n e1+i n cross e2]/sqrt(1+n3^2)

is transverse, obeys `F(-n)=conj F(n)`, and has unit Hermitian norm. The
Gaussian field `q exp(-sigma^2 k^2/2)F(khat)/|k|` is real, smooth, in every
finite Sobolev space, and finite energy. Its same-frame translated cross
energy is exactly

    rho_0 q1 q2 erf(d/(2 sigma))/(4*pi*d).

The proper rotation `R_z(pi)` sends the symbol to its negative. The target's
averaged-inner-product argument correctly shows that a nonzero transverse
deterministic symbol cannot give the same positive scalar coefficient for all
independent relative rotations: equality would force a rotation-equivariant
field `c n`, which transversality then kills. This establishes a locked-frame
or orientation-carrying channel, not scalar electric charge. It does not
exclude a physical frame lock, independent sectors, or an autonomously
controlled orientation average.

### Stationarity, evolution, and observable domain

The complete radial/toroidal `l=1` homogeneous family has the displayed
nonzero `curl((u.grad)u)` components, so no nonzero member is steady. For the
minimal fixed-frame Gaussian realization, the necessary stationary stress
condition follows from the `k -> 0` transverse Fourier projection:
`P_n M n=0` for every `n`, hence the symmetric matrix `M` must be scalar. The
computed angular stress is not scalar. This refutes that radial form factor,
not every possible faster core.

For the declared weighted classical interval, `u=O(r^-2)` makes
`u tensor u` integrable with borderline logarithmic first moment. The Newton
split gives pressure `O(r^-3)`, gradient `O(r^-4)`, and transport `O(r^-5)`,
so `u_t=O(r^-4)` and the complete inertial-coordinate `r^-2` coefficient is
constant. The two pressure derivatives are load bearing; there is no pressure
monopole or dipole. `0072` correctly states this result conditionally because
unweighted Sobolev theory alone does not construct the weighted interval.
`0075` independently supplies that missing propagation theorem.

The ordinary velocity and vorticity integrals, toroidal angular momentum, and
absolute vorticity impulse diverge for the representative; its helicity is
zero and kinetic energy finite. Independent translation of one of two
overlapping noncompact summands is not a tangent of the standard
identity-at-infinity compact-carrier orbit. A finite or renormalized
asymptotic KKS phase space, or material labels with a joint action, remains a
real additional construction.

**Unit A verdict: established as scoped.** Its strongest exact result is a
smooth finite-energy, orientation-carrying degree-minus-two Euler boundary
channel with exact `d^-1` kinetic cross energy and a conditionally conserved
asymptotic label. Its mechanical/KKS and scalar-particle exclusions are
necessary and correctly retained.

## Unit B — final 0075

### Stationary Route A

The originally preregistered explicit inverse-VSH coefficient summation was
not completed as a separate infinite table with a truncation remainder.
Instead, Route B obtains the same full inverse exactly as an anisotropic-Riesz
kernel and resolves the scientific stationarity question without truncation.
Accordingly, Route A is not independent evidence and earns the route-scoped
verdict **blocked as an unexecuted representation, with its target question
resolved exactly by Route B**. No claim depends on completing that redundant
coefficient table.

### Stationary Route B

Writing

    a(k)=sqrt(k1^2+k2^2+2 k3^2)

turns the complete symbol into

    F(khat)/|k|=P_k e1/a+i(k cross e2)/(|k|a).

This is an exact inverse, not a harmonic truncation and not an illegitimate
treatment of the oscillatory `i^l` multiplier as a classical symbol. The
anisotropic inverse of `a^-1` is

    phi=s^-2 C/(1+t^2/2),  C=1/(2*pi^2*sqrt(2)),

and the Feynman-parameter formula for `psi` gives the stated even profile
`g`. The equation `-Delta chi=phi` gives
`-(1+t^2)h''-t h'=f`; its logarithmic coefficient
`c=-C*sqrt(2)*atanh(1/sqrt(2))` is negative. The signs
`g0=I0/(4*pi^2)>0` and `g2=-I1/(4*pi^2)<0`, with `0<I1<I0`, follow directly
from differentiating the integral.

An independent symbolic reconstruction from the even equatorial jets, without
calling either target verifier or public API, produced

    x1^6 [curl((u.grad)u)]_2
      =g0(-10c+7C)+g2(-5c+3C).

Substituting the exact `g0,g2` rewrites this as

    [(I0-I1)(-5c+3C)+I0(-5c+4C)]/(4*pi^2)>0.

Thus the full fixed-frame degree-minus-two texture has a nonzero homogeneous
degree-minus-six stationary residual. No faster-decaying core can cancel that
leading order. This verdict is exactly limited to the fixed-frame texture;
it is not a general nonaxisymmetric `alpha=2` no-go and does not affect Route
C.

**Route B verdict: fixed-frame stationary texture refuted by an exact positive
residual.**

### Affine weighted-vorticity Route C

For integer `s>=4`, `0<alpha,gamma<1`, the attempt fixes a smooth global
representative `U_*` and writes `u=U_*+R`,
`q=curl R`. The phase space requires `R in H^(s+1)`,
`q in Z_(3+gamma)^(s-1,alpha) intersection L1`, and `integral q=0`. The last
row is both the Hodge far-field cancellation and a propagated constraint.

The zero-mean Biot--Savart estimate

    B: H^s intersection Z_(3+gamma)^(s-1,alpha) intersection L1
       -> H^(s+1) intersection Z_(2+gamma)^(s,alpha)

is supported. At low Fourier frequency, splitting the integral at
`|y|=|k|^-1` gives `|qhat(k)| <= C|k|^gamma`; high frequency is the ordinary
order-minus-one multiplier estimate. In physical space, subtraction of
`L(x)` uses `integral q=0`. The inner region combines the kernel difference
`O(r^(-3-j)|y|)` with the truncated first moment `O(r^(1-gamma))`; the local
ball uses odd-kernel Schauder cancellation; comparable and exterior regions
give the same or faster weight.

The top Hölder row also closes. For point pairs separated by a fixed fraction
of the annulus radius, the pointwise derivative bounds divided by
`|x-x'|^alpha` give the required scaled seminorm. For nearby pairs, the
far-source part uses the mixed difference in both `x` and `y`, gaining the
extra kernel derivative after the zero-mean subtraction, while the ball about
`x` is the standard `C^(s-1,alpha) -> C^(s,alpha)` Biot--Savart estimate.
This yields the claimed `r^(-2-gamma-s-alpha)` top weight without assuming an
unavailable extra derivative of `q`.

The exact remainder equation

    q_t+u.grad q-q.grad u=-u.grad Omega_*+Omega_*.grad u

contains no pressure. After `s-1` differentiated commutators, each top term
uses at most `D^(s-1)q` and `D^s u`, the latter supplied by the Hodge gain.
The weighted Hölder difference quotient has `grad u` multiplying the top
quotient plus lower products. All fixed-tail source terms decay at least
`r^-6`, faster than the propagated `r^(-3-gamma)` row. There is therefore no
hidden derivative loss.

Integrating each component is licensed by the declared decay. Transport,
stretching, and both fixed-tail source terms vanish after integration by parts,
so `integral q=0` is conserved. Pressure is reconstructed independently as a
lower-order consequence with the stress multipole and logarithmic first
moment; it is not used to repair the top vorticity estimate and supplies no
finite-speed claim.

The approximation is also admissible. With
`q_N=curl(chi_N R)`, every approximant is divergence free, compactly
supported, and exactly zero mean. The cutoff curl term has the asserted
`N^(-3-gamma-j)` transition size. Compact convolution preserves divergence
and the integral; no independent scalar cutoff of `q` is used. Since
`R_N=Bq_N=P_L(chi_N R)`, the data converge in `H^(s+1)` with uniform Sobolev,
weighted-Hölder, and `L1` bounds. Standard Sobolev Euler theory therefore
gives a lifespan depending on those uniform bounds, not the support radius or
mollification scale. Local compactness plus the uniform tail estimate passes
the weighted and `L1` rows to the unique Sobolev limit.

The conclusion is consequently a genuine datum-dependent local Euler theorem:
the full `r^-2 U(n)` coefficient is fixed in inertial coordinates on the
constructed interval. The theorem does not rotate `U(n)` with the flow, split
overlapping tails into material carriers, or infer stationarity from local
persistence.

**Route C verdict: established as stated.**

## Oracle and API evidence

The exact target receipts report `0072`'s completed oracle passing and ten
focused API tests passing, and `0075`'s completed oracle passing and twelve
focused API tests passing. The recorded full checkpoint reports 2775 tests
passing. Inspection confirms that the verifiers derive transform constants,
kernel homogeneities, representative curls, stress anisotropy, and the
equatorial residual rather than hard-coding a final boolean. The public module
is consistently limited to algebraic helpers and domain/order ledgers. These
receipts support the algebraic surface only; the functional-analytic Route C
verdict rests on the written proof audited above.

No production numerical run was needed or performed for this review. The
independent residual check was exact symbolic algebra and did not invoke the
target verifier or API.

## Final joined verdict and next dependency

No bounded correction is required. Unit A is established at its exact
finite-energy, oriented-channel, conditional-mechanics scope. Unit B
independently refutes the fixed-frame stationary texture and establishes the
affine weighted local propagation theorem. These conclusions are compatible:
nonstationarity of one texture does not remove its locally persistent
inertial-coordinate tail label.

The strongest next construction is either a different angular texture or a
controlled periodic tail satisfying the full `alpha=2` sphere equation,
together with a finite or renormalized asymptotic KKS phase space (or material
labels) that makes relative translation an admissible mechanical coordinate.
That is the missing bridge from exact boundary-channel energetics to a
persistent carrier interaction; it remains outside this review.
