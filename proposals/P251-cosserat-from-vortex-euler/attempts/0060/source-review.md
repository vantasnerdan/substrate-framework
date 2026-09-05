# Gaussian Beltrami import: independent source-applicability review

Reviewer: `/root/smooth_core_review`, distinct from the 0057 constructor.
Date: 2026-09-05. Contract: `README.md`, frozen before the source body was
read. Review checkpoint: `6b037a5e6e099b5c295ccc17eaabbbca5db7b029`;
parent release v0.171.0 and 0035 objective unchanged.

## Result

The primary source supports a smooth stationary Euler Gaussian ensemble
with finite pointwise velocity and derivative variances, translation
ergodicity, full local Beltrami support, positive-probability robust bounded
knotted-torus patches, and strictly positive lower spatial densities of
knotted invariant structures. This is a useful positive import for a
finite-energy-density alternative to the decaying EPS realization.

Its specified Gaussian law is not spatially SO(3)-isotropic, despite its
isotropic one-point velocity covariance. A global Haar rotation mixture
is stationary and isotropic but is not translation ergodic and is generally
not Gaussian. Time reversal alone leaves the original centered Gaussian
law unchanged. Reflection changes the curl eigenvalue sign; an equal
mixture of those two laws is not translation ergodic. These distinctions
matter when equating a spatial average in one realization to an isotropic
ensemble expectation.

This is an import/applicability receipt, not a review of an unfrozen
assembly or a completed Euler-to-Cosserat construction.

## Primary source and exact law

The source is Enciso, Peralta-Salas and Romaniega,
[Beltrami fields exhibit knots and chaos almost surely](https://arxiv.org/abs/2006.15033v1),
version 1 dated 26 June 2020. I inspected Theorem 1.2; Propositions 2.1,
2.5, 3.2, 3.4, 3.7, 3.8; Corollary 3.5; the torus definitions and
Propositions 4.6, 4.7 / Corollary 4.8; and Theorem 6.2 with its proof.
These cover the requested propositions and their actual support and
robustness bridges. They are one related applicability group, not multiple
independent validations of a parent claim.

The paper fixes curl u=u on R³. Its Gaussian field is obtained by applying
the explicit polynomial polarization p to scalar spherical white-noise
coefficients: u=U_(phi p), with independent real Gaussian harmonic
coefficients in phi and the required Hermitian symmetry. Its spectral
matrix density is `p(n) p(n)^* d sigma(n)` on the unit sphere. The complex
conjugate in that expression is essential; text extraction can obscure
the overbar in the paper's notation.

Every sample in its full-measure regularity set solves the same fixed
curl equation. A spatial rescaling and amplitude scaling can change the
positive curl eigenvalue and velocity normalization, but the length,
variance, derivative and density constants then rescale too. A curl-sign
change is a reflection, not the time reversal u -> -u.

## Smooth stationary Euler field and finite variance

Proposition 3.2 proves almost-sure convergence in C^k on every compact
set, for every finite k, and hence C-infinity regularity. Its spectral
random distribution is not L² on the sphere; the field is not a finite
spherical-harmonic realization. The exact curl relation implies div u=0
and `(u.grad)u=grad(|u|²/2)`. Thus it is a smooth time-independent Euler
velocity with pressure `p=-rho |u|²/2+constant`. Time stationarity of this
solution and statistical spatial stationarity of its law are different
facts, and both hold.

Proposition 3.4 gives translation-invariant covariance. Corollary 3.5
normalizes `E[u(x)u(x)^T]=I`, so E|u(x)|²=3 and the expected kinetic
energy in a bounded measurable D is `3 rho |D|/2`. Compact spectral
support also gives finite pointwise second moments of every derivative:
the differentiated covariance integrates bounded frequency monomials
against this finite spectral measure. These are exact finite-variance
statements, not a bound on every sample over all of R³.

Proposition 3.7 proves translation ergodicity of this Gaussian law from
its atomless spectral measure. Applied to |u(0)|², it gives the almost-sure
and L¹ spatial energy-density limit 3 rho/2. In particular the total
whole-space kinetic energy is infinite almost surely. Finite variance
does not mean finite whole-space energy, and the source supplies neither
the O(1/r) decay nor the finite angular expansion of the older EPS field.
Previous proofs using those far-field hypotheses require a new bridge
before application to this random field.

## Full support and bounded good patches

Proposition 3.8 states full support in the C^k compact-open topology on
the space of global curl-one fields: for a specified global Beltrami
template v, compact K and epsilon>0, the event
`||u-v||_(C^k(K))<epsilon` has strictly positive probability. Its proof
combines finite-harmonic approximation with a positive-probability finite
Gaussian coefficient event and an independently controlled tail.

For a template defined only near K, Proposition 2.5 supplies the needed
global approximation when R³ minus K is connected. In particular a
closed ball is an admissible K. One cannot replace this geometry by an
arbitrary locally defined template on a compact set with disconnected
complement without checking another approximation theorem.

Any strict finite local inequalities continuous in C^k persist in a small
enough template neighborhood. Such an event simultaneously bounds the
C^k norm by `||v||_(C^k(K))+epsilon` and has positive probability. Thus
bounded good patches are supported, provided their entire stated test is
indeed local and robust in that topology. The source does not make an
unproved global response, interface compatibility or nonlocal assembly
condition local merely by calling it a good-patch criterion.

For a measurable good-patch event A of positive probability, Proposition
3.7 applied to its bounded indicator gives

```
|B_R|^-1 integral_(B_R) 1_A(tau_y u) dy -> P(A)>0
```

almost surely and in L¹. This is positive volume density of good patch
centers, not a theorem that different patches are independent. A bounded
radius permits a separated-subfamily packing lower bound if needed, but
a chosen translation-equivariant point process or material partition is
additional construction, not part of the quoted support statement.

## Knotted invariant tori and their spatial density

The source's tori are invariant embedded surfaces with nonvanishing flow.
Its torus robustness uses Diophantine dynamics and nonzero twist, not an
arbitrary invariant surface. Proposition 4.7 obtains a finite-harmonic
Beltrami template containing a positive inner-measure family of such tori
of any prescribed knot type. Its imported deterministic construction is
Enciso--Peralta-Salas, *Existence of knotted vortex tubes in steady Euler
flows*, Acta Mathematica 214 (2015), especially Theorem 1.1; the local
approximation in Proposition 2.5 invokes that work's Theorem 8.3. This is
the EPS source already used in the campaign, not an unrelated torus model.

Proposition 4.6 and Corollary 4.8 supply the crucial C^k robust open event
for k>=4 among divergence-free fields. It includes an invariant closed
solid torus whose boundary has the stated dynamical properties and whose
interior contains a positive inner-measure family of invariant tori.
Consequently full support produces a strictly positive-probability bounded
patch with a genuine material invariant tube. Positivity is not inferred
from the probability of one exact prescribed field, which would be zero.

Theorem 1.2 gives an almost-sure strictly positive lower bound for the
large-ball ratio of the inner volume occupied by ergodic invariant tori
of a fixed isotopy type to the ball volume. Here ergodic on a torus means
its trajectories are dense; it is distinct from translation ergodicity
of the random-field law. The paper deliberately uses inner measure for
the union of tori. Its claim is a liminf lower bound, not an exact density
limit or a prescribed volume fraction.

The stronger Theorem 6.2 additionally gives a positive liminf density of
pairwise disjoint invariant closed solid tori, each containing a fixed
positive inner volume of the prescribed invariant-torus family, with
frequency/twist ranges chosen by the construction. This is the appropriate
source for a count of finite material tubes; counting individual nested
invariant surfaces need not be finite. Countably many knot types can be
handled simultaneously on a common probability-one set. The source does
not assert that all fluid is partitioned by those solid tori.

## Source law is not spatially isotropic

The one-point covariance I implies a rotationally invariant one-point
Gaussian distribution, but does not imply rotational covariance of the
entire field. In the prescribed law, p vanishes at the distinguished
two poles of the sphere and is nonzero elsewhere. Therefore the continuous
spectral trace |p(n)|² is not constant on the sphere. A spatially isotropic
vector-field law would require its spectral trace measure to be rotation
invariant. Thus this particular source law is anisotropic.

This is a direct consequence of the source's explicit p and spectral
measure, rather than an inference from an absent isotropy claim.
Remark 3.1 permits certain nonvanishing scalar changes of polarization;
it does not identify the source's law with its Haar mixture or prove
ergodicity for arbitrary mixtures.

## What rotations, reflection and time reversal preserve

For fixed Q in SO(3), set `u_Q(x)=Q u(Q^T x)`. Its law remains Gaussian,
smooth, stationary, translation ergodic and curl one, with covariance I
at each point. Conjugating translations by Q proves its ergodicity.
The topology and ball-volume density assertions are also preserved.

If a single global Q is drawn from Haar measure and then an entire field
u_Q is drawn, the resulting mixture is spatially SO(3)-isotropic by a
Haar change of variables and is stationary. It retains finite variance,
smooth stationary Euler samples, full local support and almost-sure
positive-density knot properties. However it is a mixture of distinct
Gaussian covariance laws, not in general a Gaussian law.

It is also **not translation ergodic**. For example consider the empirical
spatial two-point statistic of `u(x).u(x+h)`. Conditional on Q, its limit
is `tr kappa(Q^T h)` by componentwise ergodicity. Anisotropy of the spectral
trace implies that for some h this is a nonconstant function of Q.
The spatial limit is then a nonconstant translation-invariant random
variable under the mixture. Consequently a single realization's general
spatial averages cannot be replaced by the mixture's isotropic expectation.
For an orientation-dependent good event its conditional spatial density
can likewise depend on Q, although it remains positive. Rotationally
invariant knot-density statements survive without needing mixture
ergodicity because they already hold for every fixed rotated component.

Time reversal u -> -u preserves the original centered Gaussian law itself.
It leaves its translation ergodicity unchanged and its invariant sets
unchanged, merely reversing trajectories. For a frozen microscopic
realization it still has the separate physical role of circulation
pairing; equality of probability laws does not impose coherent coordinate
or reaction-momentum ties.

For an orthogonal F with det F=-1, the physical polar velocity pushforward
`u_F(x)=F u(F^T x)` satisfies curl u_F=-u_F. Its law is still stationary,
smooth, finite-variance and translation ergodic as an individual component.
An equal mixture with the positive-eigenvalue law is not translation
ergodic: the translation-invariant event curl u=u has probability one-half.
It also need not be Gaussian. Orthogonal Haar averaging gives a stationary
O(3)-isotropic mixture with these same component distinctions. Reflection
changes knot chirality, but the source supports every knot and mirror
type, so the positive-density topology conclusions persist componentwise.

If one requires both full isotropy and translation ergodicity of a single
law, globally mixing rotated source realizations does not supply that
object. A separately constructed isotropic Gaussian spectral law, or a
carefully specified componentwise versus ensemble averaging procedure,
would need its own support/applicability argument. This review does not
silently substitute such a new law for the cited one.

## Applicability disposition and hashes

Keep this import for smooth stationary finite-energy-density Euler fields,
full local support, robust positive-probability bounded tube patches and
positive lower spatial density of knotted invariant structures. Preserve
the exact law and distinguish original translation ergodicity from global
Haar/reflection mixture invariance. No source passage by itself supplies
the common material action, compact-cell assignment, constituent moduli,
or the parent continuum theorem.

- Archived `../0057/source-2006.15033.pdf` SHA-256:
  `57d83ba2178807042fd0b0c3e9cb4cc9e4d39921c36e40cda1c82bb8c841b5a4`.
- Archived `../0057/source-2006.15033.txt` SHA-256:
  `71ef04bcbb3bf89b15cbfa79295bb085d4ab032489f3ea163f7af96ddd646d7c`.

No numerical oracle or sampled eigenvalue enters this source review.
The positive import is established with the stated applicability boundary;
there is no remaining correction check unless a downstream attachment
misstates the law or its averaging operation.
