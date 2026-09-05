# Independent review: stationary Beltrami tube angle action

Reviewer: `/root/smooth_core_review`, distinct from the author of attempt 0040.
Date: 2026-09-05. This is a new, bounded scientific transaction; it does not
reopen the previously reviewed 0036 construction. The review follows the
AGENTS contract and physics/small-ratio skills already read by this reviewer.

## Positive decision and frozen scope

The explicit stationary Beltrami tube has the claimed positive two-generator
Euler coadjoint-orbit quadratic action. The noncircular core angle is a genuine
observable of the cross-section jet, the conjugate incompressible shape has a
nonzero physical KKS pairing, and eliminating that shape gives positive angle
inertia from the same Euler action. The complete restricted Hessian, rather
than one selected diagonal, is positive in the stated open parameter regime.
The localized whole-space action and its continuity under fixed-eigenvalue
Beltrami approximation are also supported by the given analytic estimates.

This acceptance recommendation concerns the explicitly constrained affine
sector. It neither assumes nor proves that the unrestricted Euler PDE leaves
the two-dimensional trial surface invariant. Such invariance is not added as
an acceptance requirement. Common-angle transport, nonuniform material
translation and identification with a prescribed closed EPS knot remain
separate parent constructions.

The frozen input consists of `README.md`, `fourier_orbit.py`,
`verify_positive_sector.py`, the saved first scientific success `stdout.txt`,
and the candidate-generation records in this directory. Base provenance is
checkpoint 3626fbf and the proposal's v0.171.0 accepted release. No accepted
registry entry or unrelated implementation changes in this transaction.

## Evidence and independent derivation

The strongest oracle is exact variation of the Euler orbit energy combined
with finite Fourier convolution and orthogonal projection. The Fourier script
is an exact algebraic implementation, not a mesh or truncated approximation
to an unknown field: every field used here has explicitly finite Fourier
support. The saved 30/30 receipt is reused after code inspection; the prior
structural-equality diagnostics do not change a scientific equation.

The following checks close the load-bearing sign, normalization and scope
questions.

1. **Coadjoint tangent and full second variation.** Advecting vorticity by a
   volume-preserving generator gives `delta omega=curl(xi cross omega)` and
   the coadjoint velocity variation `delta u=P(xi cross omega)`. Repeating the
   same generator gives `delta²u=P(xi cross curl(delta u))`. At fixed
   `curl u=lambda u`, stationarity of the energy on the orbit follows from
   `u cross omega=0`. Pairing the second variation with u and moving the
   gradient part through the divergence-free curl gives
   `E''=rho integral (|delta u|²-delta u dot curl(delta u)/lambda)`.
   This derivation includes the second-order change of the advected field;
   it is not the positive kinetic norm alone.

2. **KKS sign and canonical convention.** With
   `Omega(xi,eta)=rho integral omega dot (xi cross eta)`, the identity
   `Omega(u,eta)=rho integral u dot (eta cross omega)=dE(eta)` gives the
   Euler generator under `i_X Omega=dE`. For `Omega=B dq wedge ds`, choosing
   p=Bs makes `p qdot-H` the corresponding action. Direct integration of the
   displayed trigonometric generators gives the KKS coefficient without
   reusing the Fourier result: the terms even in z reduce to
   `2(b-a) sin²x sin²y cos²z`, while the other terms average to zero. Its
   cell average is `(b-a)/4`, yielding the stated physical
   `B=-rho ell(a-b)/4`.

3. **Projection and complete Hessian.** The Fourier representation convolves
   all modes exactly, differentiation multiplies by ik, and the nonzero-mode
   projection is `I-kk^T/|k|²`; the zero harmonic is retained. The real-field
   product average correctly pairs k and -k without a spurious complex
   conjugation. The independent direct-energy route in the verifier expands
   both ordered second variations for mixed entries. Its sign mutation and
   compressible-tangent mutation expose the intended wrong conventions.
   The derived diagonal matrix has
   `Hqq=rho[26ab-17(a-b)²]/240` and `Hss=5rho(a²+b²)/96`.
   Thus the stated inequality and a unequal to b supply precisely positive
   definiteness and a nondegenerate orbit form. At a=b the angle pairing
   degenerates; the statement correctly excludes that circular limit.

4. **Physical dimensions and elimination.** The generator displacement is
   ell times its dimensionless Fourier expression, while vorticity scales
   as velocity/ell. Consequently H has energy-density units and B has
   action-density units. In general the derived inertia is
   `6rho ell²(a-b)²/[5(a²+b²)]`; at a=2b this is
   `6rho ell²/25`. Eliminating s from the actual first-order action gives
   this positive coefficient with the unchanged restoring potential.
   The frequency `175b²/(288ell²)` follows from K/I rather than entering
   any generator or coefficient as an input.

5. **Observable angle and same-fluid reaction.** At the axis the q-generator
   has horizontal Jacobian `J cos z`. Therefore it rotates the principal
   axes of the noncircular quadratic cross-section at z=0 by q and at z=pi
   by -q. This supplies the physical relative-angle normalization chi=2q,
   dividing both K and I by four. The generator is not a rigid rotation of
   every point in a finite-radius cross-section, and that stronger statement
   is unnecessary. The s-generator is divergence-free including its axial
   component and horizontal return. Its reaction and energy are retained
   before its elimination; there is no imposed external strain reservoir.

6. **Translation boundary.** The zero translation rows computed in the
   extended restricted Hessian and KKS matrix are consistent with uniform
   translation symmetry and this particular internal pair. They do not
   manufacture a positive translational mass from a zero orbit spring. The
   separate Galilean increment requires restoring material momentum; the
   document correctly leaves the nonuniform coupled action to parent work.
   The zero harmonic of the selected internal tangents is the assertion
   directly verified here; it is not a general claim that every coadjoint
   variation preserves mean velocity.

## Localization and approximation audit

The large-cutoff construction preserves divergence exactly because it uses
`curl(chi_R A)` instead of multiplying a displacement by a cutoff. For the
finite Fourier fields the periodic vector and scalar potentials used in the
decomposition exist: the generators and projected tangents have zero mean.
The identity for the projected localized tangent follows by expanding
`curl(chi_R B)` and `chi_R grad pi`. Each error includes at least one cutoff
derivative, so its L2 norm and the L2 norm of its curl are O(R^{1/2}), compared
with the O(R^{3/2}) bulk norm. Cross terms in the energy and helicity are
therefore O(R²), while the leading quadratic matrix is O(R³). Direct local
KKS integration has the same relative O(1/R) bound. Integration by parts on
the finite nonzero Fourier modes proves the asserted weighted averaging.
This is an analytic sign-preservation estimate, not numerical evidence of a
soft mode. It establishes some sufficiently large finite localization radius
without claiming a numerical minimum radius.

For fixed supported generators in R³, the whole-space Leray operator is a
fixed bounded L2 operator. In a fixed-lambda Beltrami class,
`omega=lambda u` and `grad omega=lambda grad u`; hence the local C1 norm
controls all quantities entering the helicity form as well. This explains
why C1 background continuity is sufficient here although the uneliminated
formula contains a vorticity derivative. The same statement would need
additional regularity for arbitrary non-Beltrami backgrounds.

The inspected [EPS Theorem 8.3](https://arxiv.org/pdf/1210.6271) preserves the
curl eigenvalue and approximates on a compact set with connected complement.
A ball containing the fixed supports meets that hypothesis, so the positive
localized action transfers to a decaying global R³ Beltrami field. This does
not identify the selected axis segment with a prescribed closed knot.

The separately cited [compact-manifold Theorem 2.1](https://arxiv.org/pdf/1505.01605)
is a rescaled local approximation theorem. Changing from R³ to a torus or
sphere also changes the ambient Leray operator in the kinetic term. Local
background approximation alone does not prove convergence of that nonlocal
term. Its use as a compact-domain action-preservation theorem therefore needs
a distinct Green/Leray rescaling estimate; otherwise it is a future route,
not a consequence already supplied by the present R³ continuity argument.

## Findings and correction boundary

The positive periodic action and localized R³ theorem have no blocking
mathematical finding. Two minimal exposition corrections were requested:

- Make the already intended local core-jet angle explicit, avoiding a claim
  of exact rigid rotation over an entire finite cross-section.
- Distinguish proved fixed-ambient R³ action preservation from the additional
  nonlocal-operator transfer needed for compact-manifold transplantation.

Neither changes a computed coefficient, generator or parent objective. The
second protects against a specific invalid implication: replacing a nonlocal
projection solely from local background convergence. A sufficient repair for
the flat torus is a growing-torus Fourier/Green convergence argument for
compactly supported sources; the sphere requires its own corresponding
metric and projection analysis. These are concrete continuation routes, not
reasons to discard the established whole-space action.

## Four-axis decision and continuation

- Verification: exact analytic and symbolic-verified action and localized
  R³ preservation; no numerical spectrum or stability of the full PDE.
- Review: audited; scientific acceptance recommended for those statements.
- Compatibility: compatible classical smooth-Euler extension under the stated
  constrained ensemble and fixed physical conventions.
- Epistemic: positive child route established; the parent continuum campaign
  remains active.

The reusable Fourier operations currently live in this attempt, so canonical
module extraction and affected-consumer testing remain implementation work
for the parent's promotion transaction. This mathematical review does not
itself perform a release change or claim that those integration steps passed.
No additional scientific review is needed for unchanged coefficients or
prose-only corrections at this boundary.

Initial substantive-review SHA-256 identifiers:

```
aa0078b747fbe657ce1357c861bc62e6a0ba250945e5c34e4560b02b12cae98b README.md
8fab3eea98d2d96bc42c43650baa840ce1e537413198c3358ee560c73a8ddf40 fourier_orbit.py
fff6b6c90f4047bd1e68be5d0fc5c759f04edb8f796ab826cebe8488dd330751 verify_positive_sector.py
3312d340edbabbfd52121f830af73fa955784175229b51958a5e43b70f193d5a stdout.txt
```

The corrected claim in this review is specifically rotation of the
noncircular core's quadratic cross-section jet. That statement is already
proved by the authored Jacobian calculation and changes no coefficient.
The authored receipt has not been edited by this reviewer. A parent summary
can consume this precise statement directly.

## Review-derived repair: flat-torus localization transfer

The following closes the nonlocal projection issue for an explicit flat-torus
construction, without asserting a transfer to the sphere or to a prescribed
knotted tube. It uses the known finite-Fourier seed itself; it does not need
an extension of EPS Theorem 2.1 from its stated unit ball.

Fix the finite localization radius R for which the whole-space two-generator
matrices are positive and nondegenerate. Write the resulting compact smooth
generators as eta_i and their compact sources as
`f_i=eta_i cross u_seed`, using curl eigenvalue one. Let T_L have side 2pi L,
with L an integer sufficiently large that the supports fit inside its
fundamental cube. The periodic seed is an exact eigenvalue-one Beltrami field
on T_L. Regard each compact source as periodized on this cube.

Use `fhat(k)=integral_R3 exp(-ik dot y) f(y) dy` and
`P(k)=I-kk^T/|k|²` for k nonzero. Parseval gives exactly

```
||P_L f||²_{L2(T_L)}
  = (2pi)^(-3) L^(-3) sum_{n in Z3} |P(n/L) fhat(n/L)|².
```

At n=0 retain the harmonic velocity, as the full coadjoint Leray projection
does: its contribution is `|fhat(0)|²/(2pi L)^3` and vanishes in the limit.
No additional fixed-mean constraint is imposed on the localized family.
For a compact smooth source fhat is
Schwartz. The remaining multiplier integrand is bounded near zero, smooth
away from zero and rapidly decreasing at infinity. These facts prove the
Riemann-sum limit, not merely a formal replacement of sums by integrals:
first remove a ball of radius delta around zero, whose integral and lattice
contribution are bounded by `C(delta³+L^(-3))`; on a fixed annulus use uniform
continuity; outside a large ball use a uniform summable Schwartz bound.
Take L to infinity, then delta to zero and the outer radius to infinity.
Consequently

```
||P_L f||² -> ||P_R3 f||².
```

Polarization supplies every mixed kinetic entry. The helicity identity
`integral P_L f_i dot curl(P_L f_j)=integral f_i dot curl f_j` is exact on
both domains, and these compact-source integrals are identical once supports
fit in the cube. The KKS entries are also identical compact local integrals.
Thus the full finite matrix H_L converges to the already positive H_R3, and
B_L equals its nonzero whole-space value. Positivity survives for sufficiently
large finite integer L. No convergence of an infinite-dimensional spectrum
or unrestricted invariant mode is claimed or needed.

Finally set y=Lx on a fixed physical flat torus of side 2pi. The explicit
stationary field is `u_L(x)=u_seed(Lx)`, with `curl_x u_L=L u_L`.
The physical generators are `xi_i(x)=eta_i(Lx)/L`; they are smooth,
divergence-free and supported in a contractible ball of radius O(R/L).
The exact physical matrix scalings are

```
H_physical = rho L^(-3) H_L,
B_physical = rho L^(-4) B_L.
```

The first comes from the volume element, since the displacement factor 1/L
cancels the vorticity factor L in delta u. The second additionally contains
two displacement factors and one vorticity factor. These positive rescalings
preserve Hessian positivity and KKS nondegeneracy. They imply inertia scaling
`rho L^(-5)` and frequency-squared scaling L², as expected from the physical
length scale 1/L. The fixed support contains the same local core-jet angle
observation before and after rescaling.

This is an established positive localized action on an explicit stationary
smooth torus Beltrami field, now with the nonlocal kinetic projection
controlled. A transplanted nonperiodic EPS background or sphere geometry
still requires its own matching background approximation and, for the
sphere, metric/projection analysis. The new result does not identify this
contractible support with a prescribed closed knotted vortex tube.

## Transaction close

The authored action and its coefficients are retained. The angle wording is
made precise in this review, and the potentially overbroad compact-domain
implication is replaced by the explicitly proved flat-torus construction
above. No author script or saved output changed. The parent may incorporate
these statements and continue its common-angle, translation and EPS geometry
obligations; the review does not imply campaign completion.
