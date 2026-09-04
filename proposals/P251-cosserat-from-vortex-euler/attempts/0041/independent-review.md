# Independent review: complete finite-core Biot--Savart gradients

Reviewer: `/root/smooth_core_review`, distinct from the author of attempt 0041.
Date: 2026-09-05. This is one bounded scientific review under AGENTS.md and
the physics/small-ratio skills. The reviewed transaction is 0041 only;
previously reviewed smooth-core existence in 0036 is an unchanged input.

## Decision and strongest positive statement

The complete three-dimensional Euler kinetic-energy second variation for the
declared affine twist and common-bend families is established. In particular,
the two-triangle twist matrix is an exact positive-definite transverse-
projector Gram matrix, including finite-core and mutual-core contributions.
Its local zero-wavenumber coefficient is infrared finite and has the stated
controlled leading logarithm. The radial conjugate momentum carries its own
positive gradient coefficient; retaining it gives the stated optical
normalization correction. The constrained thin-core singular ratio is
1129/220, not the free-relative-sector factor ten.

No load-bearing defect was found, and no scientific correction or scope
reduction is requested. The common-bend response correctly retains its
infrared logarithm; this is a positive nonlocal response, not a falsely
finite local modulus. Stationary EPS identification and completion of the
parent continuum action are separate obligations, not conclusions of this
review. An unrestricted Euler invariant manifold or all-wavenumber local
Cosserat equality is not added as an acceptance requirement.

## Frozen boundary and evidence

The base provenance is checkpoint 3626fbf and the proposal's v0.171.0
accepted release. The reviewed statements and implementation are the README,
`full-biot-savart-twist.md`, `radial-momentum-gradient.md`, and the three
scripts `full_twist_energy.py`, `radial_gradient_elimination.py` and
`cage_gradient_matrix.py`. The saved exact receipts are 21/21, 8/8 and 7/7.
The label-only replay of the first script preserves the initial output and
does not alter an asserted equation. These receipts were inspected and reused,
not rerun to validate the same unchanged boundary.

The analytic documents are the exact proof of continuum positivity,
integrability and remainder control. The scripts corroborate the complete
variation, kernel differentiations, projector identities, Legendre
elimination and finite matrix algebra. They are not numerical discretizations
or independent proofs of the stated functional-analytic estimates.

The [NIST K0 integral representation](https://dlmf.nist.gov/10.32.E6) and
[K0 series](https://dlmf.nist.gov/10.31.E2) were inspected directly. They
support the axial Newton-kernel transform and the logarithmic expansion used
here. No empirically calibrated line tension or core cutoff is imported.

## Complete second variation and projector

The full pushforward retains more than the first vorticity variation. The
map `x=R_{q(z)}a, z=a_z` has unit Jacobian. Applying its derivative to the
vertical initial vorticity gives the stated exact horizontal component
`q' f(R_-q x) Jx` and vertical component `f(R_-q x)`. Differentiating twice
reproduces both `delta omega` and `delta² omega`, including the horizontal
`-2qq' B Jx` and vertical `q² L_z² f` terms. Their divergences cancel.

Since the background vorticity is vertical and independent of z, only the
vertical second variation pairs with it through the diagonal Newton
resolvent. It contributes the displayed k-independent diagonal correction
to the straight angular Hessian. A common constant rotation is an exact
energy symmetry; the second-variation term is required for that cancellation.
The separate resolvent difference in the first-variation longitudinal term
produces the negative subtraction in C(k). These are distinct physical
contributions, and the corrected mutation label identifies them accurately.

With `B_i=div A_i`, transverse Fourier transformation gives

```
|A|²/(|xi|²+k²)
 -(xi dot A)^*(xi dot A)/[|xi|²(|xi|²+k²)]
 = |P_T A|²/(|xi|²+k²).
```

Thus the gradient matrix is a Gram matrix of the complete projected
variation, not a positive expression selected by dropping a term. Its strict
positivity follows from compact support and disjoint triangle supports: a
vanishing norm would make `sum c_i A_i` curl-free separately on each
triangle. But `curl(f_i Jx)=div(x f_i)` has the nonzero moment
`integral |x|² curl A_i=-2 integral |x|² f_i`. This also excludes a hidden
twist relabeling direction. The k=0 claim uses the same argument after its
integrability is established, rather than inferring strictness from a
finite numerical sample.

## Infrared normalization and real-space kernel

The zero-mode logarithmic constant changes the straight energy only by the
fixed total-circulation constant. Every admissible variation preserves that
circulation. For twist, the stronger triangle first-moment cancellation gives
`integral A_i=0`, hence `Ahat_i=O(|xi|)` and `Bhat_i=O(|xi|²)`. This is
enough for both the k=0 Gram integral and the `G0²` pairing to be infrared
finite. It also removes the relevant polynomial ambiguity in the
biharmonic kernel. Smooth compact cores control high Fourier momentum.

Splitting the momentum integral at |xi|=|k| gives
`C(k)-C(0)=O(k² |log(kd)|)` for a fixed core profile. The resulting omitted
energy contribution is O(k^4 |log(kd)|), as stated. This derivative expansion
does not assert a uniform local law at arbitrary finite k.

The real-space kernel can be independently derived by writing
`t=|x-y|²`, `D=x dot y`, `A=x cross y` and
`Q=t log(t)/(16pi)`. Direct differentiation gives

```
L_x t=-2A, L_y t=2A, L_x L_y t=-2D,
L_x L_y Q=-A²/(4pi t)-D(log t+1)/(8pi).
```

Subtracting this from the logarithmic `A_i G0 A_j` term yields precisely

```
[-D log|x-y| + D/2 + A²/|x-y|²]/(4pi).
```

The D/2 integral and any change of logarithmic reference length vanish by
the two first-moment constraints. The remaining fraction is bounded by
`min(|x|²,|y|²)`; the only singularity is an integrable logarithm. This
checks the sign, factor and finite term independently of the scripted
two-angular-derivative calculation.

## Leading logarithm, mutual interactions and bend

For a self pair, the coefficient of log(d/epsilon) is exactly
`Gamma² |X_a|²`: all first-moment cross terms vanish. No finite self part
has been replaced by a Rankine constant. For normalized nonnegative g, the
negative logarithmic moment is bounded by
`M0 integral_|r|<1 (-log|r|)d²r=pi M0/2`; the positive part is bounded by
`log(max(2Lcore,1))`. This reproduces Mlog. The stated support separation
bounds all mutual logarithms by log(5/2), and the row count is three self
terms plus six diagonal mutual and nine off-diagonal terms. These estimates
give the displayed operator remainder bound. They do not delete any mutual
interaction and do not choose a fitted outer radius.

For common bend, `A_e=f e` has nonzero mean. Its projector matrix is positive
for nonzero k, and sixfold symmetry makes it scalar; its trace fixes the
factor one-half. Threefold symmetry forces every bend/twist mixed transverse
vector to vanish. The coefficient of the ultraviolet logarithm counts six
self circulations squared, while the infrared logarithm counts the square
of total circulation, giving 6Gamma² and 36Gamma² respectively. A constant
local bend modulus would require an additional fixed outer-coherence or
return-vorticity prescription, which the authored statement explicitly
retains as a separate construction.

## Radial momentum and constrained optical normalization

The relative-radius variation translates each core by `s_a X_a/2`; it does
not compress its area or alter its vorticity distribution. Its first moment
vanishes triangle by triangle. A nonzero constant translation of a compact
core has nonzero curl of `f_a X_a`, proving strict positivity of its Gram
coefficient. Reflection distinguishes the radial and angular generator
parities and removes the gradient mixed term. Squaring the six displacement
lengths d/2 gives the stated common leading logarithm for C_x and C_q,
without asserting equality of their finite parts.

The local orbit symplectic term remains `-P0 x qdot`: the generators are
horizontal, so their cross product pairs only with the vertical vorticity
in the Euler orbit form. There is no omitted axial-gradient contribution
to that quadratic symplectic pairing. The conjugate radius equation is

```
(h-C_x partial_z²)x=-P0 qdot.
```

Its Fourier solution gives exactly `I(k)=P0²/(h+C_x k²)` for the stated
leading-gradient action. Consequently the optical k² coefficient contains
both `h C_q` and `K C_x`. The differential field normalization changes the
field at order two and agrees with the physical angle at k=0; that map is
explicit, not hidden inside a renamed stiffness. The free-relative factor
ten follows only in its specified singular limit K/h=9 and C_x/C_q=1.

For the affine cage, eliminating the complete momentum vector first gives
`M(k)=M0-k² M0 C_p M0+...`. The constant map has
`s=(a+b)/a`, `t=-b/a`; its off-diagonal kinetic entry is
`s(at+b)=0`, and `theta_1-theta_2=s(Psi-beta)`. The pure-Psi momentum-gradient
entry is `s²(a² Cp11+2ab Cp12+b² Cp22)`. Normalizing that entry produces the
authored exact finite-core correction. This retains the cage reaction rather
than fixing its conjugate momentum incorrectly. Its positivity follows from
the two Gram matrices and the positive zero-gradient Hessian.

For an independent arithmetic check of the singular ratio, let
`a=11I/10`, `b=-9I/10`. Then

```
(a²+b²)/a = 101I/55,
(C_p)11/(C_theta)11 = 1/(9rho² Gamma² d^4),
K I = 81rho² Gamma² d^4/4.
```

It follows without inserting a target ratio that

```
C_effective/C_twist-only
 = 1 + 101 K I/(495rho² Gamma² d^4)
 = 1 + 909/220
 = 1129/220.
```

The script's matrix calculation independently yields this same result.
Finite-core coefficients remain the exact integrals and are not forced to
have the asymptotic ratio. Terms with an affine cage spatial gradient carry
the additional macroscopic displacement derivative; the pure-Psi correction
does not and therefore contributes at optical order k².

## Scope, small-ratio prescriptions and promotion

There are no blocking findings or requested corrections. The small-ratio
issues are resolved analytically at this boundary: the complete second
variation is present, symmetry directions are identified, the relevant
positive coefficient is a strict Gram norm, and an explicit remainder bound
controls the leading logarithm. No numerical soft eigenvalue, energy
subtraction, unmeasured quadrature floor or fitted finite radius supplies a
sign. Numerical eigenpair residuals and mesh/box extrapolation are therefore
not additional requirements for these exact statements.

- Verification: exact analytic proof with symbolic-verified identities and
  matrix normalization; no numerical evidence is promoted to exactness.
- Review: audited, scientific acceptance recommended at the stated scope.
- Compatibility: compatible smooth-Euler affine extension, with fixed
  circulation and explicit renormalized energy convention.
- Epistemic: the reviewed gradient and normalization routes are established;
  the complete parent P251 objective remains active.

Canonical APIs, tests and any affected coefficient consumers belong to the
parent's promotion/replay transaction. This review does not claim that a
later continuum API already uses the corrected coefficient or that release
promotion has occurred. Correction check: not needed.

## Content-addressed review boundary

The reviewed SHA-256 identifiers are:

```
3d63f9543deb6af7b2b4ce2f88e5341e1a0bb561c1d24aa191c27cd9da3d59ee README.md
6af2fee5b35b041eaab6f812315829d247cd10d8ef582abfa7bbe7278f3aac07 full-biot-savart-twist.md
3856f1b43505fd3978cc456261ca9501743be0e1f8e58b6ca90ff5fe9ede70cd radial-momentum-gradient.md
38cf68247beb6e2a7ef5e58b2d7c5fe4949400e2e7d4f9116b25c64fb65e2eae full_twist_energy.py
b88708873c2ac3f3d805c911b4df64ddbd1afe8e2f6c7917837946ae0fcd4d3a radial_gradient_elimination.py
63b38e5ef6d9913def76c253307d8860d094269b79dfc09bfc860fbf92fe2a21 cage_gradient_matrix.py
```
