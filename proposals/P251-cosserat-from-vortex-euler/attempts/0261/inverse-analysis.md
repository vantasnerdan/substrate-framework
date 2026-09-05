# The finite-radius inverse in defining-function variables

This file tests the failure-derived Keldysh/Hörmander representation and
connects it to the constrained radial construction.  It establishes the exact
linearization, its physical weighted realization, its bracket structure, the
free-boundary trace mechanism, and the finite-dimensional nonlinear reduction.
It also identifies the one estimate not supplied by the transferred theorems.

## 1. A global defining function and exact conjugacy

The amplitude of the constructed profile may be chosen with
`0<Phi<1`.  Define throughout its positivity disk

\[
 T_0=\frac1{\log(1/\Phi)}.                              \tag{1}
\]

This is smooth and positive at the Bessel center.  Choose the normalization of
the exact edge so that `Phi=exp(-1/d)` on an outer collar; then `T_0=d` there.
Thus no artificial interior/collar interface is required at the base profile.

For `epsilon=R^{-1}`, put

\[
 a_\epsilon(x)=(1+\epsilon x)^2,
 \quad
 \Delta_\epsilon^*=\Delta-\frac{\epsilon}{1+\epsilon x}\partial_x,
 \quad
 \phi=e^{-1/T}.                                        \tag{2}
\]

The physical equation and its defining-function form are

\[
 E_\epsilon(\phi)=-\Delta_\epsilon^*\phi-\lambda^2\phi
                         +a_\epsilon g(\phi)=0,         \tag{3}
\]

\[
 \mathcal F_\epsilon(T)=T^2\Delta_\epsilon^*T
 +(1-2T)|\nabla T|^2-a_\epsilon\mathcal G(T)
 +\lambda^2T^4=0,                                      \tag{4}
\]

where `mathcal G(T)=T^4g(e^{-1/T})/e^{-1/T}` is smooth at
`T=0`.  The identity

\[
 \boxed{\mathcal F_\epsilon(T)
       =-\frac{T^4}{\phi}E_\epsilon(\phi)}              \tag{5}
\]

is exact.  At a solution, with

\[
 A_\epsilon=-\Delta_\epsilon^*-\lambda^2
                   +a_\epsilon g'(\phi),               \tag{6}
\]

it gives the exact linear conjugacy

\[
 D\mathcal F_\epsilon(T)v
 =-\frac{T^4}{\phi}
      A_\epsilon\left(\frac{\phi}{T^2}v\right).         \tag{7}
\]

The natural physical Hilbert measure is

\[
 d\mu_\epsilon=(1+\epsilon x)^{-1}dx\,dz,              \tag{8}
\]

equivalent up to the irrelevant constant `R` to `r^{-1}dr dz`.
With flat physical boundary behavior, (6) is symmetric in this measure.  At
`epsilon=0`, its edge potential is

\[
 -\lambda^2+g'(\Phi)=d^{-4}-6d^{-3}+O(d^{-2}),          \tag{9}
\]

so the physical quadratic-form realization is semibounded with compact
resolvent.  The radial construction removes its `m=0` zero; the 0255
factorization leaves precisely the two `m=1` translations in the full space,
or only the radial translation after restricting to functions even in `z`.

Equation (7) is the required mode-gluing bridge.  It ties the degenerate
defining-function operator to the same physical spectral domain used in the
radial proof; a separate unweighted Dirichlet spectrum is not imported.

## 2. Exact edge operator and sign

On the straight circular collar let `rho=a-d`.  There

\[
 \mathcal G(T)=1-2T-\frac{T^2}{a-T}+\lambda^2T^4.      \tag{10}
\]

Linearizing (4) at `epsilon=0,T_0=d` gives

\[
 \mathscr L v=d^2v_{dd}
 +\left(2(1-2d)-\frac{d^2}{a-d}\right)v_d
 +\frac{d^2}{(a-d)^2}(v_{\theta\theta}+v).             \tag{11}
\]

The transverse sign is therefore

\[
 \mathscr L=2\partial_d+d^2\Delta+O(d\partial_d)+O(d^2),
                                                               \tag{12}
\]

where `d` increases into the physical support.  For Fourier mode `m`,

\[
 \mathscr L_m v_m=d^2v_m''+\left(2-4d-
              \frac{d^2}{a-d}\right)v_m'
       +\frac{d^2(1-m^2)}{(a-d)^2}v_m.                 \tag{13}
\]

For `m=1`, the moderate constant solution is the defining-function image of a
translation.  The second homogeneous derivative is proportional to

\[
 e^{2/d}d^4/(a-d),                                     \tag{14}
\]

and is excluded by the physical weighted domain.  This reproduces, rather
than assumes, the one-sided regularity selection at the edge.

The first curvature derivative of (4) at fixed `T=d` is

\[
 \left.\partial_\epsilon\mathcal F_\epsilon(d)
 \right|_{0}
 =\left[d^2-2(a-d)\mathcal G(d)\right]\cos\theta.      \tag{15}
\]

Thus `L_T tau = -dF/depsilon` gives at `d=0`

\[
 2\tau_d(0,\theta)=2a\cos\theta,                       \tag{16}
\]

which is exactly the 0255 result `v'(a)=a` after its sign convention is
translated.  The Keldysh representation therefore agrees with both the
physical curvature equation and the boundary eikonal; there is no sign
discrepancy.

## 3. Hörmander structure, with its exact scope

At a general nearby solution of (4), set

\[
 X_1=T\partial_x,\qquad X_2=T\partial_z,               \tag{17}
\]

\[
 X_0=(2-5T)\nabla T\mathbin\cdot\nabla
   -\frac{\epsilon T^2}{1+\epsilon x}\partial_x.       \tag{18}
\]

Then

\[
 D\mathcal F_\epsilon(T)=X_1^2+X_2^2+X_0+c_T          \tag{19}
\]

for a smooth zeroth-order coefficient `c_T`.  On `T=0`, equations (4) and
`mathcal G(0)=1` give

\[
 |\nabla T|^2=a_\epsilon>0.                            \tag{20}
\]

Although `X_1=X_2=0` there, the drift brackets satisfy

\[
 [X_0,X_j]|_{T=0}=2|\nabla T|^2\partial_j,
 \qquad j=x,z.                                         \tag{21}
\]

They span the tangent plane uniformly for `epsilon` small.  At the straight
edge this is equivalently

\[
 X_1=d\partial_d,\qquad X_2=\frac d{a-d}\partial_\theta,
 \quad X_0=\left(2-5d-\frac{d^2}{a-d}\right)\partial_d,
                                                               \tag{22}
\]

with `[X_0,X_1]=2\partial_d` and
`[X_0,X_2]=(2/a)\partial_\theta` at `d=0`.

This verifies the proposed step-two drift bracket.  Hörmander's theorem and
the local drift-Schauder estimates cited in `source-transfer.md` therefore
license local hypoellipticity and intrinsic estimates for a **two-sided smooth
extension** of (19).  They do not supply a one-sided free-boundary trace map,
a global Fredholm inverse, or a smoothly parameterized right inverse.

## 4. Why the boundary graph cannot be discarded

A fixed circular boundary would require the Eulerian perturbation `v` to have
zero boundary trace.  This is not a harmless gauge.  Under (7), a moderate
physical solution `u` determines

\[
 v=T_0^2u/\Phi,                                        \tag{23}
\]

and its trace contains one independently selected coefficient in every
Fourier mode.  Only the two `m=1` coefficients can be changed by adding the
translation kernel.  Requiring all other trace coefficients to vanish would
impose infinitely many extra conditions on the source.

The first-order toy reduction exposes the mechanism.  Dropping lower-order
and `d^2v_{dd}` terms from (13) gives

\[
 2v_m'-\frac{m^2d^2}{a^2}v_m=F_m.                     \tag{24}
\]

Matching to an interior value at `d=delta` selects

\[
 v_m(0)=e^{-m^2\delta^3/(6a^2)}v_m(\delta)
 -\frac12\int_0^\delta e^{-m^2s^3/(6a^2)}F_m(s)\,ds.  \tag{25}
\]

Thus the trace is selected by interior regularity and source matching, and is
generically nonzero.  Formula (25) also displays the characteristic tangential
scale `d\sim |m|^{-2/3}`.  It is diagnostic only; the omitted normal second
derivative and lower coefficients prevent importing its kernel as the full
Poisson operator.

The full normal second derivative can still be retained in an exact scalar
Volterra representation.  Divide (13) by `d^2` and use

\[
 \mu(d)=(a-d)d^{-4}e^{-2/d}.                            \tag{26}
\]

Then every moderate Fourier coefficient obeys

\[
 (\mu v_m')'+\mu\frac{1-m^2}{(a-d)^2}v_m
       =\frac{\mu}{d^2}F_m,                            \tag{27}
\]

\[
 \mu(d)v_m'(d)=\int_0^d\mu(s)\left[
   \frac{F_m(s)}{s^2}+
   \frac{m^2-1}{(a-s)^2}v_m(s)\right]ds.              \tag{28}
\]

The omitted integration constant is exactly the explosive branch (14), so the
physical quadratic-form domain selects (28).  Laplace localization of
`mu(s)/mu(d)` immediately recovers `2v_m'(0)=F_m(0)`.
Equations (27)--(28) establish the correct one-mode edge selection with the
full operator.  They do not establish an `m`-uniform trace estimate strong
enough to sum the modes and perturb to a nonradial boundary graph.  That
uniform estimate, rather than existence of each scalar moderate solution, is
the remaining analytic issue.

The compatible unknown is therefore a boundary graph `b(theta)`.  A Hanzawa
map sends its moving disk back to the reference disk, and the pulled-back
defining function is written `T=dS`.  At the linear level, the map

\[
 (\dot S,b)\longmapsto v=d\dot S+\chi(d)b(\theta)      \tag{29}
\]

is an isomorphism between a zero-trace interior correction plus the boundary
graph and an arbitrary smooth Eulerian `v`; `chi=1` at the edge.  The cosine
coefficient of `b` is fixed as the radial-center gauge, and evenness in `z`
removes the sine translation.  All other Fourier traces remain unknowns and
are selected by the global solve.

## 5. The exact nonlinear Lyapunov--Schmidt reduction

Let `beta` be the core-preserving transverse balance coordinate constructed in
`radial-construction.md`; its family changes the inverse-designed label
`g_beta` and satisfies

\[
 Q(\Phi_0)=0,
 \qquad \partial_\beta Q(\Phi_\beta)|_0\ne0.           \tag{30}
\]

For every `beta`, the straight `T_beta` is an exact solution at
`epsilon=0`.  Use the even-in-`z` free-boundary chart (29), remove the cosine
kernel by its center gauge, and project the target off the one translation
cokernel.  If the global tame right inverse described in the next section is
available, the range equation has a smooth solution

\[
 U=U(\epsilon,\beta),                                  \tag{31}
\]

where `U` denotes the interior defining function and all unrestricted boundary
graph modes.  The remaining scalar equation `S(epsilon,beta)=0` vanishes
identically when `epsilon=0`, so Hadamard's lemma gives

\[
 S(\epsilon,\beta)=\epsilon\widetilde S(\epsilon,\beta).
                                                               \tag{32}
\]

The exact 0256 projection gives, with a fixed nonzero normalization,

\[
 \widetilde S(0,\beta)=Q(\Phi_\beta).                  \tag{33}
\]

Equations (30) and (33) allow the ordinary scalar implicit-function theorem to
choose `beta=beta(epsilon)`.  This is the correct desingularized use of the
balance parameter: directly differentiating the unscaled equation with
respect to `beta` at `epsilon=0` would give zero because the entire straight
profile family consists of solutions.

If (31) is earned, (32)--(33) close the finite-dimensional part of the
nonlinear problem.  Positivity of `S`, smallness of the boundary graph, and
`R>a+norm(b)_infinity` then give one axis-separated disk.  The exact
`exp(-1/T)` composition makes the physical streamfunction and all its spatial
derivatives flat.  Smooth convergence in the fixed inner region implies the
`C^4` convergence needed by 0258 without an additional twist parameter.

## 6. The remaining analytic estimate

At the straight solution, (7) and the physical spectral gap give an `L^2`
range inverse for (19), modulo translation, in the measure (8).  Equations
(19)--(21) give local intrinsic regularity.  These facts do **not** yet prove
that the physical form inverse has a smooth moderate quotient (23), nor that
its boundary trace and inverse depend tamely on the domain graph, profile, and
iterate.

The missing load-bearing lemma is a one-sided collar Poisson/trace estimate for
the **full** operator (13), uniform for nearby coefficients.  In one acceptable
form it would state that every compatible smooth source and smooth interior
matching datum have a unique moderate collar solution, with

\[
 \|v\|_{k+2,X;\,0<d<\delta}
 +\|v|_{d=0}\|_{k+\sigma;\,\partial D}
 \le C_k\bigl(\|\mathscr Lv\|_{k,X}
       +\|v\|_{k_0,X;\,d=\delta}
       +\|\partial_dv\|_{k_0,X;\,d=\delta}\bigr),      \tag{34}
\]

for a declared intrinsic scale, some fixed finite loss/gain `sigma`, and tame
constants uniform under the Hanzawa/profile perturbations.  It must also prove
that exclusion of the explosive branch (14) is equivalent to the physical
quadratic-form domain and that the trace in (34) glues to the uniformly
elliptic interior inverse.  Formula (25) suggests the correct smoothing scale
but does not prove (34).

Once (34) is proved, the local drift-Schauder estimates patch it to the
interior spectral inverse; the Hanzawa map and the smooth function
`mathcal G(T)` give a smooth tame nonlinear map.  Hamilton's theorem can then
consume the smooth tame family of right inverses and validate (31).  Hamilton's
theorem does not itself provide (34), and local hypoellipticity cannot replace
it.

## Verdict and minimum repair

`inverse_route_verdict: blocked with the missing construction named`.

The defining-function representation is compatible and materially stronger
than an ordinary elliptic-IFT slogan: its exact sign, lower-order coefficients,
physical conjugacy, bracket generation, free trace, mode gluing, translation
border, and scalar nonlinear reduction all check.  It does not yet establish
an actual finite-`R` solution because no transferred theorem proves the global
one-sided tame estimate (34) for this quadratically degenerate drift operator.

The minimum repair is not another existence citation.  It is a proof of (34),
preferably by a Fourier/Volterra construction that retains the full
`d^2v_{dd}` term and then a perturbative intrinsic estimate for the nonradial
Hanzawa coefficients.  With that lemma, the already established spectral gap,
free boundary graph, and (30)--(33) provide the complete nonlinear inverse.
Without it, claiming Nash--Moser or a compact Euclidean Euler carrier would be
an unsupported extension.
