# Exact flat precursor and the finite-torus transfer problem

## Result at this attempt boundary

Candidate B now has an exact `C-infinity` physical compacton in the straight
cross-section model, including a nonzero nondegenerate constant-curl core, a
fully explicit logarithmic edge law, and smooth physical pressure.  Replacing
the straight invariant direction by a finite circular direction is **not** an
exact substitution: the cylindrical curvature residual is nonzero and cannot
be absorbed by the two fixed stream-label functions on circular level sets.
The finite Euclidean carrier therefore remains open at one precise degenerate
free-boundary inverse stated below.

Candidate A remains open in the smooth category, but an analytic version of
the proposed core-to-collar matching is excluded.  Its minimum repair is a
genuinely nonanalytic `C-infinity` mixed matching theorem; analytic
Cauchy--Kowalevskaya continuation cannot furnish it.

These are route verdicts only.  Neither result changes the independent fixed
Beltrami-background plus 0250 response route or the parent objective.

## 1. An exact flat logarithmic edge

Normalize the outer stream value to `s_*=1`.  For inward distance `d>0`, set

\[
 q(d)=e^{-1/d},\qquad q(0)=0.                           \tag{17}
\]

The zero extension is `C-infinity`, and direct differentiation gives

\[
 \frac{q''}{q}=d^{-4}-2d^{-3}.                          \tag{18}
\]

Thus the `p=4` proposal is not only a leading asymptotic.  It has an exact
one-dimensional label law.  For a radial straight cross-section of outer
radius `a`, put `d=a-rho` in an outer collar.  Then

\[
 \frac{\Delta_2q}{q}
 =L^4-2L^3-\frac{L^2}{a-1/L},
 \qquad L=\log(1/q)=\frac1d.                            \tag{19}
\]

The last term is the exact polar-curvature cost.  It is lower order than
`L^4`; no physical derivative has been dropped.

Let `lambda>0`.  Define on the small-label interval

\[
 g(s)=s\left[\lambda^2+L^4-2L^3
                  -\frac{L^2}{a-1/L}\right],
 \qquad L=\log(1/s).                                    \tag{20}
\]

Then `g(s)~sL^4>0`.  If `G(s)=int_0^s g(t)dt`,

\[
 G(s)\sim\tfrac12s^2L^4,
 \qquad \int_0 \frac{ds}{\sqrt{G(s)}}
       \sim \sqrt2\int^\infty\frac{dL}{L^2}<\infty.    \tag{21}
\]

Equation (21) is the exact compact-support-principle threshold.  For the
general `sL^p` law it is finite precisely for `p>2`, agreeing with the frozen
edge exponent `alpha=2/(p-2)`.

## 2. Completing the straight nonzero core without making labels smooth

Choose `rho_0>0` before the first zero of `J_0(lambda rho)`.  On
`0<=rho<=rho_0`, take

\[
 \Phi(\rho)=cJ_0(\lambda\rho),\qquad c>0.               \tag{22}
\]

Choose `rho_0<rho_1<a`.  On `rho_1<=rho<a`, take a positive multiple of
`exp[-1/(a-rho)]`, with the multiple small enough that the endpoint values
remain strictly decreasing.  A standard smooth monotone splice on
`[rho_0,rho_1]` gives a positive `C-infinity` function which is strictly
decreasing for `0<rho<a`, agrees with (22) and (17) on open endpoint
neighborhoods, and is flat at `rho=a`.  This interpolation is a design freedom,
not a field equation.

Because `Phi` is strictly decreasing away from its center, define for
`0<s<Phi(0)`

\[
 g(s)=\lambda^2s+Phi''(\rho(s))
                   +\frac1{\rho(s)}\Phi'(\rho(s)),      \tag{23}
\]

where `rho(s)=Phi^{-1}(s)`.  Equation (23) is smooth for every positive label.
It vanishes identically on the inner stream interval by Bessel's equation and
agrees with (20) near zero.  At label zero it is continuous but not Lipschitz.

In straight Cartesian coordinates `(x,y,z)`, with
`rho=sqrt(x^2+z^2)` and no `y` dependence, define

\[
 v=(-\Phi_z,\lambda\Phi,\Phi_x),\qquad
 \mathcal P=G(\Phi)-\tfrac12|v|^2.                     \tag{24}
\]

Equations (23)--(24) give

\[
 -\Delta_2\Phi=\lambda^2\Phi-g(\Phi),                  \tag{25}
\]

and the direct Grad--Shafranov calculation gives stationary Euler.  On the
inner interval `g=0`, `curl v=lambda v`.  At the center

\[
 v(0)=(0,\lambda c,0),\qquad
 D^2\Phi(0)=-\frac{\lambda^2c}{2}I_2,                  \tag{26}
\]

so the core is nonzero and nondegenerate.

The abstract derivative `g'(s)` diverges at zero, but this is not physical
nonsmoothness.  In the outer collar, `g(Phi)=lambda^2Phi+Delta_2Phi` is a flat
`C-infinity` spatial function.  Moreover
`G(Phi)~Phi^2 log(1/Phi)^4/2`; it and all its spatial derivatives are flat.
Consequently the zero extensions of `v`, `curl v`, and the full pressure
`mathcal P` are `C-infinity` and carry no boundary sheet or traction.

This is an exact compact **cross-section**, not a compact Euclidean 3D field:
it is invariant in `y`.  Periodizing `y` would change topology and is not a
solution of the Euclidean obligation.

## 3. Exact finite-radius equation and why a circular bend fails

For a ring radius `R`, write `r=R+x`, set `psi_R=R phi_R`, and prescribe

\[
 F(\psi)=\lambda\psi,\qquad
 B_R'(\psi)=\frac1R g(\psi/R).                          \tag{27}
\]

The physical Bernoulli composition has no hidden `R` factor:
`B_R(R phi)=G(phi)`.  The exact finite-radius equation is

\[
 -\Delta_R^*\phi_R
 =\lambda^2\phi_R-\frac{r^2}{R^2}g(\phi_R),
 \qquad
 \Delta_R^*=\partial_x^2+\partial_z^2
                -\frac1{R+x}\partial_x.               \tag{28}
\]

Suppose one simply inserts the straight radial profile
`phi_0=Phi(rho)`.  Since

\[
 \Delta_R^*\Phi
 =\Phi''+\frac{R}{(R+x)\rho}\Phi',                     \tag{29}
\]

the exact residual of (28), using
`g(Phi)=lambda^2Phi+Delta_2Phi`, is

\[
 \mathcal E_R
 =\frac{x}{(R+x)\rho}\Phi'
  +\left(\frac{2x}{R}+\frac{x^2}{R^2}\right)g(\Phi).  \tag{30}
\]

It is `O(R^{-1})` on a fixed cross-section but is not zero.  More strongly,
the first term contains `x/(R+x)`.  On a nontrivial circular level its third
`x` derivative is `6R/(R+x)^4`, whereas a correction by `F F'` and
`r^2B'` with fixed stream labels is affine in `r^2` on that level.  Therefore
no redefinition of the two label values on the same circular levels absorbs
(30).  A deformed finite-radius profile and free boundary are load-bearing.

## 4. Removing abstract label singularity from the free-boundary equation

In the exact logarithmic outer collar write

\[
 \phi_R=e^{-1/T_R},\qquad T_R>0.                        \tag{31}
\]

Only `phi_R` is zero-extended; `T_R` is an interior defining function.  The
identity

\[
 \frac{\Delta_R^*e^{-1/T}}{e^{-1/T}}
 =\frac{(1-2T)|\nabla T|^2+T^2\Delta_R^*T}{T^4}        \tag{32}
\]

turns (28) into

\[
 T_R^2\Delta_R^*T_R+(1-2T_R)|\nabla T_R|^2
 =\frac{r^2}{R^2}\mathcal G(T_R)-\lambda^2T_R^4,       \tag{33}
\]

where

\[
 \mathcal G(T)=T^4\frac{g(e^{-1/T})}{e^{-1/T}}
 =1-2T-\frac{T^2}{a-T}+\lambda^2T^4                   \tag{34}
\]

on the explicit outer interval.  Unlike `g'(s)`, (34) extends smoothly to
`T=0`, with `mathcal G(0)=1`.  At the free boundary, (33) gives the exact
eikonal datum

\[
 T_R=0,\qquad |\nabla T_R|=r/R.                        \tag{35}
\]

Thus the non-Lipschitz stream-label profile is compatible with smooth physical
fields: it becomes a smooth right side coupled to a degenerate second-order
operator in a nonflat defining function.  This representation is the concrete
advance over the one-dimensional heuristic.  It does not prove solvability of
(33).

The minimum repair for Candidate B is now precise: prove, for at least one
sufficiently large finite `R`, a disk-like domain
`D_R compactly_contained_in {r>0}` and a positive solution of (28) which
converges to the straight profile strongly enough to preserve (26), while its
outer defining function solves (33)--(35).  The proof must establish the
Fredholm/invertibility statement for the degenerate linearization after fixing
the `z`-translation gauge and allowing the radial center/free-boundary graph to
move.  It must also connect the outer logarithmic interval to the exact inner
`g=0` interval and expose one profile parameter whose exact return-map
flux-action derivative is nonzero.  Standard fixed-domain Dirichlet existence
or the compact support principle alone does not provide these topology, core,
twist, or inverse rows.

If that theorem is supplied, (27) immediately gives a compact Euclidean
stationary Euler field.  Equations (31)--(34) make velocity, vorticity and
pressure flat at the boundary; (26) persists by `C^2` closeness; choosing
`R` larger than the support radius gives axis separation.  Exact disjoint
assembly and the positive density formulas in the frozen README then apply
with no far-field or interaction term.  The 0250 action/current transfer is
still a separate consumer.

For clarity, the one-template physical costs would then be computed directly
from the same finite-radius solution, not inherited from the straight limit:

\[
 V_{\rm supp}=2\pi\int_{D_R}r\,dr\,dz,\qquad
 E_{\rm one}=\rho\pi R^2\int_{D_R}
       \frac{|\nabla\phi_R|^2+\lambda^2\phi_R^2}{r}\,dr\,dz,               \tag{37}
\]

\[
 p_R=G(\phi_R)-\frac{R^2}{2r^2}
          (|\nabla\phi_R|^2+\lambda^2\phi_R^2).                            \tag{38}
\]

For `N` separated copies in a cell, `E_vol=N E_one/V_cell` and
`phi_supp=N V_supp/V_cell`.  The constant material density is still `rho` on
the entire cell, including the quiescent exterior.  A claimed measured tube
density additionally needs the finite solution's positive tube integral
`J_T` and nonzero twist; neither is manufactured by (37).

## 5. Candidate A: an analytic subroute is excluded

Let

\[
 Q=\frac{|\nabla\psi|^2+F(\psi)^2}{r^2}-A(\psi).       \tag{39}
\]

If the fields and label functions in a connected neighborhood from the outer
collar through the nondegenerate core are real analytic, then `Q` is analytic.
Vanishing on the open localizable collar forces `Q=0` throughout that connected
neighborhood.  The pressure is then a stream function there, hence
`u dot grad p=0` at the core.  The reviewed 0252 theorem forces the central
swirl to vanish, contradicting the retained nonzero core.

This refutes only the analytic collar-to-core continuation.  A smooth solution
can change equations through a flat, nonanalytic stream-label transition, so
Candidate A is not refuted.  Its minimum repair is an actual `C-infinity`
mixed-system construction that matches all physical jets, maintains a
pressure-value interval reserved for `supp chi'`, and proves the speed first
integral only there.  An analytic local existence theorem is inapplicable.

## Route verdicts and licenses

`candidate_A_route_verdict: blocked with missing genuinely nonanalytic
C-infinity mixed core/collar matching; analytic continuation is refuted by
identity theorem plus the 0252 core theorem`.

`candidate_B_route_verdict: blocked with missing finite-radius degenerate
free-boundary inverse; exact C-infinity straight compacton, inner constant-curl
core, physical pressure, p=4 threshold, transformed boundary equation, and
non-absorbable circular-bend residual are established`.

`evidence_scope: exact analytic precursor and representation-scoped transfer
obstruction; no compact Euclidean field, density supplier, response channel,
parent completion, or exhaustion conclusion`.

`parent supplier license`: a future theorem solving (28) or (33)--(35) with
the stated convergence licenses exactly one compact stationary Euler template,
smooth pressure, a nonzero nondegenerate constant-curl core, axis separation,
disjoint-copy assembly, and the density normalization of README equations
(15)--(16), provided the same theorem also supplies the stated nonzero
flux-action twist.  It does not license the 0250 response/action/current join.
