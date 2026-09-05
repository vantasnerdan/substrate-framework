# 0261 — admissible radial nondegeneracy and finite-radius inverse

Owner: `herdr geometry-review pane w3:p4`, continuing the independent 0253/0255
finite-radius geometry route as a fresh bounded construction attempt.  Write
scope is attempt 0261 only.  Base checkpoint is `47b1d25`; the parent is the
unchanged issue200/P251 objective.  No prior attempt, source module, central
proposal entry, registry, release, generated file, or commit is changed here.

This contract is frozen before opening any new source body or carrying out new
spectral/free-boundary analysis.  The coordinator must add the central 0261
entry and a successful schema receipt before substantive work begins.  The
0253, 0255, 0256, and 0258 artifacts are accepted branch inputs only at their
stated scopes; in particular, a heuristic in 0255 is not promoted here merely
by repetition.

## Four frozen levels

The parent objective remains one stationary smooth constant-density Euler
ensemble supporting actual coupled histories, inherited physical action and
current, and compact Euclidean nonzero finite-core tubes at usable positive
measured density.  The periodic supplier and the independent 0250 response
integration remain active.  This attempt cannot reduce the parent to stationary
existence or geometry.

The active obligation is an actual finite-`R` solution of the axisymmetric
Grad--Shafranov equation with compact **velocity**, a single axis-separated
disk-like meridional positivity set, a physically `C-infinity` flat edge, and
`C^4` convergence on one fixed regular inner annulus to the nonzero
constant-curl Bessel core.  The two new load-bearing parts are:

1. construct a translation-balanced admissible straight profile whose radial
   linearized operator has no zero eigenvalue, while preserving the exact
   core, edge, support radius, and balance constraints; and
2. construct the compatible nonlinear defining-function/free-boundary inverse,
   including the degenerate collar, domain graph, kernels, cokernel, gauges,
   and matching to the interior.

The concrete attempt compares two radial-nondegeneracy routes before using the
winner in the nonlinear inverse:

1. a bordered transversality construction using the explicit eigenvalue
   variation supplied by the coordinator;
2. a constructive profile/Sturm--Liouville route using exact ODE segments and
   a Wronskian, comparison, or Evans-function sign, if the bordered
   independence cannot be proved on the admissible constraint manifold.

A verdict attaches only to the route actually established or refuted.  Failure
of either radial route does not refute all compactons or finite-radius Euler
carriers, and no exhaustion conclusion is proposed.

## Inherited exact finite-radius problem

Use `r=R+x`, `epsilon=R^{-1}`, `psi_R=R phi_R`,
`F(psi)=lambda psi`, and `B_R'(psi)=R^{-1}g(psi/R)`.  The exact equation is

\[
 \mathcal E_\epsilon(\phi)=
 -\Delta\phi+\frac{\epsilon}{1+\epsilon x}\phi_x
 -\lambda^2\phi+(1+\epsilon x)^2g(\phi)=0.             \tag{1}
\]

For a radial straight profile `Phi(rho)`, 0256 gives the first curvature
forcing

\[
 f_1=\Phi_x+2xg(\Phi)                                  \tag{2}
\]

and its translation pairing

\[
 Q(\Phi)=\int_{\mathbb R^2}\bigl(\Phi_x^2-
                  \lambda^2\Phi^2\bigr)\,dx\,dz.       \tag{3}
\]

The fixed-support profile family in 0256 brackets both signs of `Q`, and
quadraticity gives a simple balanced root in its family.  This licenses the
translation scalar, not radial nondegeneracy or a nonlinear inverse.

At `Q(Phi)=0`, the first `m=1` correction recorded in 0255 has the exact form

\[
 h=\Phi'(\rho)v(\rho)\cos\theta,
 \qquad
 (\rho\Phi'(\rho)^2v'(\rho))'
 =\rho\Phi'(\rho)(\Phi'(\rho)+2\rho g(\Phi(\rho))).    \tag{4}
\]

The zero total integral in (3) is precisely the compatibility condition making
`v'` bounded at the flat edge; for `Phi\sim C e^{-1/d}`, its required leading
edge behavior is `v'\to a`.  Equation (4) is an explicit first-order shape
correction only.  It does not supply the full inverse or nonlinear solution.

## Frozen admissible profile manifold

Fix radii

\[
 0<a_{\rm core}<a_{\rm buf}<a_{\rm edge}<a,            \tag{5}
\]

before choosing perturbations.  On `[0,a_core]`, require the same exact
constant-curl Bessel core

\[
 \Phi=cJ_0(\lambda\rho),\qquad c>0,                    \tag{6}
\]

including its core value and all matching jets.  On `[a_edge,a]`, require the
same logarithmic flat-edge law and the same support radius.  Perturbations may
be supported only in the declared transition interval
`I=(a_core,a_edge)`; no later argument may shrink the retained core or move
the flat edge silently.  A Bessel subinterval inside `I` may be used as a
sacrificial spectral buffer only if its existence is part of the constructed
balanced profile and it is disjoint from the frozen core.

Write

\[
 y=-\Phi'>0,qquad
 \Phi(\rho)=\int_\rho^a y(s)\,ds.                      \tag{7}
\]

For a multiplicative transition perturbation

\[
 y_t=y e^{t h},\qquad h\in C_c^\infty(I),              \tag{8}
\]

preserving the full inner profile is equivalent, after the unchanged outer
edge is fixed, to the exact scalar constraint

\[
 A(t,h):=\int_{a_{\rm core}}^{a_{\rm edge}}
             y(e^{th}-1)\,d\rho=0.                    \tag{9}
\]

The exact translation constraint is `Q(Phi_{t,h})=0`.  At a balanced base
profile its linearizations are

\[
 A_1(h)=\int_I yh\,d\rho,                              \tag{10}
\]

\[
 B_1(h)=2\pi\int_I yh\left[
     \rho y-2\lambda^2\int_0^\rho s\Phi(s)\,ds
                         \right]d\rho,                 \tag{11}
\]

up to an irrelevant common nonzero normalization of `B_1`.  All trial
directions must also preserve positivity/monotonicity, the signed transition
needed by 0256, and inverse-design regularity of `g`; these are open conditions
only after the exact equalities (9) and (3) are solved.

## Candidate R1 — explicit bordered transversality

For the radial linearized operator write

\[
 L_0=-\partial_\rho^2-\rho^{-1}\partial_\rho+V,
 \qquad
 V=\frac{y''}{y}+\frac{y'}{\rho y}-\frac1{\rho^2}.      \tag{12}
\]

Assume provisionally that `L_0 z=0` has a regular radial solution satisfying
the actual flat-edge domain condition, and put `q=z/y` and
`M=\rho y^2`.  Then the supplied exact identities are

\[
 (Mq')'=-\frac{M}{\rho^2}q,                            \tag{13}
\]

\[
 \delta V=h''+\frac{M'}M h',                           \tag{14}
\]

and, with the correct radial measure and boundary terms verified,

\[
 C_1(h)=\int_I2M\left[(q')^2-\frac{q^2}{\rho^2}\right]
                       h\,d\rho                       \tag{15}
\]

is the zero-eigenvalue quadratic-form derivative.  Candidate R1 must prove,
for an actually constructed balanced admissible profile, that there exist
three declared transition directions `h_0,h_1,h_2` for which

\[
 \det\begin{pmatrix}
 A_1(h_0)&A_1(h_1)&A_1(h_2)\\
 B_1(h_0)&B_1(h_1)&B_1(h_2)\\
 C_1(h_0)&C_1(h_1)&C_1(h_2)
 \end{pmatrix}\ne0.                                   \tag{16}
\]

It must then apply the finite-dimensional implicit-function theorem to the
**exact** constraints (9) and `Q=0`, leaving a tangent direction on which the
simple radial eigenvalue moves.  Continuity and spectral isolation must turn
that derivative into an admissible nearby profile with `0` outside the radial
spectrum.  If the nullspace has dimension greater than one, the full crossing
matrix replaces the scalar (15).

Equation (16) is the central independence obligation.  A genericity slogan,
three arbitrary bumps, or an affine/nonaffine weight heuristic without an
actual admissible base profile is insufficient.  The Bessel-buffer expansion
suggested in 0255 may prove (16) only after showing simultaneously that the
buffer is outside the immutable core, the relevant radial zero has the stated
Wronskian representation there, the three physical weights are not affinely
dependent, and the profile still satisfies the signed transition and exact
balance.

## Candidate R2 — constructive profile and spectral exclusion

If R1's determinant vanishes or cannot be made compatible with the exact
constraints, build `y` and `Phi` in explicit ODE segments:

1. the unchanged Bessel core on `[0,a_core]`;
2. a finite collection of smooth signed-transition segments furnishing an
   independent balance border;
3. the unchanged logarithmic edge on `[a_edge,a]`, with all physical spatial
   jets flat at `a`.

Inverse-design `g` from the resulting strictly decreasing `Phi`.  Propagate
the regular radial solution of `L_0z=0` from the center and the admissible
edge solution from the flat boundary.  Prove their Wronskian is nonzero at a
fixed matching radius, or establish an equivalent Sturm-comparison or
Evans-function sign, while a separate segment parameter solves `Q=0`.
The proof must retain the same exact core and edge and show that smoothing the
segments preserves both the spectral sign and the exact scalar balance after
retuning its declared border.

A monotone-`V` shortcut is not presumed compatible with the signed `g`
transition required by (3).  Candidate R2 succeeds only with an actual
profile and a rigorous spectral exclusion, not with a statement that such a
profile should be generic.

## Compatible nonlinear defining-function/free-boundary inverse

For the logarithmic collar write

\[
 \phi_R=e^{-1/T_R}.
\]

The exact regular physical equation becomes

\[
 T_R^2\Delta_R^*T_R+(1-2T_R)|\nabla T_R|^2
 =\frac{r^2}{R^2}\mathcal G(T_R)-\lambda^2T_R^4,       \tag{17}
\]

with the boundary jet

\[
 T_R=0,\qquad |\nabla T_R|=r/R
       \quad\hbox{on }\partial D_R.                    \tag{18}
\]

The positivity domain is an unknown graph near the straight disk.  Candidate
I must fix it by an explicit domain map, use an even-in-`z` or barycentric
gauge for vertical translation, use the balanced-profile parameter for the
broken radial translation row, and incorporate the R1/R2 radial
nondegeneracy.  The linear construction must account for:

- the exact `m=1` translation mode and the quadrature correction (4);
- the radial `m=0` block, with zero excluded by R1 or R2;
- all `m>=2` blocks, including the factorization/coercivity and their actual
  edge-domain conditions;
- the boundary graph and collar equation, whose principal coefficient
  `T^2` vanishes at `T=0` and therefore is not covered by an ordinary uniformly
  elliptic implicit-function theorem;
- matching of collar `T` and interior `phi`, with enough derivatives to make
  velocity, vorticity, and pressure spatially `C-infinity` after zero
  extension.

Success requires a stated Banach or tame scale on a fixed reference disk, a
bounded right inverse for the full bordered linearization, compatibility of
its collar and interior traces, and nonlinear composition estimates for
`e^{-1/T}` and the non-Lipschitz label law.  The estimates must be strong
enough at one finite large `R` to close a contraction, Nash--Moser, or other
fully stated inverse theorem.  An `O(R^{-1})` residual, modewise formal
solutions, or an inverse that projects away the free-boundary condition does
not establish finite-radius existence.

The completed solution must satisfy

\[
 D_R\Subset\{(r,z):r>0\},
\]

with `D_R` one disk, and its exact physical fields and pressure must be

\[
 u=-\frac{\psi_z}{r}e_r+\frac{\lambda\psi}{r}e_\theta
       +\frac{\psi_r}{r}e_z,
 \qquad
 p=G(\phi)-\frac{R^2}{2r^2}
       (|\nabla\phi|^2+\lambda^2\phi^2).               \tag{19}
\]

Distributional substitution across the flat edge, axis separation, nonzero
core velocity, and the pressure normalization are part of the construction,
not downstream conventions.

## Retained twist and physical density costs

Attempt 0258 supplies, for the exact inner profile
`Phi=cJ_0(lambda rho)`,

\[
 q=\frac{J_1(\lambda\rho)}{\rho J_0(\lambda\rho)},
 \qquad J_u=c\rho J_1(\lambda\rho),
 \qquad \left.\frac{dq}{dJ_u}\right|_{\rm core}
       =\frac{\lambda^2}{8c},                           \tag{20}
\]

and the finite-`R` normalized contour formula

\[
 q_R(s)=\frac{2\pi}{\lambda s
       \int_{\{\phi_R=s\}}(R/r)\,d\ell/|\nabla\phi_R|},
 \quad
 J_{u,R}(s)=\frac{\lambda}{2\pi}
       \int_{\{\phi_R>s\}}(R/r)\phi_R\,dr\,dz.        \tag{21}
\]

Consequently no extra twist-tuning parameter is demanded: `C^4` convergence
on one fixed regular inner annulus is the precise conditional input under
which the nonzero local twist persists by contour quadrature.  This source
scope does not prove the finite-`R` solution, density, response, action, or
current.

For any successful solution the one-copy costs remain

\[
 V_{\rm supp}=2\pi\int_{D_R}r\,dr\,dz,
 \qquad
 E_{\rm one}=\rho_{\rm mat}\pi R^2\int_{D_R}
 \frac{|\nabla\phi_R|^2+\lambda^2\phi_R^2}{r}\,dr\,dz. \tag{22}
\]

For disjoint copies, measured support density and energy density are
`N V_supp/V_cell` and `N E_one/V_cell`; the constant material density occupies
the whole cell.  Positive measured tube density still requires a positive
tube observable on this same solution and exact later assembly.  Copy counting
does not supply the 0250 field-changing response or the physical
action/current join.

## Frozen comparison, claim ladder, and oracle

Compare R1 and R2 by: preservation of the exact core/edge and balance;
explicit proof of radial zero exclusion; independence and economy of profile
borders; compatibility with the full degenerate inverse; physical smoothness;
and robustness adequate for `C^4` core convergence.  R1 is tried first because
it uses the exact variation already derived.  R2 is a genuine constructive
spectral continuation, not a restatement of transversality.

| Step | Positive intent | Pass licenses | Does not license | Failure scope | Unlocks |
| --- | --- | --- | --- | --- | --- |
| 1 | Admissible balanced radial nondegeneracy | One profile preserving the exact core, edge and `Q=0`, with `0` excluded from the radial spectrum | Any finite-`R` solution | R1 bordered family or R2 segmented spectral family | Full bordered linear inverse |
| 2 | Degenerate collar/interior right inverse | A right inverse with all modes, gauges, boundary graph and matching controlled | A nonlinear solution before estimates close | Declared defining-function atlas and spaces | Nonlinear inverse theorem |
| 3 | Actual finite-`R` compact carrier | One classical solution of (1), one disk, axis separation, flat zero extension, (19), and `C^4` core convergence | 0250 response/action/current or parent completion | Constructed large-`R` branch | Individual geometry review and same-field integration |

`requires`: activated 0253 compacton equations; 0256 exact translation balance
at its stated scope; 0255 first-correction identities after direct checking;
0258 twist only as the conditional consequence (20)--(21); and any primary
degenerate/free-boundary theorem opened only after schema activation.

`pass_licenses`: at maximum, one actual finite-radius, axis-separated,
physically smooth compact-velocity stationary Euler template retaining the
specified core and its local flux-action twist, with exact pressure and
one-copy costs.

`does_not_license`: an extra twist parameter; radial genericity without
constraint independence; ordinary elliptic inversion of (17); a translating
far field; compact vorticity in place of compact velocity; response gain,
action/current normalization, actual coupled histories, arbitrary topology,
stability, parent completion, promotion, or scientific exhaustion.

`maximum_verdict`: exact analytic radial nondegeneracy plus a complete
finite-`R` defining-function/free-boundary existence proof with all geometry,
regularity, core, pressure, and normalization rows above.

`failure_scope`: the explicit R1 bordered profile manifold, the explicit R2
segmented spectral construction, or the declared nonlinear inverse spaces.
Each negative result names its mechanism and activates the next route; none is
an obligation-level no-go.

The strongest oracle is analytic: an explicit nonzero constraint determinant
or Wronskian/Sturm sign, exact nonlinear constraint restoration, a proved
bordered right inverse in the degenerate collar/interior spaces, exact PDE and
pressure substitution, a distributional zero-extension check, disk/axis
geometry, and quantitative `C^4` core convergence.  No numerical remainder is
frozen.  Therefore the small-ratio numerical skill does not bind and no
numerical production is licensed.  Any future spectral computation requires
an append-only verifier contract with the relevant eigenvalue residual,
zero-mode gauge, discretization floor, observed-order extrapolation, and sign
robustness prescriptions before it is run.

`route_verdict: pending central/schema activation and substantive analytic
construction`.

`evidence_scope: preregistered constrained radial transversality, constructive
spectral alternative, compatible degenerate nonlinear inverse, retained 0258
twist scope, physical costs, and parent boundary only`.

Next coordination achievement: `w3:p1` adds the central 0261 manifest entry and
a schema-valid receipt.  Until then this pane opens no new source body and
performs no new determinant, spectral, collar, or nonlinear-inverse analysis.
