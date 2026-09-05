# Perturbation, physical gluing, and finite-radius closure

## 1. What glues automatically once collar regularity is known

At the straight profile the physical operator

\[
 A_0=-\Delta-\lambda^2+g'(\Phi)                       \tag{1}
\]

has the 0261 form-domain range inverse modulo the one even translation mode.
For a smooth defining source `F`, the conjugate physical source is

\[
 f=-\frac{\Phi}{d^4}F.                               \tag{2}
\]

It is flat at the support edge.  Solve the global projected form problem for
`u`, and define `v=d^2u/Phi` on `d>0`.  The branch calculation in
`trace-estimate.md` excludes the explosive local solution and proves that
`v` is precisely the moderate solution selected by the one-sided exit
formula.  There is no separate conormal matching condition: `u` is one
global form solution, so its collar and uniformly elliptic interior
restrictions already have the same Dirichlet and conormal traces.

This repairs a possible overdetermination in the preregistration.  For an
independently solved collar problem, inner Dirichlet data determine the inner
conormal value.  For the global problem, both values are inherited from the
same physical solution.  They are never prescribed independently.

On `d>=delta_0/2`, ordinary elliptic estimates control the low inner trace in
0263 equation (26).  The exact radial collar estimate and the physical
spectral gap therefore yield a global smooth right inverse at `epsilon=0`,
modulo the translation projection, in the intrinsic projective-limit space.

## 2. The no-loss perturbation estimate

After flattening a nearby Hanzawa graph and using `T` as the normal variable,
the linearization has the same structural form

\[
 P_q=\sum a_q^{ij}X_iX_j+X_{0,q}+c_q,                 \tag{3}
\]

with `a_q` positive in the diffusion directions, normal drift bounded below,
and coefficients close to the radial ones.  Its diffusion has the exact
one-sided representation (9) of `trace-estimate.md`, uniformly in a small
smooth neighborhood.  The invariant-half-collar/two-sided-resolvent
construction there proves:

- the characteristic edge remains entrance-only;
- the inner exit is reached with a uniform exponential moment;
- the actual outer trace is selected from the source and inner exit datum;
- the solution and trace obey a uniform intrinsic estimate on every shrunken
  collar;
- differentiating the extended equation gives smooth tame dependence on the
  Hanzawa graph, profile, and iterate; and
- no arbitrary two-sided boundary value is introduced.

Use fixed flattened model fields and write

\[
 P_q=P_0+K_q,\qquad
 R_q=(I+R_0K_q)^{-1}R_0,                              \tag{4}
\]

where `R_0` is the global radial right inverse modulo translation.  Estimate
(27) of `trace-estimate.md` controls every `X_iX_jv`, `X_0v`, and selected
trace in the same intrinsic graph scale.  Hence `K_q` maps
`mathcal X^{s+2}` to `mathcal Y^s`, with the tame product estimate

\[
 \|K_qv\|_{\mathcal Y^s}
 \le C_s\|q\|_{s_0}\|v\|_{\mathcal X^{s+2}}
 +C_s\|q\|_{s+2}\|v\|_{\mathcal X^{s_0+2}}.         \tag{5}
\]

Smallness in the fixed low norm makes (4) converge there; the higher estimates
and the resolvent identity make the inverse family smooth tame.  The
finite-dimensional translation kernel and cokernel persist as smooth spectral
projections.  The cosine boundary-graph coefficient is the center gauge and
evenness removes the vertical translation, exactly as in 0261.

The global construction is performed on the physical form solution and then
conjugated.  Consequently the collar and interior restrictions are one
solution and conormal matching is automatic.  There is no independently
specified Cauchy datum at the interface and no boundary condition at the
characteristic edge.

## 3. Downstream finite-radius closure

The 0261 scalar reduction is

\[
 S(\epsilon,\beta)=\epsilon\widetilde S(\epsilon,\beta),
 \qquad \widetilde S(0,\beta)=C Q(\Phi_\beta),        \tag{6}
\]

and the frozen profile family has a simple zero of `Q`.  The smooth tame range
inverse above supplies `U(epsilon,beta)`.  Hadamard division by `epsilon` and
the ordinary scalar implicit-function theorem then choose
`beta=beta(epsilon)` for all sufficiently small positive `epsilon=1/R`.

For such an `R`, the small Hanzawa graph has a single disk-like positivity set
strictly inside `|x|<R`; revolution gives an axis-separated solid torus.  Put

\[
 \psi_R=R\phi_R,\qquad F(\psi)=\lambda\psi,
 \qquad B_R'(\psi)=R^{-1}g_{\beta(1/R)}(\psi/R).      \tag{7}
\]

Then the exact Grad--Shafranov equation holds and the full pressure is

\[
 p_R=G_\beta(\phi_R)-\frac{R^2}{2(R+x)^2}
       (|\nabla\phi_R|^2+\lambda^2\phi_R^2).          \tag{8}
\]

Because `phi_R=exp(-1/T_R)` and all polynomial `T_R^{-1}` losses are dominated
by that exponential, velocity, vorticity, and (8) have smooth zero extensions;
there is no vortex sheet or pressure traction.  Smooth convergence on the
unchanged inner Bessel neighborhood gives the `C^4` hypothesis of 0258, so
the same finite ring has nonzero normalized flux-action twist on a fixed
positive-width regular annulus without another parameter.

`perturbation_route_verdict: established as stated in the intrinsic
projective-limit scale`.

`finite_R_candidate_I_verdict: established as stated, consuming the 0256
balance, 0261 constrained radial nondegeneracy, and 0258 twist implication at
their recorded scopes`.

This establishes the stationary compact Euclidean carrier and its retained
twist.  It does not establish the 0250 same-field response/action/current
join or actual coupled histories.  Those remain independent parent
dependencies.

## 4. Same-field stationary density and normalization cost

Let `A_core` be the fixed positive-width regular meridional annulus used by
0258.  On the finite ring its physical core volume is

\[
 V_{\rm core}=2\pi\int_{A_{\rm core}}(R+x)\,dx\,dz>0. \tag{9}
\]

Enclose the whole compact support in a ball of radius `D` and choose a lattice
spacing `L>2D`.  Translated and rigidly rotated copies centered on that lattice
have disjoint supports for every orientation.  Because each velocity and its
pressure are flat and zero off its own support, the locally finite sums

\[
 u=\sum_j u_j,\qquad p=\sum_j p_j                 \tag{10}
\]

solve the same constant-density stationary Euler equations pointwise: at each
point at most one summand and its derivatives are nonzero, so there is no
discarded Biot--Savart, convective, or pressure cross term.

Occupy each lattice site independently with probability `p_occ>0` and apply a
uniform random lattice shift; optional independent rigid rotations give an
isotropic orientation law without changing stationarity.  The exact mean
number density, measured core-volume fraction, and energy density are

\[
 n=\frac{p_{\rm occ}}{L^3},\qquad
 f_{\rm core}=\frac{p_{\rm occ}V_{\rm core}}{L^3}>0,
 \qquad
 e=\frac{p_{\rm occ}E_1}{L^3}.                       \tag{11}
\]

Here `E_1` is the finite kinetic energy of one template.  Thus no amplitude,
pressure, or long-range interaction renormalization is spent to obtain usable
positive measured stationary tube density.  The density can be varied through
`p_occ` and `L` below the explicit nonoverlap packing bound.  Equation (11)
does not create the field-changing response or the physical action/current;
those must be evaluated on this same ensemble by the independent parent
construction.
