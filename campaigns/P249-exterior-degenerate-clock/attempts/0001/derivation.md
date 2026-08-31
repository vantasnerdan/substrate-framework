# Attempt 0001 exact derivation

This attempt separates the compact clock license from the existence of a
stationary clock. The first is an exact group/action result. The second is
tested first on the smallest commuting core-split profile, where algebra and
calculus decide the route before a numerical solver is eligible.

## Field, action, and compact subgroup

Use the mostly-plus metric `eta = diag(-1,1,1,1)` and the conditional completed
P239 density

\[
  \mathcal L=-\frac12\langle F_{\mu\nu},F^{\mu\nu}\rangle_h
  -V(M)-\kappa\eta^{\mu\nu}q(\partial_\mu P_t,\partial_\nu P_t),
\]

where `X=eta^{-1}M`, `P_t(X)` is the simple timelike spectral projector,
`h^{-1}=eta^{-1}-2P_t eta^{-1}`, and

\[
  F_{\mu\nu}=(\partial_\mu M)\eta^{-1}(\partial_\nu M)
  -(\partial_\nu M)\eta^{-1}(\partial_\mu M).
\]

Let `A_23=-1`, `A_32=1`, with all other entries zero, and
`R(q)=exp(qA)`. Then `R(q)^T eta R(q)=eta`, `R(q+2pi)=R(q)`, and

\[
  \Phi_q(M)=R(q)M R(q)^T
\]

is a smooth global SO(2) action on the complete symmetric-field space. Because
`R` commutes with `eta`, `X`, `P_t`, `h^{-1}`, every derivative, and every
curvature transform by the matching similarity/congruence laws. Trace
potentials, the double two-form contraction, and the projector-current
bilinear are therefore invariant. No field-dependent or tapered generator is
used.

For the exterior vacuum

\[
 M_\infty=\operatorname{diag}(-g,1,0,0),\qquad g>0,
\]

the last two eigenvalues coincide and `Phi_q(M_infinity)=M_infinity`
pointwise. The infinitesimal generator is

\[
  \zeta(M)=AM+MA^T,
\]

which vanishes on `M_infinity` and is nonzero wherever the 2-3 degeneracy is
lifted or the core leaves the fixed set.

## Noether current and Legendre map

Let `Pi^mu` be the variational derivative of the full density with respect to
`partial_mu M`, regarded as a cotangent on symmetric matrices with the trace
pairing. Localizing the exact global transformation to `q=q(x)` gives

\[
  j^\mu=\langle\Pi^\mu,\zeta(M)\rangle.
\]

Global invariance gives

\[
 \left\langle\frac{\partial\mathcal L}{\partial M},\zeta\right\rangle
 +\langle\Pi^\mu,D\zeta[\partial_\mu M]\rangle=0.
\]

Combining this identity with the Euler-Lagrange equation yields
`partial_mu j^mu=0`. This is the exact Noether bridge; a finite quadratic norm
alone is not being reinterpreted as charge.

On a relative equilibrium `M(t,x)=Phi_{q(t)}(M_0(x))`, curvature and projector
time derivatives are linear in `omega=dot q`, while the action is independent
of `q`. The parity-even P239 action therefore reduces exactly to

\[
  L[M_0,\omega]=-E_{\rm stat}[M_0]+C[M_0]\omega^2,
  \quad J=2C\omega,
  \quad E_J=E_{\rm stat}+\frac{J^2}{4C}.
\]

The Legendre map is nondegenerate precisely where `0<C<infinity`. A smooth
Gaussian core split supplies an explicit nonempty finite-positive-inertia
sector; stationarity is a later O3 obligation.

## Smallest commuting core-split route

Take

\[
 M_0=\operatorname{diag}(-g,1,d(x),-d(x)),\qquad d(x)\to0.
\]

At `beta=c=1` the projected M5.17 potential is exactly

\[
  V(d)=3d^2+4d^4.
\]

All static derivatives are proportional to the same tangent-split matrix, so
their eta-commutators vanish. The clock/spatial curvature is nonzero and the
exact reduced coefficient is

\[
  C[d]=32\int d^2|\nabla d|^2\,d^3x.
\]

Writing `u=d^2` gives the especially transparent functional

\[
 C[u]=8\int|\nabla u|^2,\qquad
 E_J[u]=\int(3u+4u^2)+\frac{J^2}{32\int|\nabla u|^2}.
\]

For `u>0`, its exact fixed-J Euler equation is

\[
 16\omega^2\Delta u+3+8u=0,
 \qquad \omega=\frac{J}{16\int|\nabla u|^2}.
\]

The width/amplitude question is decided by a collapse path. For any nonzero
smooth nonnegative shape `f`, put `u_{a,R}(x)=a f(x/R)` and define positive
shape integrals `A1=int f`, `A2=int f^2`, `K=int |grad f|^2`. Then

\[
 E_J(a,R)=R^3(3aA_1+4a^2A_2)+\frac{J^2}{32a^2RK}.
\]

Along `R=1/a`, every term tends to zero as `a` tends to infinity:

\[
 E_J(a,a^{-1})=\frac{3A_1}{a^2}
 +\frac{4A_2}{a}+\frac{J^2}{32aK}\longrightarrow0.
\]

For nonzero `J`, every finite admissible profile has strictly positive energy,
so the infimum zero is not attained. The commuting one-profile split cannot be
the requested stationary fixed-width minimum. This route verdict does not
touch noncommuting spatial cores, bundle/gauged representations, or minimal
kinetic/topological completions.
