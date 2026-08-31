# Attempt 0004 exact auxiliary-amplitude seed

## Scope boundary

This attempt tests the smallest canonical vacuum-vanishing phase field before
any coupled M5 solver is eligible.  It deliberately does **not** identify a
decoupled nonlinear Klein--Gordon soliton with the requested M5 clock.  The
load-bearing question here is whether the auxiliary seed itself is exact and
well posed, so that the next route has one sharply named missing construction:
a local Lorentz-covariant same-action bridge to the M5 core split together with
the complete cross-sector variations.

Use the mostly-plus metric and a complex Lorentz scalar `psi` with

\[
 {cal L}_\psi=-\frac12\eta^{\mu\nu}
 \partial_\mu\bar\psi\,\partial_\nu\psi-W(|\psi|),
 \qquad
 W(f)=3f^2-6f^4+4f^6.
\]

The action `psi -> exp(iq) psi`, `q in R/(2 pi Z)`, is smooth, faithful away
from the vacuum, and fixes `psi=0` pointwise.  Its current is

\[
 j^\mu=\operatorname{Im}(\bar\psi\,\partial^\mu\psi).
\]

For `psi(t,x)=exp(i Omega t) f(x)`, with real `f`, define

\[
 N[f]=\int f^2,\qquad Q=\Omega N,
\]

so the Legendre map is regular for every nonzero finite profile.  The exact
fixed-charge functional is

\[
 E_Q[f]=J_0[f]+\frac{Q^2}{2N[f]},\qquad
 J_0[f]=\int\left(\frac12|\nabla f|^2+W(f)\right).
\]

## Exact existence and stability bridge

The primary mathematical source is Benci and Fortunato,
arXiv:0903.3508, Theorem 13 and its hypotheses `(W-i)`--`(W-iiii)`.
Every hypothesis is discharged by exact polynomial calculus:

\[
 W(f)=f^2\left[4\left(f^2-\frac34\right)^2+\frac34\right]\geq0,
 \quad W(0)=W'(0)=0,\quad W''(0)=m^2=6.
\]

Writing `W(f)=m^2 f^2/2+N_{nl}(f)` gives

\[
 N_{nl}(f)=-6f^4+4f^6.
\]

At `f_0=sqrt(3)/2`, `N_nl(f_0)=-27/16<0`; and
`N_nl'(f)=24f^3(f^2-1)`, so `f_1=2>f_0` has positive derivative and satisfies
the alternative growth condition `(W-iiii)(b)`.  The theorem therefore gives
a nonempty open charge set of hylomorphic standing waves and orbital stability
in the scalar phase space.

The fixed-frequency existence interval is also exact.  Its lower edge obeys

\[
 \Omega_0^2=\inf_{f>0}\frac{2W(f)}{f^2}=\frac32,
 \qquad m^2=6.
\]

Thus `Omega=2` lies strictly inside `Omega_0<Omega<m`.  At that frequency

\[
 W(f)-\frac12\Omega^2f^2
 =f^2(1-6f^2+4f^4),
\]

which is negative for
`(3-sqrt(5))/4 < f^2 < (3+sqrt(5))/4`.  This is an exact hylomorphy witness,
not a fitted frequency or numerical shooting result.

## Complete fixed-charge calculus

For any real variation `eta`, with `Omega=Q/N[f]`,

\[
 \delta E_Q
 =\delta J_0-\frac{\Omega^2}{2}\delta N,
\]

and hence a critical profile satisfies

\[
 -\Delta f+W'(f)-\Omega^2f=0.
\]

The complete second variation is

\[
 \delta^2 E_Q
 =\delta^2J_0-\frac{\Omega^2}{2}\delta^2N
 +\frac{\Omega^2}{N}(\delta N)^2.
\]

Equivalently,

\[
 \delta^2E_Q[\eta,\eta]
 =\int\left(|\nabla\eta|^2+[W''(f)-\Omega^2]\eta^2\right)
 +\frac{4\Omega^2}{N}\left(\int f\eta\right)^2.
\]

The last positive rank-one term is required by the Legendre reduction.  The
envelope identity gives `dE/dQ=Q/N=Omega` on a minimizing branch.

At `Omega=2`, the radial equation and its linear far field are

\[
 f''+\frac2r f'=2f-24f^3+24f^5,
 \qquad f(r)\sim A\frac{e^{-\sqrt2 r}}r.
\]

The tail is natural and unsupported.  Derrick scaling of the exact
fixed-charge functional gives

\[
 E_Q[f(x/R)]=R T+R^3V+R^{-3}E_{\rm rot},
 \quad T+3V-3E_{\rm rot}=0,
\]

so the canonical amplitude route has competing scale powers rather than the
commuting M5 collapse of attempts 0001/0003.

## Route verdict

The auxiliary scalar seed is established exactly: compact normalized phase,
positive Hamiltonian, regular charge-frequency map, a nonempty open set of
orbitally stable localized standing waves, complete fixed-charge Hessian, and
an action-derived exponential tail.  The campaign route is nevertheless
**blocked at the M5 bridge**, not established at O3.  This action fragment does
not dynamically split the M5 tangent eigenvalues and its scalar theorem says
nothing about cross-sector criticality or stability.  Those are the named
constructions for the next attempt.
