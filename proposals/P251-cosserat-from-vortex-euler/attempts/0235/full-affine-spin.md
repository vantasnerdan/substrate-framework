# All three acoustic rotations and the complete finite-cut spin

This uses the actual C015/C016 first laboratory spatial jet on its fixed
curl u=-u cell, in its unit curl-length and density-rho convention. It
does not infer a full vector spin from the normal-plane detector. Let
u=(psi,v), psi=cos Y+a cos Z, a=1/100, v=J grad psi, and T=v.grad.
All ambient pressure is retained. The passive C016 returns have the
symmetric input factor kappa_Y V_X+kappa_X V_Y and consequently vanish
in EACH antisymmetric affine contraction, not just the X-axis one.

## Exact full-axis Euler columns

For unit Cartesian e, write A_e q=e cross q and R_e=A_e x. The constant
curl identity and the actual periodic pressure projection give

    grad u_e=partial_e u-A_e u,
    P(A_e u)=partial_e u,
    w_R=[u,R_e]=A_e u-(R_e.grad)u,
    L w_R=0,    L partial_e u=0,
    L R_e=w_R-2partial_e u.                            (1)

Here L=-P(u.grad+Du), and the affine terms are actual laboratory-Bloch
derivatives, not periodic rotations of a fixed cell. In particular
partial_X u=0, so 0231's special axial cancellation is included.

The complete C015 API antisymmetric contraction
prepared(e_j,e_k)-prepared(e_k,e_j), for cyclic (e,j,k), has zero
symmetric d/current and zero complete material rate. Its raw lift is
q_e=A_e u-partial_e u=-grad u_e, and its actual negative-helicity
return is +partial_e u. Their complete first velocity is A_e u;
the return vanishes only on X. Thus the actual D column is the rigid
stationary rotation xi_D=R_e,w_D=w_R for ALL three axes. The complete
C016 z remainder has forcing -partial_e u; together with the raw lift
it gives the actual first-affine solution

    w_e(t)=R_e+t(w_R-2partial_e u),
    p_e(t)=2rho u_e-t(R_e.grad)p-2t partial_e p.         (2)

Direct Euler substitution proves (2); its extra pressure and velocity
terms are not discarded. This is a finite-dimensional neutral Euler
response, not an unresolved spectral integral. Its material motion is
nevertheless not t R_e. With W_t=exp(-t ad_u), the actual Lin solution
with xi_e(0)=0 is

    xi_e(t)=t R_e+eta_e(t),
    eta_e(t)=2t e-2 integral_0^t W_s e ds.              (3)

Indeed (partial_t+ad_u)eta_e=-2t partial_e u. The initial tag variation
is zero for the V column. The D column is the genuine rotation of the
entire transported tag. These facts fix both phase preparations and
their physical material observations; no homogeneous Lin term is chosen
after seeing a desired spin history.

## A positive finite axial-cut material tag

Let chi(psi) be the fixed positive normal tag profile of 0231 and choose
a smooth nonnegative axial profile alpha(X_0), of finite support and
positive integral L_alpha, centered at zero. This is an actual finite
material tag, not a wall. Its reference transport is

    g_t(X_0,q_0)=(X_0+t psi(q_0),F_t(q_0)),
    chi_t=chi(psi) alpha(X-t psi),                     (4)

where F_t is the normal steady-flow map. Write bar psi for the
chi-weighted mean, M=rho L_alpha integral chi, and

    Sigma=rho L_alpha integral chi(psi-bar psi)^2>0.

The strict inequality follows because a smooth nonzero positive tag
has support of nonzero area on which this nonconstant analytic psi
cannot be constant. The tag centroid travels with bar psi e_X; its
normal centroid stays zero by the actual reflection symmetries. Its
central mass trace obeys exactly

    tr I(t)=tr I(0)+Sigma t².                          (5)

The normal marginal and its second moment are stationary, while the
axial material cut shears by t psi. Equations (4)--(5) retain the
complete material boundary; a reset Eulerian axial cell would miss
this term. The continuous ambient remains the complement of these
actual tags in the hybrid observation and pressure-bond balance of
0232. Point-fluid and finite-cut currents are distinct, with their
exact cut-cell or finite-radius comparison retained.

## Exact physical spin and its nonaffine coefficient

For an actual xi,w, use the full moment variation from 0232,

    delta S=rho integral chi_0 [
       (xi-delta X) cross (u-V)
       +r cross(D_t xi-delta V)].                     (6)

Integrals can equivalently be evaluated on the transported tag. Let
S_e(t) be (6) for the rate column (2)--(3), and let
s(t)=sum_e e.dot S_e(t). The whole-field isotropic axial response is
s(t) Omega_dot/3; it is the FULL three-axis trace. The projected
detector would instead use only its body-X row.

At t=0, xi=0, w=R_e and p_e=2rho u_e, giving

    s(0)=2 tr I(0),       s'(0)=0,
    s''(0)=2rho integral chi_0 |u-V|²
                     -2 integral chi_0 r.grad p
           =d_t² tr I(0)=2Sigma>0.                   (7)

Here is the independent full-boundary derivation of the last row. The
material torque variation is

    delta tau=-integral_D [(xi-delta X) cross grad p
                     +r cross(grad p_e+(Hess p)xi)].

Differentiate at zero using D_t xi=R_e, delta V=0 and (2). The two
rigid rotation terms cancel in the e component. The sum over axes
annihilates the symmetric Hessian pressure derivative and uses

    sum_e e.dot[(u-V) cross grad u_e]=(u-V).curl u,
    sum_e e.dot[r cross((Hess u_e)u)]=r.(u.grad)curl u.

Since curl u=-u and rho(u.grad)u=-grad p, (7) follows. This calculation
includes pressure, moved normals and lever arms; it is not a guessed
rigid geometric inertia. The initial first derivative also vanishes
by the fixed centered axial profile and invariant normal measure.

The D column rotates the actual reference spin vector. Its isotropic
mean static spin is zero, as expected, but that does not remove (7).
Time reversal leaves Sigma and this even-time curvature unchanged.
The stationary angle-null b return of 0231 adds only a constant to
s(t). It can normalize s(0), but cannot cancel s''(0)>0 on these actual
finite-cut tags. This is a refutation of that one proposed *complete-
current repair*, not of 0231's projected result or the parent objective.

## Explicit full memory rather than an unspecified remainder

The trace can be written in normal-flow correlation functions, with
the exact finite-cut shear separated. For a fixed scalar weight w(psi),
define

    C_w(t)=rho L_alpha integral w(psi(q_0))
                         F_t(q_0).q_0 d²q_0,
    d(psi)=(psi-bar psi)chi'(psi).

These are actual defining integrals; d is a derived response weight,
not a negative material tag. The exact result is

    s(t)=2 tr I(t)+2[C_chi(t)-C_chi(0)]
                +2[t C_d'(t)-C_d(t)+C_d(0)].          (8)

To derive it, the flow Jacobian is
Fcal_t=[[1,t grad psi(q_0)^T],[0,DF_t]]. In (3),
W_s e(g_t(a))=Fcal_t Fcal_(t-s)^(-1)e. For the matrix
E with columns eta_e, define ax E_i=epsilon_ijk E_jk. Summing (6)
for its eta part gives

    rho integral [r.D_t ax E-(u-V).ax E].

The normal part of ax E is -2 integral_0^t(t-s)v(F_s)ds.
The normal-normal antisymmetric Jacobian part satisfies, by integration
by parts in its source coordinate and invariance of psi,

    integral w ax(DF_tau)_X
         =-d_tau integral w'(psi) F_tau.q_0.

Applying this with w=chi(psi)(psi-bar psi), and adding the normal
part, gives (8). Fixed center terms cancel in (6). The correlations
are even functions of time by invariance of their real scalar weight.
Their second derivatives are
C_w''(0)=-rho L_alpha integral w|v|². The independent divergence
identity

    integral [(chi(psi)(psi-bar psi))']|grad psi|²
          =integral chi(psi)(psi-bar psi)psi

uses Delta psi=-psi and the compact normal support; it reduces the
second derivative of (8) exactly to (7). All correlations and their
finite time derivatives are finite on the fixed smooth tag. No
acoustic-time or constant-clock approximation was made.

## Failure-generated actual angle-null transverse controls

The same Euler cell has an exact passive axial family

    w=(exp(-tT)g,0,0),
    xi=(exp(-tT)(h+t g),0,0).                          (9)

Here h,g are smooth fixed normal profiles. Its pressure is zero and
the full Euler/Lin equations hold. It leaves the integrated normal tag
density and its axial angle unchanged, but moves the actual axial cuts.
For mean-free odd profiles the centroid terms vanish. Define the
vector P_f(t)=integral chi q exp(-tT)f. Direct substitution in (6) gives

    delta S_perp=rho L_alpha J[
                    P_h'(t)+t P_g'(t)-P_g(t)],
    delta S_X=0.                                     (10)

The required identity is P_f'=integral chi v exp(-tT)f, proved by
transport integration by parts. Thus (10) is a concrete physical
transverse-current control, rather than an assigned time-varying
coefficient multiplying stationary b. Its corresponding orbital,
shape and ambient-current changes and its full initial action must
accompany any use in the 0227 interface. The sign-changing response
profiles do not change the fixed positive reference tag.

Equations (8) and (10) now name the exact kernel-matching problem:
match the actual three-axis trace and its current with sources h,g,
including the finite-cut shear polynomial, while preserving the actual
macro U and inherited full forms. An arbitrary dense-range assertion
without accuracy-versus-norm control is not substituted for that
construction. A different current representative may instead place
the transverse axial-cut moments into an explicitly derived orbital/
face current, as allowed by 0232, but it cannot simply discard them.

The proved response and its failure-generated physical control are
durable progress on both registered routes. Full kernel matching or
the complete alternative orbital-current transformation is the next
construction; neither is asserted by the stationary axial repair.
