# Fixed finite generalized-force-free profile: an actual acoustic theorem

This is the registered stronger continuation, not an extrapolation of the
ordered delta/k^2 estimate. The planar geometry and its entire-cell
Arnold operator are those of 0161. Choose a sufficiently large but FIXED
finite C, then let k approach zero. All averages use the whole fluid.

## 1. Profile, exact frame and state space

In laboratory coordinates use W_lab=sqrt(C+lambda^2 psi^2) and the
stationary generalized-force-free u=(v,W_lab) already constructed.
Hereafter use the exactly specified axial mean frame, and write

    W=W_lab-<W_lab>,  <W>=0,
    c(psi)=W'/zeta'=-W'/lambda^2,
    a(psi)=1-W'^2/lambda^2,
    delta=||W||_{C^m}=O(C^(-1/2)).

This changes neither the physical background nor its pressure; laboratory
frequencies retain the known axial Doppler shift k<W_lab>. The full state
is

    eta=curl_h w_h, r=b-c eta,
    w_h=m+K eta+i k d b, b=c eta+r.

At k=0 the zero-mean group is exactly the direct combination

    eta_t=L eta,  r_t=-A r,
    L=-A H, H=I-lambda^2 G, A=v.grad.

Its boundedness in ||eta||_2+||r||_2 follows from 0161 and unitary
transport. The k-dependent full Euler/pressure operator is an analytic
bounded O(k) perturbation in this fixed space. No transport inverse on
all L2, and no ambient spectral gap, is assumed.

## 2. Full first-order translation source

Put T_X=-X.grad zeta=lambda^2 X.grad psi. The complete zero-k velocity
translation is (K T_X,c T_X)=(-Dv X,-X.grad W). In the Hodge coordinates,
its first-order horizontal-vorticity source is -i F_X, where

    F_X=lambda^2 grad psi.grad G(X.grad W)
         +lambda^2(W-psi W')X.grad psi
         -W' X.grad(|v|^2/2).                           (1)

The exact r equation is

    r_t+A r=-i k S,
    S=W r+pi-c(grad W cross w_h)_z+c zeta b.

On the full translation its leading source reduces to

    S_T=a A(v.X).                                      (2)

In deriving (2), both the transported position/axial translation and the
full pressure pi_T=-X.grad p are retained. Consequently the r corrector
is explicitly

    r1_X=-i a v.X,   -A r1_X-i S_T=0.                  (3)

At W'=0 this is the planar pressure-response corrector. At the exact
constant-curl slope W'=-lambda it is zero, explaining why that previous
candidate has a different acoustic sign.

## 3. The conserved-row drift is exactly the axial mean

Let ell be 0161's bounded conserved translation projection, ell T=I.
It is essential to check its value on F instead of declaring a frame.
The result for an uncentered profile is

    ell F=<W> I.                                       (4)

Here is a whole-cell boundary proof. Write

    F_X=curl_h V_F,
    V_F=W T_h X+(A+Dv)grad G(T_z X).

Its complete mean vanishes, since the mean of its second term is
<v T_z X>, cancelling <W T_h X>. On a straight separatrix edge,
psi=psi_s, zeta=zeta_s and W=W_s are constant. Use

    (A+Dv)grad phi=grad(A phi)+zeta J grad phi,
    phi=G T_z X=-X.grad G W,
    Delta G W=-W+<W>.

Integrating along a complete periodic edge of length l and normal n gives

    integral T_h X.dot t=-zeta_s(X.n)l,
    integral partial_n phi=(W_s-<W>)(X.n)l,
    integral V_F.dot t=-zeta_s <W>(X.n)l.                (5)

Tangential total derivatives integrate to zero. Thus V_F-<W>T_h has
zero period on each fundamental edge. On each bounded polygon its curl
integral also vanishes: W_s,zeta_s are constant, and the boundary
integrals of T_h and the volume integral of T_z vanish. Replacing x by
its polygon-centered R therefore does not change the curl moment.
Integration by parts over the fundamental parallelogram, retaining its
edge periods, now proves (4). In the chosen mean frame ell F=0.

The exact polynomial calculations in fixed_profile_source.py provide
independent exposing examples, including profiles with nonzero mean:
they are not a substitute for (5) or an assumption that every drift is
the mean. No residual O(k) axial-advection term has been discarded.

## 4. Actual transport primitives on every polygon

Each hexagonal cell is invariant under C6 about its center; each of the
two triangular cells is invariant under C3 about its own center. The
affine C3 rotations preserve the GLOBAL periodic psi, not only a local
quadratic jet. In reciprocal coordinates the triangle centers have
(a,b)=(+/-2pi/3,+/-2pi/3); the three phases there agree modulo 2pi.
The periodic Green operator and H commute with those affine rotations.

There is only one critical point in each polygon interior, its center.
Solving sin a=sin b=-sin(a+b) gives that fact directly; all other critical
points lie on the common saddle network. Every regular contour in a
polygon is consequently a single invariant loop preserved by that
polygon's C3 or C6 rotation. A vector-covariant source has zero orbital
average on EACH such loop, since neither C3 nor C6 fixes a nonzero
two-dimensional vector. This applies to the column family F_X and to
the adjoint sources below. Global oddness alone would not suffice on
the two different triangular cells.

The transport inverse used here is only the zero-orbital-mean primitive
of those particular smooth sources. Near the nondegenerate saddle level,
the period obeys T(s)<=C(1+|log|s-psi_s||). Direct integration along a
loop and subtraction of its mean gives

    |A^{-1}g|<=2T(s)||g||_infinity,
    ||A^{-1}g||_2^2<=C integral T(s)^3 ds <infinity.       (6)

Away from saddles, periods are bounded. These primitives are elements
of D(A)={f in L2: A f in L2 in the distributional sense}; they need
not be globally H1 or bounded across separatrices. No derivative across
a jump is omitted: the weak transport domain uses zero normal flux.

By (4) and skew-adjointness,
<t_i,A^{-1}F>=-<A^{-1}t_i,F>=0. The whole-cell inverse of H on its
first-shell complement is bounded in every Sobolev order. Hence solve

    H eta1_X=-i A^{-1}F_X,
    L eta1_X=i F_X,   ell eta1_X=0,                      (7)

adding the unique finite translation needed for the last condition.
This is legitimate even though eta1 is initially only L2: A G is
bounded, so (7) also puts eta1 in D(A)=D(L). Its norm is O(delta),
with constants determined by the fixed polygon geometry and the exact
Arnold complement, not by k.

## 5. The exact adjoint current and the computed stiffness

The first-shell Green identity now gives the full physical mean

    m_t=-i k[<g eta>+<v r>]+k^2<W d b>,
    g=J grad G W+c v.                                   (8)

In particular <g T_X>=0 is the exact complete-translation flux
cancellation. The vector function g is smooth, in the vector sector,
orthogonal to ker H, and O(delta). Solve the ACTUAL adjoint row

    A h=H^{-1}g,    <h,L eta>=<g eta>.                   (9)

The source H^{-1}g is smooth and vector-covariant on each polygon, so
(6) constructs this h in L2 and D(A), with ||h||_2=O(delta). There is
no omitted cohomological period or assumed bounded inverse of A.

Define the improved physical current with the translation row removed,

    p_c=m+i k[<h(eta-T_X)>+<R r>].                       (10)

Write the exact eta equation as eta_t=L eta+T_m+i k E_k, retaining
the complete pressure/Hodge formula in E_k. Equations (2),(8),(9) give

    (p_c)_t=k^2[-<h E_k>+<R S>+<W d b>].               (11)

Both zero-k fast response rows have become actual derivatives in
(10); their finite-k production terms remain in (11).

The resulting acceleration matrix is the following SAME-field integral:

    K_W=<h F>-<a v tensor v>+<W d(c T)>,
    C_W=-K_W.                                          (12)

Equivalently K_W=-i<g eta1+v r1>+<W d(cT)>, by (3),(7),(9).
This independent form identifies (12) with the complete mean of the
prepared Euler tangent, not merely a formal current. No coefficient is
supplied from a target wave equation.

The matrix is real. Proper C6 rotations and the corresponding reflected
source/transport transformation make it dihedrally invariant; it is a
scalar multiple of I. In particular the reflection sign of A^{-1}
cancels the pseudovector sign of g. With W in the mean frame,

    ||F||=O(delta), ||h||=O(delta), c=O(delta),
    a-1=O(delta^2),
    C_W=C_v+O(delta^2).                                (13)

All constants in (13) are finite source and primitive estimates from
(6)-(9) at the fixed planar field. Thus for all sufficiently large
but FIXED finite C, C_W=c_W^2 I with c_W^2>0. This positivity is an
analytic bound on the full reaction, not a sampled small eigenvalue.

## 6. Smooth actual histories on the acoustic interval

The exact divergence-free slow tangent is specified in Hodge variables by

    eta=T_X+k eta1_X, r=k r1_X,
    b=c eta+r, w_h=K eta+i k d b.                        (14)

Equations (3),(7) cancel its order-k Euler residual, leaving O(k^2)
in the fixed response norm. The full k=0 group is bounded and the
remaining generator perturbation is O(k). For x=|k|X, X_t=m, the
actual remainder y=w-m-T_k X therefore satisfies the same bounded-group
Gronwall estimate as 0161, uniformly on 0<=t<=T/|k|.

To obtain SMOOTH data, approximate eta1 in the generator graph norm by
smooth periodic functions, project onto the vector sector, and subtract
its finite conserved translation rows. Smooth functions form a core
for transport by a smooth divergence-free vector field (periodic
mollification and the Friedrichs commutator identity give this directly);
L differs from -A by the bounded lambda^2 A G. Choose the approximation
so that both its L2 error and its L-residual error are at most |k|.
The sector projection commutes with L, and the translation correction
has zero L-image. Consequently the smooth version of (14) still has
O(k^2) residual, uniformly bounded response norm, and changes (12) only
by O(k). Its high Sobolev norms may depend on k; none are presumed
uniform. Hodge reconstruction is exact, including every pressure tail.

Explicitly, for a periodic mollifier J_e the commutator is

    [A,J_e]f(x)=integral [v(x)-v(y)].grad rho_e(x-y) f(y)dy.

Its L2 norm is uniformly bounded by c||Dv||_infinity||f||_2; on smooth
f it tends to zero by direct differentiation, and density extends that
limit to every L2 f. For f in D(A), A J_e f=J_e A f+[A,J_e]f therefore
converges to A f. The bounded A G perturbation gives the asserted L
graph convergence. The smoothing scale is selected for each nonzero k
before launching its actual Euler history; no uniform high-derivative
estimate for the logarithmic primitive is asserted.

Using these actual smooth preparations, (10) differs from m by O_T(k):
eta-T_X and r stay bounded after the |k|X scaling. Equation (11) then
gives, on tau=|k|t,

    x_tau=p_c+O_T(k),   (p_c)_tau=-C_W x+O_T(k).

For the unchanged actual common-V data and the smooth displacement
phase, this proves

    m(t)=cos(|k| c_W t)V0-|k|c_W sin(|k|c_W t)X0
            +O_T(|k|)(|V0|+|k X0|),
    0<=t<=T/|k|,                                       (15)

on each sufficiently large FIXED-C background. Laboratory observations
retain the mean axial Doppler shift stated in section 1. The O(k)
constants may depend on this fixed C. No fixed background is changed
while taking k->0 in (15).

## 7. Same action and physical phase normalization

Choose the smooth velocity preparation (14) and initial material
xi0=(X0,0)exp(ikz). Lin reconstruction gives the actual initial momentum
pi_D=rho(w0-T_full X0). Its horizontal cell mean is zero. The common-V
phase has xi0=0, pi0=rho(V0,0). Hence the full symplectic phase matrix
has mass rho exactly, with zero displacement/displacement block. These
are explicit full-Euler circulation data, not a claim that an arbitrary
well-prepared phase lies on one fixed coadjoint leaf.

On the actual phase family, (15) and its corrected-current derivative
make the (|k|X,p_c) chart C1-close in slow time to its invertible
oscillator matrix. Pulling the conserved full material symplectic form
and the exact moving-frame connection through that chart gives the
leading action with mass rho and stiffness rho k^2 C_W, with O(k)
normalized corrections. This is the same actual material action, not
an appended kinetic energy or a free variation of a prescribed on-shell
cell profile. The finite-k connection and explicit current (10) remain
part of the theorem.

Nonlinear amplitude is chosen small last, for each smooth preparation
and its finite T/|k| window. The fixed-C acoustic theorem supplies no
uniform nonlinear-amplitude or all-time invariant-subspace claim.
The same-field optical construction, generic spatial directions,
global closed-torus embedding and the parent's complete coupled
continuum remain their own active obligations.
