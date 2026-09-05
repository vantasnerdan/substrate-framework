# Positive actual force-free acoustic histories: a finite-C ordered family

This executes the registered failure-generated repair after the
constant-curl sign calculation. It keeps the SAME triangular planar
geometry, entire cell, density and Euler action. The new axial profile
is not claimed to have constant curl.

## 1. Exact background and explicitly specified observation frame

For C>0 set

    W_C=sqrt(C+lambda^2 psi^2),  u_C=(v,W_C),
    f_C=-lambda^2 psi/W_C.

Then curl u_C=f_C u_C exactly: its horizontal component is -W_C' v
and its vertical component is zeta=-lambda^2 psi. The Bernoulli
constant is C/2, since

    W_C^2-lambda^2 psi^2=C,
    p=-(|v|^2+lambda^2 psi^2)/2=-|u_C|^2/2+C/2.

The stationary field is smooth and periodic with finite energy density
for every finite C. There is no finite-energy limiting claim as C tends
to infinity. Observe it in the specified axial Galilean coordinates
s=z-U0 t, U0=sqrt(C). The same actual Euler histories then have base
(v,h), h=W_C-U0. Laboratory phases retain the Doppler factor
exp(ik(z-U0 t)); a new arbitrary optical clock is not being chosen.

For every fixed finite derivative order m,

    delta=||h||_{C^m} <=D_m/sqrt(C),                     (1)

with D_m depending only on the fixed planar field. In particular
h=lambda^2 psi^2/(W_C+sqrt(C)) and |psi|<=3Psi. Higher derivatives
follow by differentiating sqrt(C) sqrt(1+lambda^2 psi^2/C).
The constant advection U0 is removed by an exact coordinate change,
not placed inside a growing energy-estimate constant.

## 2. Full operator difference, including harmonic pressure

Use the planar variables of 0146/0161 in the moving frame:

    w_h=m+a+i k d b,  a=T_h X+y,  X_t=m,
    Z=b+i k v.X,     x=epsilon X,  epsilon=|k|.

Fix the planar length/time units, so epsilon<=1 below; physical units
are restored by epsilon=|k|/lambda. The response norm is

    N=|x|+|m|+||y||_Y+||Z||_2,  ||y||_Y=||curl_h y||_2.

The planar evolution has the uniform operator bound
||S_0(t)||<=M exp(c epsilon |t|) in this norm, by 0161's proved
horizontal group and 0146's exact augmented system. This holds for
arbitrary compatible initial augmented states in the vector sector,
not only one common-V trajectory.

The actual lifted-minus-planar velocity operator is

    Delta w_h,t=-i k h w_h-grad_h Delta pi,
    Delta b_t=-i k h b-w_h.grad h-i k Delta pi,
    Delta pi=2 i k(-Delta_h+k^2)^(-1)(grad h.w_h).         (2)

All operators in (2) use the complete periodic Green function. Its
harmonic pressure is exactly

    <Delta pi>=-2<h b>,                                 (3)

because div_h w_h=-i k b. Thus the inverse k^2 on that harmonic
does not create a hidden pole. In the excited representation (3)
vanishes by symmetry, but the estimate retains it. On other modes
the horizontal Green denominator is uniformly bounded below.

The full mean difference is

    Delta m_t=-i k<h w_h>,                              (4)

not zero. Curl eliminates the pressure gradient in the y equation,
but pressure in the vertical equation remains. From

    ||w_h||_{H^1}<=c N/epsilon,    ||b||_2<=c N,

equations (2)-(4) give

    ||Delta pi||_2<=c delta N,
    |Delta m_t|+||Delta y_t||_Y<=c delta N,
    ||Delta Z_t||_2<=c (delta/epsilon) N.

In particular the full difference is a BOUNDED operator on the
augmented response space with norm at most c delta/epsilon. The
apparently largest term is the actual axial shear
w_h.grad h: translation amplitudes may be O(1/k), and this power
has not been discarded.

Bounded-perturbation Duhamel, with no new spectral assumption, yields
for 0<=t<=T/epsilon

    ||S_C(t)||<=C_T exp(C_T delta/epsilon^2),
    ||S_C(t)-S_0(t)||
       <=C_T (delta/epsilon^2) exp(C_T delta/epsilon^2).  (5)

Both are bounds on actual full Euler evolution, not a truncated Bloch
matrix. They apply to the same phase and common-V data. Equation (5)
is obtained after controlling the rapid planar group; no exp(c/epsilon)
is silently beaten by a small coefficient.

For displacement data, choose the FULL Kelvin preparation of the new
background with the same material xi0=(X0,0)exp(iks). The difference
from the planar prepared velocity is

    P_k[-(X0.grad h)e_z].

Its contribution to N(0) is O(delta/epsilon)|epsilon X0|;
its horizontal pressure return is retained. This initial error is
bounded by the right side of (5) for epsilon<=1. Common-V data coincide
exactly. Each phase has its prescribed circulation data as in 0161.

## 3. Positive physical acoustic response and its exact quantifiers

For A0=|V0|+|epsilon X0|, (5) and the actual planar theorem imply

    sup_{0<=t<=T/epsilon}
      |m_C(t)-cos(epsilon c_b t)V0
          +epsilon c_b sin(epsilon c_b t)X0|
      <=C_T [epsilon+(delta/epsilon^2)
                         exp(C_T delta/epsilon^2)] A0,
    c_b^2=3Psi^2 lambda^2/4.                             (6)

The analogous integrated-position and compensated-current estimates
hold. The current is still the actual complete-fluid point mean plus
the explicit 0146 correction; (4) and the pressure/shear production
are not set to zero. Its slow derivative error is controlled by the
same right side: the additional current production is O(delta), hence
O(delta/epsilon) on slow time, which is no larger than (5).

The order of choices is concrete. Fix the planar field, T, a desired
relative error, and then a sufficiently small NONZERO epsilon. Choose
a finite C so that delta/epsilon^2 is sufficiently small using (1).
For example C proportional to epsilon^(-6), with a sufficiently large
fixed coefficient, gives delta/epsilon^2=O(epsilon). This produces
actual stationary generalized-force-free backgrounds and actual
histories with arbitrarily accurate positive acoustic response over
the full T/epsilon window. It is an ordered family, not a fixed-C
homogenization theorem or a claim that C has been measured or fitted.

## 4. Same material action, physical mass and positive phase energy

In the moving coordinates the material derivative is exactly
partial_t+v.grad_h+h partial_s. The Galilean transformation preserves
this derivative and the pressure Hessian; the known uniform axial
kinetic/background terms are not appended to the perturbation mass.

For xi0=(X0,0)exp(iks), phi=v.X0, the full fixed-Kelvin material rate
and convective rate are

    A xi0=-i k P_k[h(X0,0)+phi e_z],
    B xi0=i k h(X0,0).

The actual Jacobi Hamiltonian of this initial displacement phase is

    H_D=rho k^2/2 [||P_k(h(X0,0)+phi e_z)||^2
                         -<h^2>|X0|^2].                (7)

This follows by substituting pi=rho(A+B)xi into the COMPLETE material
Hamiltonian; the spatial mean pressure Hessian is zero. At h=0 it
reduces to the exact planar pressure-return stiffness in 0161/0167.
For k<=k_max and c0=lambda^2/(lambda^2+k_max^2), orthogonal projection
and the triangle/Cauchy inequalities give the useful explicit lower bound

    2H_D/(rho k^2)
      >=[c0 c_b^2-2delta c_b-delta^2]|X0|^2.             (8)

Thus sufficiently small delta gives a genuine positive phase stiffness
from the SAME Euler action. All added axial terms, including their
kinetic contribution, are present in (7).

The common-V phase has pi0=rho(V0,0), xi0=0, hence mass rho.
The full symplectic X0,V0 pairing is exactly rho X0.V0: the horizontal
mean of pi_D is zero because the mean-preserving Leray projection
cancels the h X0 mean between A and B. A possible initial cross energy
-rho Re[i k<h> V0*.X0] is retained; it is a known small axial-advection
term, not an additional mass. It vanishes on a common real standing
quadrature, and for the complete complex phase is O(delta) in normalized
(epsilon X0,V0) coordinates. Reducing delta further controls that cross
by the positive diagonal energy. One could instead use the specified
mean axial frame <W_C>; no such frame change is needed for (6)-(8).

On the actual solution family the full moving-frame action and its
connection are pulled through the observed (epsilon X,p_c) chart. By
(5)-(6) this chart and its slow derivative are close to the planar
invertible chart. Its exact symplectic/action pullback therefore retains
positive leading mass rho and stiffness rho k^2 C_v with the same
ordered error; it does not vary an on-shell profile as a freely chosen
static Cauchy--Born field. The physical mean and current map remain
explicit throughout.

Take nonlinear amplitude small LAST on each resulting finite interval.
No all-k fixed-background or uniform nonlinear-amplitude assertion is
needed to obtain these actual smooth Euler histories.

## 5. Optical interface and next stronger route

At the same core psi0=3Psi,

    W0=sqrt(C+9lambda^2 Psi^2),
    W'_0=3lambda^2 Psi/W0>0,
    W''_0=lambda^2 C/W0^3>0.

These are the actual inputs for the separate same-field local optical
construction. Global closed-torus/EPS embedding and an autonomous
coupled continuum are not inherited from this periodic array.

For ANY W=F(psi), the exact k=0 variable
r=b-(W'/zeta')eta is passively advected, since its coefficient is a
first integral. This suggests a stronger W-dependent slow corrector
that could improve (5) to fixed-C continuity. It is a distinct next
achievement; it is neither presumed by (6) nor necessary to the
established finite-C ordered family.
