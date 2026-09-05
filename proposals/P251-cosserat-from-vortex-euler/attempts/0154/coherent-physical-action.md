# The actual coherent Euler phase action in its physical mean chart

## 1. Full material phase space, including the pressure constraint

Use the actual stationary mean-zero field and actual physical initial
data of 0151. Pressure p is physical pressure divided by rho. At fixed
nonzero macro K, let H_K be the periodic Bloch-solenoidal displacement
space, P=P_K its full Leray projector, and A=u·∇_K. The real physical
encoding includes the conjugate −K sector. Cell inner products below
use unit-normalized real macro modes, so density rho is not changed by
a hidden cosine-average factor.

The exact Euler/Lin material action is

    L=rho/2[||η_t+Aη||²−〈η,Hess(p)η〉], η∈H_K.

Since η_t∈H_K, put C=P A on H_K and N=(I−P)A. The ACTUAL constrained
canonical momentum and Hamiltonian are

    π=rho(η_t+Cη),
    H(η,π)=||π||²/(2rho)−〈π,Cη〉
                +rho/2〈η,[Hess(p)−N* N]η〉.             (1)

The Hessian is understood as a quadratic form on H_K. C is skew-adjoint;
the ambient symplectic form is the standard

    Ω[(η1,π1),(η2,π2)]=〈η1,π2〉−〈π1,η2〉.

In particular the normal-pressure contribution −N* N remains in (1).
Replacing the full constraint by an unconstrained vector oscillator
would lose it. Completing the square in the action proves (1), without
assuming a constant-curl identity or a numerical pressure inversion.

## 2. Exact finite-K pairing of the actual initial columns

Let a fixed solenoidal macro mode D exp(iK·x) be the material displacement
column and let V exp(iK·x) be the independent common-velocity column.
For the first, the actual Kelvin preparation gives

    η_D(0)=D exp(iK·x),
    η_D,t(0)=−P[Aη_D+(Dη_D)^T u],
    π_D(0)=−i rho P[K(u·D) exp(iK·x)].                  (2)

For the second, η_V(0)=0 and π_V(0)=rho V exp(iK·x). The cancellation
of PAη_D in the third line of (2) is important. It is the constrained
material momentum, not the Eulerian velocity or its circulation label.

The field u has zero microscopic mean. Every term of π_D is therefore
on a nonzero microscopic sideband. For the fixed 0151 cell and |K|<q/8,
there is also no real-mode resonance with twice the macro wave. Hence,
on orthonormal physical transverse macro polarizations,

    Ω_DD=0,       Ω_DV=rho I,       Ω_VV=0              (3)

EXACTLY at finite K. The same is true for the small fixed-cell
perturbation in 0151B below its minimum microscopic wave number.
Each trajectory is a valid Euler trajectory with its own conserved
circulations; (3) does not put arbitrary common V on one fixed Kelvin
leaf. No new invariant or fitted inertia is supplied.

## 3. Average the actual phase action before choosing the mean coordinate

For realization R of the whole-field rotation/phase/time-reversal law,
let E_R(t) map the common initial coefficients z=(D,V) into the exact
material phase solution (η_R,π_R). Its columns are actual solutions:

    E_R,t=B_R E_R,     B_R=−Ω_R^(-1) H_R,
    E_R*Ω_R E_R=Ω0,    Ω0=rho [[0,I],[-I,0]].         (4)

The complete moving pullback has

    H_eff,R=E_R*H_R E_R+sym(E_R*Ω_R E_R,t)=0.          (5)

The connection cancels the actual Hamiltonian because the embedding
is propagated by that Hamiltonian. This does not set the physical
Hamiltonian to zero: z consists of conserved INITIAL coefficients.
Its pulled-back equation is z_t=0. Equation (5) is the general linear
Hamiltonian identity also implemented by `moving_phase_pullback`.

Now average the actual actions on the direct-integral family {E_R z}
with the SAME z, before charting the physical mean. Equations (3)-(5)
give exactly the phase form Ω0 and zero coefficient Hamiltonian.
Each realization retains its own Euler history; no common physical
history X_R=Xbar is imposed. The physical observation is the actual
full-fluid mean followed by the specified ensemble average,

    Xbar=F(t,K)D+G(t,K)V,
    mbar=F_t D+G_t V.                                 (6)

This is not an average of separately eliminated scalar inertias or
frequencies. The finite-window estimates from 0151 apply to these
particular rows, including their time derivatives.

## 4. Proper-rotation isotropy and real physical helicity charts

Proper SO(3) isotropy alone permits chirality. For a fixed K, F and G
on its transverse plane commute with rotations about K and have the
form a I+i b Jκ, where Jκ v=κ×v and a,b are REAL. Reality follows
from F(−K)=conjugate(F(K)) together with a proper rotation taking K to
−K; the same argument applies to G. Thus their spatial-helicity
eigenspaces are fixed and their scalar eigenvalues are real. A real
physical encoding is the cosine/sine circular macro pattern plus its
translation quadrature. These patterns have unit mean-square norm.

Use either actual helicity sector, with either of its real quadratures.
Equations (3)-(6) reduce there to

    Ω0=rho J0,  J0=[[0,1],[-1,0]],
    X=f(t,K) D+g(t,K) V,   m=f_t D+g_t V.              (7)

No mean displacement is renamed or its axis rotated to select a
positive clock. Both helicities have the same proven second spatial
jet; their potentially different chiral O(K³) terms remain in the
actual f,g. Time-reversal pairing makes f even and g odd in time,
but does not remove all spatial chirality. The complete transverse
representation is the direct sum of these physical sectors.

## 5. Exact physical Wronskian, action and momentum connection

Set

    T=[[f,g],[f_t,g_t]],   W=f g_t−g f_t.

Where W is nonzero, (X,m)=T z is an exact chart of the coherent phase
family. At t=0 it is the identity and W=1. Applying the actual moving
action, including (T^(-1))_t, gives

    M=rho/W,
    Kphys=rho(f_t g_tt−g_t f_tt)/W²,
    Bphys=[[0,1],[-Kphys/M,W_t/W]],
    Ωphys=M J0,
    Hphys=[[Kphys,M_t/2],[M_t/2,M]].                   (8)

These are the existing `physical_scalar_chart` formulas with the
actual physical displacement row (7), coefficient generator zero and
actual momentum row rho(f_t,g_t). The API's generic field names
`angle` and `spin` do not make this displacement a microrotation or
this mean momentum a tagged mechanical spin.

Up to the explicit endpoint −d(M X m/2)/dt, the first-order action is

    L=M m X_t−M m²/2−Kphys X²/2.

Variation in m gives m=X_t, so the EXACT scalar action and equation
on this prepared coherent family are

    L_X=M X_t²/2−Kphys X²/2,
    M X_tt+M_t X_t+Kphys X=0.                         (9)

Every moving connection is retained. This is now a pullback of the
actual averaged Euler phase action, rather than a Lagrangian supplied
after observing a force coefficient.

The physically measured full Fourier momentum density is

    Pphysical=rho m,
    Pcanonical=M m,
    Pphysical−Pcanonical=rho(1−1/W)m.                 (10)

The actual momentum row has no additional X connection in this chart.
Its bracket is {X,Pphysical}=W, not automatically one. The mismatch
in (10) remains explicit at exact finite K; it is not repaired by
renaming the canonical variable. Its size at the claimed order is
derived next.

## 6. Positive mass and stiffness at the actual controlled order

Write a=2v²/15. From the complete 0151 physical second spatial jet,
uniformly on each fixed interval [0,T] and with the required time
derivatives,

    f=1−a k²t²/2+O_T(k³),
    g=t−a k²t³/6+O_T(k³).

Substitution into (8), rather than imposing a normalized W, gives

    W=1+O_T(k³),           M=rho+O_T(k³),
    M_t=O_T(k³),           Kphys=rho a k²+O_T(k³).      (11)

The apparent second-order W terms cancel between the two ACTUAL
initial phases. Consequently physical and canonical momentum agree
through the complete second spatial jet, and the action (9) becomes

    L_X=rho X_t²/2−rho a k² X²/2+O_T(k³)              (12)

with the exact current correction (10) and all connections still
available. Both coefficients in (12) have been derived from the same
Euler action and physical observations.

More quantitatively, let δ=k/q and τ=qvT. The bounded full-operator
argument in 0151 gives uniform constants such that

    |W−1|+|M/rho−1| ≤ C_T |δ|³,
    |Kphys/(rho a k²)−1| ≤ C_T |δ|,                  (13)

after choosing |δ| small enough to keep W, for example, between 1/2
and 3/2. The constants have the same finite-window exponential bounds
as 0151, with the additional finite number of time derivatives supplied
by its bounded operator. Thus positive physical-chart mass and positive
acoustic stiffness hold for every selected sufficiently small nonzero
k on that window. No soft numerical sign test enters this conclusion.

For 0151B's actual fixed-cell elliptic tube, the C²-in-k finite-time
perturbation comparison gives instead M/rho=1+O_T(d δ²)+O_T(δ³),
M_t/(rho qv)=O_T(d δ²)+O_T(δ³), and
Kphys/(rho a k²)=1+O_T(d)+O_T(k/q), with dimensional powers of q,v
and the fixed dimensionless time understood as in 0151. The positive action survives the ordered small
d then small k/q choice, while these nonzero second-order connections
remain part of that perturbed action. Its exact initial pairing is
still (3). No claim that its complete second jet equals the unperturbed
one is made.

## 7. Why the order of averaging matters, and what remains active

As an exposing algebraic example only, consider two canonical scalar
solution frames with frequencies 0 and w, each of individual Wronskian
one and mass rho. The average physical row has

    f=(1+cos wt)/2,  g=(t+sin(wt)/w)/2,
    W=(2+2cos wt+wt sin wt)/4.

Its physical mass is rho/W, not the average of the two masses. W even
vanishes at t=pi/w although each individual chart is regular. This is
not an oscillator model inserted into Euler; it exposes the invalid
exchange of phase-action averaging and physical-coordinate elimination.

The present construction avoids that exchange by (4)-(8). It is an
exact time-dependent action on the actual coherent initial-data family,
with a positive controlled second-order local acoustic approximation.
It is not an autonomous invariant manifold for unrestricted Euler
initial states, and its estimates do not extend to acoustic time 1/k.
The scalar-chart obstruction at larger times and 0151's retained
orientation variance remain meaningful.

0153 can compare its same-field inserted histories and action on their
separately controlled finite-K window; a large periodic quadrature cell
does not inherit a uniform kP≪1 derivative estimate from (13). The next
parent achievement is the actual same-field optical/material-angle and
hybrid-current cross-action join with this acoustic family. No accepted
claim, EPS topology, spin normalization or parent closure is silently
supplied by the present acoustic action result.
