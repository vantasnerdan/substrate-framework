# The inherited stationary Euler action with its actual internal cell

Central preregistration/schema266 was confirmed before this derivation.
Use0170/C-CST-013 conventions: stationary periodic mean-zero Euler u,
pressure p per unit density, A=u.grad, a=kappa.u, K=k kappa, with one
common laboratory kappa, and the full microscopic mean-preserving Leray P.
All cell brackets are full normalized volume averages. This construction
applies before whole-field O(3)/TR averaging, which is performed on the
same common phase inputs and actual action.

## 1. Expand the actual constrained material action

The exact Euler Jacobi action on div_K eta=0 is

    L_K=(rho/2)<|eta_t+(A+ik a)eta|²-eta* Hess(p) eta>.   (1)

Its variation is0170's pressure-constrained Jacobi equation. The pressure
Hessian is not discarded or replaced by a fitted restoring force.
Use real parity coefficients in the complex Bloch expansion

    eta=Y+ik chi+k² zeta+O(k³),
    <chi>=<zeta>=0, div chi=0, div zeta=kappa.chi,
    kappa.Y=0.                                         (2)

Y is the MEAN MATERIAL displacement, not the physical velocity-integral
mean X. The zeroth action is rho|Y_t|²/2. In the second coefficient,
<A zeta>=<zeta_t>=0 and

    <Y.Hess(p) zeta>=-<(Y.grad p)(kappa.chi)>.

Thus the unknown second cell cancels from the ACTION itself, not merely
from its mean equation. The inherited stationary second-jet action is

    L=rho|Y_t|²/2+(rho k²/2){
        <|chi_t+A chi+aY|²-chi.Hess(p)chi>
        -2< a Y_t.chi>+2<(Y.grad p)(kappa.chi)>}.        (3)

No inverse microscopic wave number appears in(3). The solenoidal chi
constraint still carries its full pressure multiplier. The neglected
terms are higher spatial jets for a chosen finite cell; this is not an
arbitrary-k replacement of Euler by a few oscillators.

## 2. Both variations reproduce the actual equations

Vary chi in its zero-mean solenoidal space. Since A is skew adjoint and
A a=-kappa.grad p, (3) gives exactly

    chi_tt+2PA chi_t+P(A²+Hess p)chi
       =-2P(aY_t)+P[(kappa.grad p)Y+kappa(Y.grad p)].    (4)

The two pressure terms have different origins: transport of a, and the
solenoidal second-cell constraint. Omitting either changes the actual
force. Put C=-A+Du and L_E=-P(A+Du). Stationary Euler gives
A(Du)+(Du)²=-Hess p. With w=chi_t-Cchi, (4) is

    chi_t=Cchi+w,
    w_t=L_E w+F_Y Y+F_V Y_t,
    F_Y Y=P[(kappa.grad p)Y+kappa(Y.grad p)],
    F_V V=-2P(aV).                                     (5)

These are the actual first-cell Euler velocity and Lin displacement,
not an auxiliary oscillator interpretation.

Vary Y in the physical transverse plane. Its conjugate momentum is

    Pi_Y=rho Y_t-rho k² P_kappa<a chi>.                 (6)

Using <a Achi>=<(kappa.grad p)chi> gives

    Y_tt=k² P_kappa{<a²>Y+2<a chi_t>
                      +<(kappa.grad p)chi+grad p(kappa.chi)>}. (7)

This is0170's actual MATERIAL mean equation. In particular the factor2
belongs here, and not in the physical Euler mean.

## 3. The physical mean and initial phase are kept distinct

Lin's actual velocity observation, not variation of a guessed effective
mass, gives

    m=Y_t-k² J chi+O(k³),
    J chi=P_kappa<a chi-u(kappa.chi)>,
    X=D+integral_0^t m=Y-k² integral_0^t Jchi+O(k³).     (8)

Consequently

    rho m=Pi_Y+rho k² P_kappa<u(kappa.chi)>+O(k³).       (9)

This is an actual canonical/measured momentum difference and its current.
The exact inherited second jet can be written locally with an integrated
current variable I_t=Jchi and X=Y-k² I, I(0)=0. This variable is a
physical accumulated current, not a new adjustable modulus.

The bare Kelvin-D/common-V preparation is

    Y(0)=D, Y_t(0)=V, chi(0)=0,
    chi_t(0)=w(0)=-P[aD+kappa(u.D)].                    (10)

Its internal cotangent is
Pi_chi=rho k² P(chi_t+Achi+aY), hence
Pi_chi(0)=-rho k² P[kappa(u.D)]. The initial internal displacement
vanishes for both D,V variations, and Pi_Y(0)=rho V. The restricted
common initial symplectic form is therefore exactly rho J through this
jet, in agreement with the full material phase result0154/0170. Replacing
(10) by zero internal rate changes a real initial circulation/preparation
row; it is not a harmless choice of origin for the same response.

## 4. Scope of positivity and time

Equation(3) has the kinetic and gyroscopic structure inherited from Euler,
but its full internal potential need not be positive: Hess p and transport
are actual operators. No positive finite optical mass is supplied by
calling chi an optical mode. The positive C-CST-011 packet phase remains
its own actual source; whether its registered angle/spin closes the stress
functional is the next test.

On any fixed finite window, (7),(8) give Y-X=O(k²). Replacing Y,Y_t by
X,X_t in the forcing(5) changes its solution by O_T(k²), hence changes
the observed k² stress only atO_T(k⁴), within a second-jet statement.
This uses the actual bounded finite-time Euler/Lin group and the fixed
cell, not a spectral gap, uniform Bloch radius or acoustic-time estimate.
The variational parent remains(3) with Y and the measured current(8);
substituting X in the response does not license dropping(9) from its action.

## 5. Corrected preparations carry a computable initial phase correction

For actual corrected cells chi0=xi(D,V), chi_t0=v(D,V), fixing the PHYSICAL
initial mean velocity m0=V requires Y_t0=V+k² Jxi. Therefore

    Pi_Y0=rho V-rho k² P_kappa<u(kappa.xi)>,
    Pi_chi0=rho k² P[v+Axi+aD].

The restricted initial one-form is consequently

    Theta0=rho V.dD+rho k²{
       -P_kappa<u(kappa.xi)>.dD
       +<v+Axi+aD, dxi>}.                              (11)

Its symplectic form is the same -dTheta0 convention used for the inherited
phase action. Formula(11) gives the exact finite correction to the common
rho J row at this order; it need not vanish. Initial optical phase/cross
conditions, when present, must use this actual form and the same physical
angle/spin rows. Merely preserving the mean displacement D while changing
the first cell does not preserve either physical initial velocity or its
phase mass automatically. The bare preparation xi=0 recovers the earlier
rho J identity exactly.
