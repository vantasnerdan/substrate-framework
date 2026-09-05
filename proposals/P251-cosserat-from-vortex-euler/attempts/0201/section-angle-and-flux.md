# A measured section angle and its literal angular-momentum flux

This companion uses the positive m=0 mode of `radial-mode-and-action.md`.
All quantities are first variations of the same Euler fluid. A section
observable is not a material parcel: the distinction is essential here.
Pressure is divided by density in the mode equations, but every angular
moment below includes the actual density rho.

## 1. The complete-fluid section current

Write V=r Omega, L0=r V and let delta L=r v_theta. Axisymmetry eliminates
the azimuthal pressure torque, not the radial and axial advective flux.
The EXACT linear conservative equation is

    (delta L)_t + (r L0 v_r)_r/r
       + (W delta L+L0 v_z)_z = 0.                 (1)

For a smooth radial window chi, define

    S_chi = 2pi rho integral chi r² v_theta dr,
    F_chi = 2pi rho integral chi r²(W v_theta+V v_z) dr.

Then

    (S_chi)_t+(F_chi)_z
       =2pi rho integral chi' r² V v_r dr.          (2)

The right side is a physical radial flux. For the complete section,
chi=1 and the actual Bessel tail removes its boundary value at infinity.
The perturbation angular momentum and flux converge even though the
unperturbed column's unsubtracted angular momentum does not. Thus this
is a RELATIVE complete-section moment, not a finite background moment.

For the real mode with phase S=kz-omega t and positive radial y,

    S_E=-S0 cos(S),
    S0=2pi rho integral r Z y dr >0,
    F_E=c S_E.                                    (3)

The last equality follows from (1), not from discarding exterior flow.
In particular F_E includes the nonzero exterior term V v_z even where
Z=W=v_theta=0. A finite window outside the vorticity core has the SAME
S_E but a nonzero radial return in (2). It cannot be used as a closed
material spin balance.

## 2. An actual nonnegative painted marker measures the angle

Choose chi supported in the open flat-Omega inner core and 0<epsilon<1.
The nonnegative passive marker

    a0(t,r,theta,z)=chi(r)[1+epsilon cos(2(theta-Omega_c t))]

is an actual advected material fraction, independent of z. The shear W
does not invalidate its base transport. Its perturbation is the actual
Lin pullback delta a=-xi.grad(a0), not a freely prescribed readout.
At each fixed Eulerian z its centroid is zero by the m=0 and m=2 angular
selection rules. The measured central Euclidean quadrupole is

    Q_E=integral r² exp(2i theta) a r dr dtheta.

Let D_chi=integral chi r³ dr>0. Its background value is
pi epsilon D_chi exp(2i Omega_c t). Angular integration gives

    exp(-2i Omega_c t) delta Q_E
       =-pi epsilon integral r³ chi' xi_r dr
          +2i pi epsilon integral chi r² xi_theta dr.             (4)

Consequently the actual incremental covariance angle, relative to the
transported background quadrupole, is

    theta_E = integral chi r² xi_theta dr / D_chi.                 (5)

This has unit response to a rigid rotation xi_theta=r theta. The radial
displacement in (4) changes the quadrupole magnitude and is retained;
it does not change the first-order angle. The observable is evaluated
at a fixed section; it is not the covariance of a fixed set of particles
followed through the axial flow.

Using xi_theta=i B exp(iS), B=-2Omega_c f/[k(c-W)], yields

    theta_E=C_theta sin(S),
    C_theta= integral chi r² 2Omega_c f/[k(c-W)] dr / D_chi >0.

Combining with (3) proves the all-time linear identity

    S_E=J_E (theta_E)_t,
    J_E=S0/(omega C_theta)>0.                     (6)

The spin in (6) is the COMPLETE-fluid relative section moment, while
the marker supplies an independently defined angle. It is not the
angular momentum of only the painted fraction. Weighting spin by that
fraction introduces its own delta a and boundary/selection terms and
does not reproduce (6).

The inherited scalar angle-action mass is

    M_E=-beta/(omega C_theta²)>0.

The actual overlap eta_E=J_E/M_E is positive and computed, not set to
one. The fixed axial phase measure in beta must also be included in
S0 when both are interpreted per the same phase volume. Equivalently
use beta/L_z and section S0 in this formula. No normalization by a
carrier-dependent period is made when differentiating carrier jets.

## 3. Material spin conservation and the full current improvement

The same mode has

    delta_material(r u_theta)=r v_theta+(r V)'xi_r=0               (7)

pointwise. This is fully compatible with (6). Eulerian angular momentum
changes because particles carrying different background angular momentum
cross the section. Since div xi=0, equation (7) implies

    delta L=-div(L0 xi),
    S_E=-partial_z R_E,
    R_E=2pi rho integral r² V xi_z dr.              (8)

The radial surface term vanishes only at infinity. Direct insertion of
the mode in (8) gives integral r V y'=-integral r Z y, verifying (3).

For the literal displacement moment and its axial transport, put

    G_E=2pi rho integral r² xi_theta dr,
    F_G=2pi rho integral r² W xi_theta dr,
    C_shape=2pi rho integral r³ Omega' xi_r dr.

The actual cylindrical Lin relation gives

    S_E=(G_E)_t+(F_G)_z-C_shape.                    (9)

Both F_G and C_shape are retained. In particular (6) does NOT license
G_E=J_E theta_E. Integrating (9) supplies the physical current with its
initial value, shape exchange and transport. Equation (1) has no
azimuthal pressure torque because m=0; its nonlocal pressure is already
present in the actual v_r,v_z and hence every flux above.

The established object is therefore a positive-action actual Euler
mode, a unit-response physical covariance angle, and nonzero positive
complete-fluid SECTION spin/rate overlap with exact flux. Literal
material torsional spin is zero. A coupled local constitutive use must
carry (8)-(9), not rename the section current as conserved parcel spin.
