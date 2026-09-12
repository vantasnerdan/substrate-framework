# Fixed-tail compatibility and an exact charged Cao steady map

## 1. A periodic dressing cannot repair a nonstationary leading tail

Fix `0<gamma<1` and an integer `s>=4`.  On a common time interval let

    u(t,x)=U_{-2}(x)+R(t,x),
    U_{-2}(r n)=r^(-2) U(n),

where `U` is fixed in inertial coordinates and the remainder and its curl
obey the derivative-weighted bounds

    partial^a R=O(r^(-2-gamma-|a|)),       |a|<=s,
    partial^a curl R=O(r^(-3-gamma-|a|)), |a|<=s-1.       (A1)

These spatial bounds are the affine weighted hypotheses established for a
local interval in P253/0075 and independently reviewed in P253/0076.  Put
`Omega_{-3}=curl U_{-2}`.  Now additionally assume that the solution is
exactly `T`-periodic in fixed inertial coordinates, with (A1) uniform over the
period.  The exact vorticity equation is

    partial_t omega+[u,omega]=0.                              (A2)

Integrating (A2) over one period cancels the time derivative exactly:

    integral_0^T partial_t omega dt=omega(T)-omega(0)=0.     (A3)

The only degree `-6` term left in the integrated bracket is

    T [U_{-2},Omega_{-3}].                                   (A4)

Indeed, every bracket containing `R` or `curl R` has degree at most
`-6-gamma`.  Taking the degree-minus-six coefficient gives

    [U_{-2},Omega_{-3}]=0.                                   (A5)

For a relative-periodic state the same proof requires

    omega(T,x)-omega(0,x)=o(r^(-6))                          (A6)

with corresponding derivative control, or an endpoint symmetry that acts
trivially through this asymptotic order.  A spatial translation or rotation
generally changes subleading coefficients and is not automatically harmless.
A compact internal phase or core deformation is covered only when it obeys
(A6).  A rotation of `U` is a different asymptotic sector.
Consequently the nonzero degree-six residual of the fixed-frame texture in
P253/0075 cannot be repaired by a bounded faster periodic dressing.  This is
a compatibility theorem for one affine asymptotic class, not a classification
of all homogeneous degree-minus-two Euler fields.

## 2. Traveling Maxwell reduction

Use physical mass density `rho_m`, permittivity `epsilon_EM`, permeability
`mu_EM`, coupling `g`, and

    c_EM^2=1/(epsilon_EM mu_EM).

Let `y=z-c t` and write all fields as functions of `(r,y)`.  With

    E=-grad phi-partial_t A,       B=curl A,

Faraday's law gives the exact comoving gradient

    E_c:=E+c e_z cross B=-grad Phi,
    Phi:=phi-c A_z.                                         (5)

In Lorenz gauge

    div A-(c/c_EM^2) partial_z phi=0,

the potentials solve

    L_c phi=(g/epsilon_EM) chi,
    L_c A=mu_EM g chi u,
    L_c=-Delta+(c^2/c_EM^2)partial_zz.                       (6)

The operator is uniformly elliptic precisely when

    a_c:=1-c^2/c_EM^2>0.                                    (7)

Its fundamental solution is

    G_c(r,z)=1/[4 pi sqrt(a_c r^2+z^2)].                     (8)

Thus `phi-(Q_g/epsilon_EM)G_c` lies in a decaying weighted space, whereas
`phi` itself need not be in unweighted `L2`.  The Gauss row

    Q_g=g integral chi dx                                  (9)

is an explicit scalar constraint.  The finite-energy fields `E` and `B`
have degree-minus-two tails.  The analogous leading coefficient of `Phi`
also contains the total axial current through `A_z`; it is fixed by (6), not
silently set to a zero-charge potential class.

## 3. Exact axisymmetric no-swirl closure

Let the relative fluid velocity be

    W=u-c e_z=(1/r) grad P cross e_theta
      =-(P_z/r)e_r+(P_r/r)e_z.                              (10)

The laboratory streamfunction is `psi=P+c r^2/2`.  Because
`Delta_* r^2=0`, where

    Delta_* P=P_rr-r^(-1)P_r+P_zz,

the potential vorticity is

    zeta=omega_theta/r=-r^(-2) Delta_* P.                   (11)

Choose a material tag `chi=chi(P)` on a regular tagged cell.  The stationary
continuity equation is `div(chi u)-c partial_z chi=div(chi W)=0`.
An axisymmetric poloidal current has a poloidal vector potential and a
toroidal magnetic field.  Write

    B=B_theta e_theta=r H e_theta.                          (12)

Equations (5), (10), and (12) give

    E+u cross B=-grad Phi+W cross B
               =-grad Phi-H grad P.                        (13)

The theta curl is not generally zero.  Direct cylindrical differentiation
gives the exact identity

    r^(-1) curl[g chi(-grad Phi-H grad P)]_theta
      =g W dot grad[chi'(P)Phi-chi(P)H].                    (14)

The forced steady vorticity equation therefore integrates on every connected
regular streamline cell to

    zeta=h(P)+(g/rho_m)[chi'(P)Phi-chi(P)H].                (15)

This is the modified Grad--Shafranov row.  It retains the electric curl
`-g grad chi cross grad Phi`; absorbing `chi grad Phi` into pressure would
lose this term.

The radial component of Ampere's law integrates in `z`.  Axis regularity and
finite field energy exclude the source-free toroidal field `B_theta=C/r`, so

    (1/mu_EM-epsilon_EM c^2) r H+epsilon_EM c Phi_r
       =g K(P)/r,             K'(P)=chi(P).                 (16)

The `z` component of Ampere's law together with Gauss's law then gives

    -Phi_rr-r^(-1)Phi_r-a_c Phi_zz
      =(g chi/epsilon_EM)
        [a_c-c P_r/(c_EM^2 r)].                            (17)

Conversely, (16)--(17), stationary continuity, and the affine weighted decay
of `Phi` and `H` imply the full axisymmetric Maxwell equations after
reconstructing

    E=-grad Phi+c r H e_r,       B=r H e_theta.             (17a)

The decay kills the `z`-integration function in the radial Ampere primitive;
axis regularity alone does not remove every such row.  The Lorentz force is
poloidal, so no swirl is generated.

Finally,

    W dot grad W=grad(|W|^2/2)+zeta grad P.

Using (15) in the momentum equation gives the Bernoulli first integral

    p+rho_m |W|^2/2+g chi(P) Phi=b(P),
    b'(P)=-rho_m h(P).                                     (18)

Equations (15)--(18) are an exact reciprocal steady fluid--field reduction,
not a supplied-force model.

## 4. The exact charged Cao map

Fix one sufficiently thin exact Cao ring with exponent `p>=6`, circulation
`kappa`, speed `c_0`, and even axial center.  Its exact uncharged row is

    -r^(-2)Delta_* P_0=epsilon_core^(-2)(P_0-mu_0)_+^p,
    psi_0=P_0+c_0 r^2/2 ->0 at infinity.                   (19)

Here `epsilon_core` is Cao's core parameter and is unrelated to
`epsilon_EM`.  Cao's auxiliary parameter system (3.36) is not an equation of
the map below.

On a compact regular band of the positive core define `I_P` to be enclosed
physical volume action and choose a fixed smooth profile `bar_chi` with

    chi_P=bar_chi(I_P),       integral chi_P dx=1.          (20)

The action label, rather than the numerical value of `P`, fixes the tag
distribution as the carrier deforms.  The support of `bar_chi` stays away
from the core center, the free boundary, and every point where
`partial_I zeta_0=0`.

For a prospective charged member set `Q_g=g` and solve for
`(P,Phi,H,mu,c)`:

    F_1:=-r^(-2)Delta_*P
          -epsilon_core^(-2)(P-mu)_+^p
          -(g/rho_m)[chi_P'(P)Phi-chi_P(P)H]=0,            (21)

    F_2:=-Phi_rr-r^(-1)Phi_r-a_c Phi_zz
          -(g chi_P/epsilon_EM)
             [a_c-c P_r/(c_EM^2 r)]=0,                    (22)

    F_3:=(a_c/mu_EM)rH+epsilon_EM c Phi_r-gK_P(P)/r=0,     (23)

together with exact circulation, axial impulse, axial center, gauge, axis,
and affine whole-space decay rows.  The laboratory velocity and pressure are
reconstructed from (10), `u=W+c e_z`, and (18).  Since `p>=6` and the tag is
supported strictly inside the positive core, no electromagnetic source or
new vortex sheet occurs at the Cao free boundary.

At `g=0`, `Phi=H=0` and (21) is exactly (19).  The Maxwell derivative in
`(Phi,H)` is the elliptic triangular block (22)--(23).  Its inverse maps the
compact current to `Phi,H=O(g)` in the affine weighted potential space and
to finite-energy `E,B=O(g)`.  Substitution in (21) is `O(g^2)`.

Equivalently, and more cleanly for the Banach map, solve (6) for
`phi=g hat_phi` and `A=g hat_A`.  Then

    L_c(hat_phi)=chi_P/epsilon_EM,
    L_c(hat_A)=mu_EM chi_P(W+c e_z),
    hat_Phi=hat_phi-c hat_A_z.                             (23a)

The reduced fluid map depends smoothly on

    tau=g^2>=0,                                           (23b)

while `phi,A,E,B` are odd in the signed charge `g`.  Applying an implicit
function theorem in `tau` proves the fluid, speed, and tag corrections are
`O(g^2)` and prevents a spurious order-`g` fluid row.  After solving in
`tau`, either sign of `g` reconstructs the corresponding signed charge.

## 5. The tag stabilizer is fixed on a nonisochronous band

Use regular volume action coordinates `(I,beta,theta)` with physical volume
`dI d beta d theta`.  Here the relative velocity is tangent to the meridional
streamline,

    W_0=omega(I) partial_beta,

whereas the physical toroidal vorticity is

    omega_0=omega_theta e_theta=zeta_0(I) partial_theta.   (24)

The last equality uses `partial_theta=r e_theta`, so
`zeta_0=omega_theta/r`; it is not a meridional vorticity component.

For a smooth compact stabilizer displacement
`xi=xi^I partial_I+xi^beta partial_beta+xi^theta partial_theta`, the first two
components and the toroidal component of `[xi,omega_0]=0` are

    -zeta_0 partial_theta xi^I=0,
    -zeta_0 partial_theta xi^beta=0,
    xi^I zeta_0'-zeta_0 partial_theta xi^theta=0.           (25)

Period integration in `theta` gives `xi^I zeta_0'=0`.  On the support chosen
in (20), both `zeta_0` and `|zeta_0'|` are bounded below, hence

    xi^I=0,       xi dot grad chi_0=bar_chi'(I)xi^I=0.      (26)

The `beta` and `theta` components are tangent to the action torus and do not
change the tag.  The claim is restricted to the declared nonisochronous
regular band;
it is false without that support restriction.  Exceptional levels, the
center, boundary, and exterior carry no tag and are handled by smooth zero
extension.  Thus fluid-vorticity stabilizers cannot rearrange the chosen
charge while leaving the tagged carrier unchanged.

## 6. Exact linear range and the remaining Schur row

Work in the even-in-`z` axisymmetric slice.  The axial translation mode is odd
and is absent.  Linearizing (21) at `g=0` in the laboratory streamfunction
uses

    varphi=delta psi,       delta P=varphi-r^2 delta_c/2.   (26a)

Thus the exact Cao operator is

    L_epsilon varphi=-r^(-2)Delta_*varphi
       -p epsilon_core^(-2)(P_0-mu_0)_+^(p-1)
          (varphi-delta_mu-r^2 delta_c/2),                 (27)

augmented by the linearized circulation and impulse rows.  Its blown-up
principal operator is

    -Delta-p U_+^(p-1),                                   (28)

whose bounded kernel is exactly the two translations by Cao Lemma 3.8.  The
even slice removes the axial translation.  The radial translation has
nonzero impulse pairing

    delta I_z=2 pi rho_m r_* kappa delta r+o(delta r),     (29)

and speed variation has a nonzero radial-translation projection.  The
circulation row removes the constant/profile column.  A contradiction and
blow-up argument of the type used in Cao Lemma 3.9 is the natural route to an
inverse on the complement of the radial translation cell, but the source does
not automatically provide that inverse for the present map.

The intended complement operator is

    L_epsilon^perp:X_epsilon^perp -> Y_epsilon^perp        (30)

on an even weighted axis/decay streamfunction space.  This attempt has not yet
constructed the exact `X_epsilon,Y_epsilon`: their norms must include the
axis weight, affine `P` versus decaying `psi` row, core-interface regularity,
whole-space Green mapping, and differentiability of
`P -> I_P -> chi_P`.  Cao Lemma 3.9 supplies a concentrated a-priori mechanism
under its stated support and norms, but it does not by itself give range
closure or surjectivity for (30).  Denote the missing theorem by

    HSE: L_epsilon^perp is a Fredholm index-zero isomorphism
         on the fully declared exact steady-map spaces.                  (30a)

The limiting Lane--Emden kernel and compact core coefficient identify a
concrete proof route to HSE; source naming does not establish it.

Conditional on HSE, let `B_epsilon=(partial_mu F_1,partial_c F_1)` and let
`R_epsilon=(D kappa,D I_z)` be the two physical constraint rows.  Eliminating
`X_epsilon^perp` gives the exact finite-dimensional matrix

    mathscr_S_epsilon
      =D_(mu,c)(kappa,I_z)
         -R_epsilon (L_epsilon^perp)^(-1) B_epsilon,
    S_epsilon=det mathscr_S_epsilon.                       (31)

The impulse row normalizes the radial generalized kernel, and the circulation
row fixes the profile/chemical-potential column.  This is the exact
speed--impulse Schur determinant; it contains no arbitrary adjoint
normalization.  If HSE holds and `S_epsilon!=0`, the augmented derivative is
an isomorphism.
The implicit-function theorem in `tau=g^2` then
produces an exact branch satisfying

    P_g-P_0=O(g^2), c_g-c_0=O(g^2),
    Phi_g,H_g=O(g),                                       (32)

for all sufficiently small `|g|`.  The centered Maxwell forcing is even in
`z`, while the axial translation cokernel is odd, so its exact pairing is
zero; this is the concrete translation solvability row rather than an appeal
to integrated self-force.

Route B1 therefore has two exact remaining achievements: prove HSE, then
evaluate (31), including its physical impulse normalization and sign.  The
source asymptotics strongly
suggest it is nonzero, but importing (3.36) as (31) would repeat the error
identified in P253/0057.  The charged branch is therefore conditional on
HSE and `S_epsilon!=0` at this attempt boundary.  These are route-specific
missing constructions, not evidence against the carrier or the gauge
extension.

## 7. Continuation ladder

Route B1 has produced the exact two-scalar steady map, the whole-space affine
Maxwell inverse, the modified Grad--Shafranov compatibility, tag-stabilizer
locking, and the axial cokernel cancellation.  Its next analytic achievement
is to construct HSE from the exact Green equation and then evaluate (31) from
the physical circulation/impulse rows.

If that determinant vanishes, Route B3 activates with its actual radial
kernel and a finite-dimensional bifurcation equation.  Independently, Route
B2 can use the exact reduced Maxwell functional on the fixed tag/circulation/
impulse class; it must compute the signed joint Hessian rather than inherit
the uncharged maximizer's stability.  Neither alternative uses the dynamical
zero-frequency Kelvin resolvent, whose `m=0` eigenvalues accumulate at zero.

Even after (32), all-time charge-profile persistence requires the joint
energy--Casimir/electric `H^-1` coercivity estimate on the same weighted orbit
topology.  Charge and action quantization, an analyzer, Born statistics, a
Lorentz cone for the incompressible pressure sector, and any electron or
neutrino identification remain separate active obligations.
