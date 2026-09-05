# An actual low-frequency Euler mode, its centroid, and its action

This calculation establishes the rotating-fluid m=1 inertial/Bessel
sector and its actual translational observation. It also determines why
its laboratory zero-frequency crossing, by itself, is **not** a
positive finite-mass acoustic sector. The physical laboratory clock and
carrier are fixed throughout. No sign is selected by changing a Floquet
logarithm or subtracting a rotation after the calculation.

## 1. Exact Euler and material reconstruction

Use the stationary reference u0=Omega e_z cross x, Omega>0. Follow any
constant axial drift by an actual translation of coordinates. In the
frame rotating at Omega, the perturbation velocity w and pressure obey

    w_t+2 Omega J w=-grad pi,  div w=0,                 (1)

where J is e_z cross on the transverse plane and zero axially. This
follows from the inertial Euler equation, including the vector-frame
rotation and convected spatial argument. Write

    pi=P J_m(kappa r) exp(i m theta+i k z-i sigma t),
    D=4 Omega^2-sigma^2.

For sigma!=0 and D!=0, the exact velocity is

    w_r=i P[ sigma*kappa*J_m'(kappa r)
                       -2 Omega*m*J_m(kappa r)/r ]/D,
    w_theta=P[2 Omega*kappa*J_m'(kappa r)
                       -sigma*m*J_m(kappa r)/r]/D,
    w_z=k*pi/sigma,                                   (2)

with the common exponential understood in the first two lines. Direct
divergence gives the Bessel equation and

    sigma^2=4 Omega^2 k^2/(kappa^2+k^2).                (3)

No pressure component has been omitted. The material displacement in
the rotating frame is zeta=i*w/sigma. Indeed the Lin velocity is
w=zeta_t. Applying the full projector to (1) gives
2 Omega P J w=i sigma w and hence

    P[zeta cross (2 Omega e_z)]=w.

Thus this is the fixed-Kelvin material solution, not just a velocity
ansatz. In the laboratory, xi=R_(Omega t) zeta(R_(-Omega t)x,t), and the
cylindrical components have frequency

    omega_lab=sigma+m Omega.                           (4)

Choose m=1 and the negative branch sigma=-2 Omega k/sqrt(kappa^2+k^2)
for k>0. At k0=kappa/sqrt(3), sigma=-Omega and omega_lab=0. The Eulerian
velocity AND material displacement are stationary laboratory patterns
at this crossing. Their carrier remains exp(i theta+i k0 z); it has
not become a spatially uniform displacement.

At fixed free-space kappa, the longitudinal sideband k=k0+K has

    omega_lab= -3 sqrt(3) Omega K/(4 kappa)
                 +27 Omega K^2/(32 kappa^2)+O(K^3),
    omega_lab^2=27 Omega^2 K^2/(16 kappa^2)+O(K^3).      (5)

The positive coefficient in the last line is a true carrier-sideband
coefficient. It is not yet a sound speed for a mean displacement field.

## 2. Normalizable action and the boundary-dependent sideband

The entire Bessel beam is not L2 transversely. To evaluate a finite
action, first state a concrete reference domain: a circular cylinder of
radius b_c, with the impermeable/slip condition w_r=0, and an axial
period (or a normalized Bloch period for nearby k). The background
rotation is tangent to this boundary. This is an exact Euler reference
problem, not an asserted boundary condition on the EPS torus.

Writing x=kappa b_c, y=k b_c and s=sigma/Omega, the boundary condition is

    F(x,s)=s*x J1'(x)-2 J1(x)=0.                        (6)

At the zero-frequency crossing s=-1 this becomes
x J1'(x)+2 J1(x)=x J0(x)+J1(x)=0. There are positive roots, by the
oscillating large-argument Bessel asymptotics. Each is simple: at a
root the derivative of x J0+J1 is -(x+3/x)J1, nonzero because J0 and
J1 cannot vanish together (uniqueness for their Bessel recurrence/ODE).

The actual fixed-cylinder branch must vary x when y varies; holding
kappa fixed would violate (6). At y=x/sqrt(3), differentiation of
F(x,s(x,y))=0, with s=-2y/sqrt(x^2+y^2), gives

    dx/dy=-3 sqrt(3)/(2 x^2+3),
    d omega_lab/dk =-Omega*b_c *
           3 sqrt(3)(x^2+3)/[2 x(2 x^2+3)].             (7)

This nonzero derivative yields a genuine low-frequency branch of this
fixed boundary problem. It differs from the free-Bessel slope in (5).
No numerical Bessel root or fitted frequency is needed to establish it.

For a fixed real frequency on this branch, pressure integration by parts
is legitimate: radial flux vanishes and opposite axial/Bloch faces cancel.
Taking the complex inner product of (1) with w gives

    2 Omega integral w_bar·Jw = i sigma ||w||_2^2.

The exact material KKS form is
Omega_KKS(xi1,xi2)=rho integral (2 Omega e_z)·(xi1 cross xi2).
For the complex mode xi_hat=i*w_hat/sigma this yields

    Omega_KKS(xi_hat_bar,xi_hat)=-i rho ||w_hat||_2^2/sigma.

If E1=Re xi_hat, E2=Im xi_hat, define

    beta=Omega_KKS(E1,E2)=-rho ||w_hat||_2^2/(2 sigma).   (8)

It is strictly positive and finite at sigma=-Omega, for a fixed
nonzero mode normalization. This is the actual complete reference-domain
symplectic form; a small observed parcel is not substituted into it.

With the convention L=-z^T Omega_phase zdot/2-z^T H z/2, put
Omega_phase=-beta J. A real laboratory mode has coefficient history
z=(cos(omega_lab t),sin(omega_lab t)), so A_lab=omega_lab J. Therefore

    H_lab=-Omega_phase A_lab=-beta omega_lab Id.         (9)

The corresponding rotating-frame coefficient is -beta sigma>0. The
difference is the physical angular-momentum/frame term m Omega, already
fixed in (4); discarding it would replace the laboratory energy by a
different generator. In the laboratory (9) changes sign LINEARLY across
the crossing. Squaring the frequency in (5) or (7) cannot make this
Hamiltonian positive on both sides.

Scalar elimination, where omega_lab!=0 and the coefficient is positive,
gives I_can=-beta/omega_lab. For a physical centroid normalization it
is divided by the squared, finite centroid amplitude computed below.
Thus it diverges as 1/K at the crossing. This is not a finite positive
mass tending to the parcel's ordinary translational mass.

## 3. The actual finite material centroid is visible, but its momentum vanishes

Choose an observed material cylinder D with radius b<b_c and axial
half-length d, centered on the rotation axis. It is invariant under the
unperturbed rotation; its boundary is transported in the perturbed
solution. Its mass is M_D=2*pi*rho*b^2*d. For m=1, exact angular
integration of (2) gives

    average_D w_perp =
      -i P J1(kappa b)/[(2 Omega+sigma)b]
                   sinc(k d) (1,i) exp(-i omega_lab t),
    delta X_D =
       P J1(kappa b)/[sigma(2 Omega+sigma)b]
                   sinc(k d) (1,i) exp(-i omega_lab t). (10)

The laboratory vector-frame factor is included in this expression.
The identities used are J1'+J1/x=J0 and
integral_0^b r*kappa*J0(kappa r)dr=b*J1(kappa b).
Choose b,d away from these explicit zeros. Then the physical centroid
amplitude is nonzero and remains finite at sigma=-Omega.

The ACTUAL material mean momentum variation is

    delta P_D=rho integral_D [v+(xi·grad)u0]
             =M_D d_t delta X_D
             =-i omega_lab M_D delta X_D.              (11)

The equality uses the Lin equation and the invariant reference tag,
not the Eulerian integral of v alone. In particular (10)'s velocity
average is generally nonzero at omega_lab=0, while (11)'s physical
momentum is zero. The displacement-induced change of the background
rotation is the indispensable other term. A velocity-only average
would falsely identify the stationary internal deformation as a moving
mass.

For the COMPLETE axial carrier period, sinc(k d)=0. Its total transverse
centroid and mean momentum vanish even away from the frequency crossing.
A shorter tag can see the internal displacement, while the rest of the
carrier period supplies the compensating translation. This is an actual
observation, not a no-go for tagged physical fields; it specifies which
field is and is not supplied by this mode.

Normalize its nonzero tagged centroid amplitude to one. Then at K=0 the
rows (X_D,P_D) have determinant zero, while (8)'s symplectic form is
nondegenerate. At small K their physical determinant is O(M_D omega_lab),
whereas beta has a nonzero limit. This is the explicit mechanical/action
mismatch. It cannot be repaired by naming omega_lab^2's coefficient a
sound speed. A construction with an independently varying momentum or
a degenerating symplectic normalization would be a different sector.

The same comparison retains intrinsic angular and symmetric shape
moments through the full material variation formula of 0124. Their
centroid/ambient current can redistribute momentum inside a carrier;
it does not change the zero complete-period mean in (10)-(11).

## 4. Executed alternative: genuine Galilean translational data

The exact Euler symmetry gives, for ANY smooth stationary u0,

    u_e(x,t)=u0(x-e(U0+Vt))+e V,
    xi_G=U0+Vt,
    v_G=V-(U0+Vt)·grad u0.                             (12)

Its pressure is the translated stationary pressure. Direct substitution
in Euler or differentiation verifies (12). An actual transported parcel
has delta X=U0+Vt, delta P=M_D V, with unchanged intrinsic spin under
this rigid translation. Thus position and momentum are independent
data at k=0. They are not the two stationary-pattern quadratures in
Section 3, which both have zero momentum at the crossing.

On R3 the boost is nondecaying and belongs to the mean-flow sector, not
the decaying compact Kelvin leaf. On a periodic cell the constant mean
velocity is retained, with the corresponding material/Kelvin periods.
The finite mean mass density gives the free kinetic term rho|Udot|^2/2.
The actual k=0 equation is Uddot=0. No positive spatial stiffness follows
from Galilean symmetry, and one cannot add an STF cage and call it the
dispersion of (12).

Consequently the neutral Bessel pattern and the genuine Galilean
generalized mode are materially different branches. Their coexistence
does not itself prove that a slow modulation dynamically couples them
into an acoustic band.

## 5. Transfer boundary and the next concrete acoustic construction

The free-space Bessel beam is not a compact finite-core perturbation,
and the exact slip cylinder in Section 2 is not a boundary supplied by
the EPS approximation. Those distinctions prevent an automatic transfer
of (7) as a torus dispersion law. A smooth compact narrow-cone packet
can instead be evolved by the exact uniform affine Fourier propagator
of 0124; each Fourier direction has (3). Its laboratory frequency spread
is at least the actual spectral width times the gradient in (5), and
the full pressure/collar terms are retained. A fixed-time actual-core
transfer follows at the scope of 0124's affine kernel estimates, not an
infinite-time spectral identity.

Resolving a sideband |omega_lab|~Omega |K|/kappa over its acoustic period
requires T~kappa/(Omega |K|). A packet whose envelope stays in a fixed
core has a finite wave-number spread. Taking K to zero while discarding
that spread or using a fixed-optical-time error bound is not licensed.
This is a time/space scale issue, separate from the already explicit
symplectic/physical-mass mismatch.

The next representation is the ACTUAL mean translation Jordan sector
of a fixed periodic Euler field, not the carrier's squared frequency.
Let L0 v=-P[u·grad v+v·grad u], and define on a fixed periodic cell

    t_j=-partial_j u,    e_j=constant unit velocity.

Stationary Euler and its derivative imply exactly

    L0 t_j=0,   L0 e_j=t_j.                            (13)

Moreover <t_j>=0 and <e_j>=e_j. Thus this Jordan chain already contains
the correct independent mean momentum and translated microstructure.
For a Bloch perturbation the complete mean equation is

    d_t <v_K>=-i P_K <u tensor v_K+v_K tensor u> K.      (14)

At K=0 the stress of t_j has zero mean by periodic integration by parts;
the stress of e_j also has zero mean if <u>=0. Acoustic order therefore
depends on the FIRST genuine cell-response correction to these modes,
not the zero-order stationary pattern. Equations (13)-(14) identify the
new constructive task: solve the actual small-K Jordan/cell response,
including its nonaffine Euler complement or controlled memory, and
derive the stress residue seen by the mean momentum. The exact chain and
vanishing stress have been executed here; a positive cell-response
coefficient has not been assumed or supplied.

## Verdicts and scientific scope

- Exact m=1 inertial/Bessel Euler/Lin field, true carrier sideband,
  slip-cylinder action and nonzero subparcel centroid: established.
- Interpretation of the zero-frequency crossing alone as a positive
  finite-mass acoustic continuum: refuted for this fixed mode family.
  Mechanism: nonzero KKS with vanishing physical momentum, sign-changing
  laboratory Hamiltonian, and zero complete-carrier mean.
- Galilean/periodic mean translation Jordan chain: established, with
  independent physical k=0 position/momentum and no inferred stiffness.
- Positive acoustic cell-response route: active at its named missing
  actual Euler stress/complement calculation. The candidate space and
  parent objective are not exhausted by the first route's failure.

The oracle differentiates the actual Euler/Lin formula, boundary branch,
centroid, action and mean Jordan identities with canonical CheckLedger.
No numerical stability inference or preselected sound speed is involved.
