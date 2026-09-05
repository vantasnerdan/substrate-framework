# Exact triangular core and the circular optical candidate

This follows the central264-claim validation of0166. The parent separately
registered the ordinary Euclidean m2 candidate after the circular core was
derived; no0155 artifact is changed. The base is v0.176.0. All imports here
are exact unpromoted evidence at their declared scope, not new canon.

## Field, cell and physical signs

Choose b1=lambda(1,0), b2=lambda(-1/2,sqrt(3)/2), b3=-b1-b2 and

    psi=Psi sum_j cos(b_j.r),  Psi>0, lambda>0,
    v=J grad psi, W=-lambda psi, u=(v,W),
    J(x,y)=(-y,x).

The normal cell is the oblique torus dual to b1,b2. The axial coordinate
s has a separately declared fixed period L_s, so a periodic optical
carrier is k=-p=-2pi N/L_s, N a positive integer. No particular axial
period or action normalization is borrowed from the normal wavelength.
The full Euler field has curl u=lambda u, div u=0 and pressure per unit
mass -|u|²/2. Its maximum of psi at r=0 gives

    Omega=3Psi lambda²/2,
    Dv(0)=-Omega J, W0=-3lambda Psi=-2Omega/lambda,
    omega_s(0)=-2Omega.                                  (1)

Small levels of3Psi-psi are invariant core-tube cross sections. The
axial flow and the full pressure do not require a wall on them.

## Exact jets and their first nonradial term

Writing R=|r| and theta=arg(x+iy), direct Taylor expansion gives

    psi=3Psi-3Psi lambda² R²/4+3Psi lambda4 R4/64
          -Psi lambda6 R6(10+cos6theta)/7680+O(R8).        (2)

The radial terms in(2) coincide with3Psi J0(lambda R).
Consequently u agrees through the radial quartic orders with the
negative-amplitude Lundquist column

    u_L=(-3lambda Psi J1(lambda R)e_theta,
                              -3lambda Psi J0(lambda R)e_s).

In particular v's first nonradial term is degree5, not degree3;
W's first nonradial term is degree6, not degree4. The Cartesian harmonic
is R6 cos6theta=x6-15x4y²+15x²y4-y6. This explicit fact, rather than a
generic circular-jet resemblance, makes the optical transfer useful.

## Selection and natural carrier scales

Both the calibrated odd-m candidate inherited as a method from0155
and a standard m2 physical material quadrupole were registered. The
latter is selected here because the circular helicity structure keeps
its unmarked angle row separate from the small reference marking:
xi_+ has angular m-1 and xi_- has angular m+1. It also has an exactly
even reference label, hence zero reference centroid by inversion.
This avoids the order-one elliptic metric of0155; it is a geometric
change, not a renamed canonical angle. The actual observable remains
a registered material shape, not an asserted absolute vortex director.

For k=-p choose

    ell4=2/(lambda p³), delta=sqrt(2lambda/p), cD=Omega delta.

The full pressure/doppler balance, derived in the companion proof,
has oscillator frequency

    omega=(2-m)Omega-(2n+m)cD+O(Omega delta²).

The measured material clock also includes the actual reference and
axial material transport. For n=2,m=2 its first coefficient is

    gamma=2Omega+(1/3)cD+O(Omega delta²),
    p² partial_p² gamma²=Omega² delta+O(Omega² delta²).     (3)

All comparisons are with the small Omega² delta scale in(3).
The actual C6 transverse correction in the mode operator is O(delta6),
and the C6 axial Doppler correction is O(delta7): lambda ell is
proportional to delta^(3/2), while p/lambda is proportional to delta^-2.
Their complete pressure returns still need the fixed-cell estimate;
local jet closeness alone is not that estimate.

## Compatibility boundary for a different acoustic lift

The same normal array with any smooth W=F(psi) is stationary Euler
with the unchanged planar pressure. Its optical leading coefficients
depend on a=|F'(3Psi)|, provided a>0: choose sign k=sign F', replace
lambda by a in ell4 and delta, and retain Omega and omega_s(0) from(1).
This is an explicit coefficient transfer condition, not a proof of its
acoustic response or of constant-curl compatibility. For
W²=C+lambda²psi², F'(3Psi)=3lambda²Psi/W0 is nonzero. construction_review
owns that acoustic comparison. A negative acoustic verdict on the
constant-curl field does not change(1)--(3), and positive optics alone
does not give a positive joint continuum on either field.

Verification design: exact background residual, independent polynomial
jet expansion, full pressure/Lin carrier balance, KKS normalization,
actual moment Laplace derivative, finite tag rank and a wrong-carrier
mutation. No numerical spectrum or comparator has been inspected.
