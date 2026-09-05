# Controlled optical packets on the actual constant-lambda Euler column

## 1. Scope and exact radial domain

Fix lambda,U>0, Omega=lambda U/2 and an optical time window Omega T<infinity.
The actual background is the global smooth Lundquist solution

    u=U J1(lambda r)e_theta+U J0(lambda r)e_z, curl u=lambda u.

It is bounded with bounded derivatives; its total transverse energy is
not claimed finite. The perturbations below have finite transverse energy
in each fixed nonzero axial Fourier fiber. There is no radial wall.
Use the physical axis-following frame with axial speed U and angular
speed Omega, fixed before choosing an optical branch. The laboratory
frequency of a radial mode is omega=kU+mOmega+sigma.

Choose k=-p<0, delta=sqrt(lambda/p), ell=(2/(lambda p³))^(1/4), m=2.
The no-critical-radius argument in analytic-route.md gives |s(r)|>=Omega
for all r when delta is sufficiently small and sigma is near -2Omega,
where s=omega-m O(r)-k W(r). The entire radial interval [0,infinity)
is used. At large r the packet is rapidly decaying; its pressure is
the complete axial-fiber pressure, not a boundary prescription.

## 2. Pressure-resolved oscillator and arbitrary finite-order construction

In the physical rotating/axial frame the exact horizontal and vertical
Euler equations in the total azimuthal-m sector are

    -i s v_perp+2O J v_perp+r O' v_r e_theta+grad_perp p=R_perp,
    -i s v_z+W'v_r+i k p=R_z,
    div_perp v_perp+i k v_z=0.                                  (1)

Here p denotes pressure divided by density. The same equations yield
the exact0140 radial system. The Cartesian representation avoids a
spurious small-r expansion of m²+k²r².

The fiber Leray inverse contains

    (p²-Delta_perp)^-1=p^-2(1-(delta/sqrt(2))Delta_R)^-1.          (2)

For a Schwartz function its expansion through order j has exact remainder
(delta/sqrt(2))^(j+1)(1-(delta/sqrt(2))Delta_R)^-1 Delta_R^(j+1).
The resolvent has L² norm at most1. Polynomially weighted derivative
versions follow by commuting R and derivatives through it. This is a
nonlocal pressure estimate on the full plane.

Let e_+=(1,i)/sqrt(2), e_-=(1,-i)/sqrt(2), J e_+=-i e_+.
The leading transverse profile is e_+ F_n(R), with

    F_n=(R_x+iR_y)^(m-1) exp(-|R|²/2)L_n^(m-1)(|R|²).

Set v_z=i div_perp v_perp/k exactly. Projection of (1) onto e_+ at
the first nontrivial order gives

    (-Delta_R+|R|²)F_n=sqrt(2)c_n F_n,
    c_n=sqrt(2)(2n+m),
    sigma=-2Omega+Omega c_n delta+O(delta²).                     (3)

The opposite-helicity part is removed at first order by adding
e_- (partial_x+i partial_y)²F_n/(4p²). This follows directly by
projecting the pressure gradient; its leading helicity gap is4Omega.
The pressure itself can equivalently be reconstructed from the exact
vertical equation, p=(s v_z+iW'v_r)/k, so no axial-shear term is lost.

For completeness, the higher-order construction is an algebraic recursion,
not an assumed eigenmode theorem. Expand the Bessel coefficients in
lambda ell R=2^(1/4)delta^(3/2)R and expand (2) to the required order.
All resulting coefficient operators map polynomial-times-Gaussian
functions to finite sums of the same type. At each order:

1. Incompressibility fixes the axial coefficient.
2. The nonresonant transverse helicity is solved using its gap4Omega.
3. The e_+ solvability condition fixes the next real sigma coefficient
   by projection onto F_n.
4. On its orthogonal complement in the fixed azimuthal sector, the
   oscillator inverse has denominators4(j-n), j!=n, and fixes the
   remaining coefficient with zero F_n component.

All coefficients are real in the radial pressure/displacement convention.
The fixed-m oscillator eigenvalue is simple, so there is no unsolved
degenerate eigenspace. Each right-hand side is a finite oscillator sum;
the stated inverse and projection therefore construct it explicitly.
This proves existence of a pressure-resolved quasimode to every prescribed
finite order N, with a real frequency polynomial sigma_N. It does not
assert an exact isolated Euler eigenvalue.

Taylor remainders are controlled in weighted Schwartz norms: Bessel
remainders grow at most exponentially in lambda ell R and hence are
integrable against the Gaussian. The exact resolvent remainder in (2)
controls the exterior pressure at the same order. Thus, after increasing
the recursion order if needed,

    ||R_N||_L² <= C_N Omega delta^N ||v_N||_L²,                  (4)

and the same statement holds for any fixed finite number of scaled
radial derivatives, weights and carrier derivatives. At fixed physical r,
p partial_p=-(delta/2)partial_delta+(3/4)R dot grad_R, which preserves
this polynomial-Gaussian class. These bounds provide the physical-carrier
derivatives; they are not inferred from a finite-box eigenvalue.

## 3. Exact Kelvin preparation and actual finite-time Euler evolution

Invert the exact Lin relation to define the displacement xi_N:

    xi_r=i v_r/s,
    xi_theta=i v_theta/s-rO' v_r/s²,
    xi_z=i v_z/s-W'v_r/s².                                     (5)

The divergence commutator for Lin transport and s!=0 imply div xi_N=0.
All fields are smooth at r=0 by their Cartesian construction.

The Kelvin defect is controlled without assuming a finite invariant
ansatz. For a covector in this Fourier sector the material Lie operator
is

    C_s=-i s I+e_r tensor(rO'e_theta+W'e_z),
    C_s^-1=(i/s)I+e_r tensor(rO'e_theta+W'e_z)/s².                (6)

The off-diagonal matrix squares to zero. The full linearized Euler/Lin
identity, modulo exact gradients, gives

    v_N-P(xi_N cross omega_0)=P C_s^-1 R_N.                     (7)

The inverse is uniformly bounded by C/Omega; exact gradients commute
with the covector transport and are removed by P. Thus use the exact
isovortical initial velocity v_K=P(xi_N cross omega_0), not v_N with
its Kelvin defect ignored.

The linear Euler generator about this bounded smooth background has
L² growth bounded by exp(T||Du||_infinity), independent of the carrier.
Duhamel applied to (4),(7), followed by the exact Lin transport equation,
gives actual linear Euler and material-displacement histories differing
from the quasimode by C_(N,T)delta^N in the corresponding normalized
L² observations. Differentiating the exact fiber evolution gives the
first two carrier derivatives as well. Two factors of p can cost at
most delta^-4, since the unscaled derivative of axial advection is
bounded by U rather than 1/p. This loss is retained below.

For genuine finite-amplitude Euler, use the area/volume-preserving
flow of the smooth xi_N for initial labels and the exact coadjoint
one-form preparation. Smooth perturbation well-posedness about the
bounded stationary background supplies the prescribed finite T when
the disturbance amplitude is chosen sufficiently small *after* p,N,T.
Its quadratic remainder is retained; no nonlinear periodic orbit or
all-time invariant Euler manifold is asserted.

## 4. The selected excited packet and positive KKS

Choose n=2,m=2. Its oscillator eigenvalue is12, hence

    sigma_N=-2Omega+6sqrt(2)Omega delta+O(delta²).                (8)

The leading radial displacement is proportional to
r exp(-g r²/2)L_2^1(g r²), g=ell^-2. Its azimuthal displacement is
i times the same profile to leading order; the axial displacement is
smaller by delta^(1/2). The exact full KKS includes both axial and
azimuthal base vorticity. Its leading term is a positive constant times

    rho Omega Lz integral_0^infinity f_2(r)² r dr>0.             (9)

All remaining terms are relatively O(delta), so beta>0 for small delta.
There is no sign selection by frame winding. The laboratory carrier
kU+2Omega is retained; the fixed-frame intrinsic sigma_N<0 gives a
positive two-phase angular action, up to the computed quasimode remainder.
The actual Euler evolution preserves this initial KKS pairing exactly.
The next file uses its actual material-angle row to derive the physical
scalar action and the stronger, positive physical carrier curvature.

Global EPS/toroidal transfer is a separate parent construction. This
file is the whole-plane Lundquist Fourier-fiber/finite-time result,
not a claim that local field closeness alone controls a global pressure
operator or that the column is already a closed knotted tube.
