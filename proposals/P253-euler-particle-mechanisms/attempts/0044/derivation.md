# Physical-domain curved Feshbach calculation

## 1. Exact domain: compact vorticity, global velocity

Let `Omega_delta` be the compact toroidal vorticity of the fixed Cao ring and
`K_delta` the whole-space Biot--Savart map. Define the physical graph domain as
the closure of

    eta=curl(xi cross Omega_delta),  div xi=0,  xi in C_c^infinity,

in

    ||eta||_X=||K_delta eta||_L2+||eta||_H-1
              +||A_delta eta||_H-1.                         (1)

The first term is kinetic energy, the second makes the distributional curl
explicit, and the last is the generator graph norm. If `eta_n` is supported in
the fixed closed core and converges in (1), then it converges in distributions.
Pairing with any test supported outside the core proves that the limit has the
same support. Thus the physical closure contains no passive exterior
vorticity. It retains the global velocity `K_delta eta`, including its harmonic
exterior tail and pressure.

The linearized coadjoint equation preserves this space: differentiating
`(g_t)_* Omega_delta` gives the tangent above, while the linearized Euler flow
is the tangent map of the Euler coadjoint flow. Equivalently, outside the
stationary material support both `Omega_delta` and its derivatives vanish, so
homogeneous exterior vorticity transport preserves zero data. This proves only
support preservation, not that every core-supported distribution is
dynamically accessible.

Consequently the Gallay--Smets essential endpoint produced by passive
vorticity at radii where `Omega(r)->0` is absent from this realization. The
exterior Hodge problem is not absent.

## 2. Independent finite-core response and its sign

For the limiting column, let

    F(r)=integral_0^r W(s)s ds,       Omega(r)=F(r)/r^2.

The actual planar `m=1` Euler pencil on the compact core is

    T(b)=-d_rr-r^-1 d_r+r^-2+W'(r)/[r(Omega(r)-b)]          (2)

with regular origin behavior and the exterior matching condition
`psi'(a)=-psi(a)/a`. Since `W=2Omega+r Omega'`, direct substitution gives

    psi_0=r Omega>0,          T(0)psi_0=0.                  (3)

For `psi=psi_0 f`, integration by parts with measure `r dr` gives

    <psi,T(0)psi>
      =integral_0^a r psi_0^2 |f'|^2 dr.                   (4)

The form boundary term `+|psi(a)|^2` is essential: it cancels
`[r psi_0 psi_0' f^2]_a` because `psi_0'(a)=-psi_0(a)/a`.
Equation (4) proves a simple zero, and compact finite-interval form embedding
gives a positive complement eigenvalue `lambda_*`.

The requested sign is fixed by the convention in (2):

    partial_b (Omega-b)^-1=+(Omega-b)^-2.

Therefore

    <psi_0,T'(0)psi_0>_(r dr)
      =integral_0^a r^2 W'(r)dr
      =[r^2 W]_0^a-2 integral_0^a rW(r)dr
      =-2F(a)<0.                                            (5)

This is not the Cao stationary elliptic inverse. It is the time-dependent
Euler translation response.

For nonzero longitudinal `k`, the exact core pencil and exterior matching give

    d(k,b)=-2F(a)b+F(a)^2 k^2 log(|k|a)+O(k^2)
      +O((|b|+k^2|log k|)^2).                              (6)

The potential DtN is `D/a`, `D=x K_1'(x)/K_1(x)`, but the scalar variable in
(2) obeys `psi'/psi=(1+x^2)/(aD)`. Hence the form coefficient is

    -1/D=1+x^2[log(x/2)+gamma_E]+O(x^4 log(x)^2),           (7)

not `D`. From (5), the scalar implicit-function theorem produces

    b(k)=-(F(a)/2)k^2 log(1/(|k|a))+O(k^2).                 (8)

Moreover, if the relative form perturbation is at most
`C(|b|+k^2|log k|)<lambda_*/2`, then

    ||[Q_1 T(k,b)Q_1]^-1||_(D* -> D)<=2/lambda_*.          (9)

Equations (5), (8), and (9) are the finite-core response and graph-domain
rank-one Riesz input for the straight `m=1` fiber.

## 3. Exact toroidal Hodge/Leray blocks

Use tube coordinates

    X(s,alpha,theta)=(R+a s cos(alpha))e_r(theta)
                     +a s sin(alpha)e_z,
    h=1+delta s cos(alpha),       delta=a/R.                (10)

The orthogonal scale factors are `(a,a s,Rh)`. After multiplying derivatives
by `a`, the exact gradient and divergence in the orthonormal frame are

    G_delta f=e_s f_s+e_alpha f_alpha/s
                 +delta h^-1 e_theta f_theta,              (11)

    D_delta V=(s h)^-1{partial_s(s h V_s)
                +partial_alpha(h V_alpha)
                +delta partial_theta(s V_theta)}.           (12)

For the fixed longitudinal sector `partial_theta=i l`, write
`G=G_0+delta G_1+delta^2 G_2+...` and similarly for `D`. Then

    G_1 f=i l e_theta f,
    G_2 f=-i l s cos(alpha)e_theta f,                       (13)

    D_1 V=cos(alpha)V_s-sin(alpha)V_alpha+i l V_theta,
    D_2 V=-s cos(alpha)D_1V.                                (14)

The scalar Hodge Laplacian `L_delta=D_delta G_delta` has

    L_0=partial_ss+s^-1 partial_s+s^-2 partial_alphalpha,
    L_1=C=cos(alpha)partial_s-sin(alpha)s^-1 partial_alpha,
    L_2=-s cos(alpha)C-l^2.                                 (15)

Let `S_0=L_0^-1` with the whole-space inner/exterior matching conditions. Its
inverse coefficients are

    S_1=-S_0 L_1 S_0,
    S_2=S_0 L_1 S_0 L_1 S_0-S_0 L_2 S_0.                  (16)

Thus the actual Leray projector `P_delta=I-G_delta S_delta D_delta` is

    P_0=I-G_0S_0D_0,                                       (17)

    P_1=-G_1S_0D_0+G_0S_0L_1S_0D_0-G_0S_0D_1,             (18)

    P_2=-G_2S_0D_0+G_1S_0L_1S_0D_0-G_1S_0D_1
        -G_0S_0L_1S_0L_1S_0D_0+G_0S_0L_2S_0D_0
        +G_0S_0L_1S_0D_1-G_0S_0D_2.                       (19)

These are the requested pressure blocks; omitting any of the inverse-Laplacian
terms changes the second-order Feshbach matrix.

For completeness, let `B_delta(V,Y)=nabla_delta_V Y` in the moving orthonormal
frame. The componentwise cross-sectional frame derivative is

    D_alpha Y=(Y_s,alpha-Y_alpha)e_s
       +(Y_alpha,alpha+Y_s)e_alpha+Y_theta,alpha e_theta,

and its longitudinal derivative is

    D_theta Y=(Y_s,theta-cos(alpha)Y_theta)e_s
      +(Y_alpha,theta+sin(alpha)Y_theta)e_alpha
      +(Y_theta,theta+cos(alpha)Y_s-sin(alpha)Y_alpha)e_theta,

so

    B_0(V,Y)=V_s partial_sY+(V_alpha/s)D_alphaY,
    B_1(V,Y)=V_theta D_theta Y,
    B_2(V,Y)=-s cos(alpha)V_theta D_theta Y.                (20)

If the Cao base has `U_delta=U_0+delta U_1+delta^2U_2+...`, define

    N_0v=B_0(U_0,v)+B_0(v,U_0),
    N_1v=B_1(U_0,v)+B_1(v,U_0)+B_0(U_1,v)+B_0(v,U_1),
    N_2v=B_2(U_0,v)+B_2(v,U_0)+B_1(U_1,v)+B_1(v,U_1)
          +B_0(U_2,v)+B_0(v,U_2).                          (21)

The full velocity generator coefficients, including pressure, are exactly

    A_0=-P_0N_0,
    A_1=-(P_1N_0+P_0N_1),
    A_2=-(P_2N_0+P_1N_1+P_0N_2).                           (22)

Equations (18)--(22) also expose a source boundary: Cao Proposition 3.11 gives
an odd first core correction via an explicit linearized Lane--Emden cell
problem and an `O(delta^2|log delta|)` remainder. It does not provide a
convergent `U_2` coefficient in the graph norm required by (22). The exact ring
determines `U_2`; the cited theorem does not yet construct that spectral jet.

## 4. The `m=0,2` audit

The statement “`C_1` shifts `m`” applies only to its curvature factors. The
`i l` pieces in (13)--(14) preserve `m`. The curvature operator in (15) acts on
the translation streamfunction exactly as

    C[psi_0(s)e^(i alpha)]
      ={(psi_0'+psi_0/s)/2}
       +{(psi_0'-psi_0/s)/2}e^(2i alpha)
      =W(s)/2+{s Omega'(s)/2}e^(2i alpha).                  (23)

Thus the `m=1 -> m=0` forcing is the nonzero function `W/2`; it is not removed
by angular orthogonality. Any cancellation in the full `A_1` would have to be
an explicit cancellation among (18), (21), and the base-profile cell term.

The `m=2` zero-wavenumber core is regular. Its scalar form satisfies

    t_m[psi]=t_1[psi]+(m^2-1)integral |psi|^2/r dr
                         +(m-1)|psi(a)|^2,                 (24)

so (4) makes it strictly positive for every `m>=2`.

The `m=0` sector is different. For `s=i sigma`, the exact Gallay--Smets radial
equation is

    H_k u:=(-partial_r partial_r^*+k^2)u
       =(k^2/sigma^2) Phi(r)u.                             (25)

With exterior decay included, the core-localized operator

    B_k=Phi^(1/2) H_k^-1 Phi^(1/2)                          (26)

is positive, compact, and infinite rank. If its eigenvalues are
`mu_n(k)>0`, `mu_n ->0`, then (25) has the accessible Kelvin frequencies

    sigma_n(k)=+/- |k| sqrt(mu_n(k)) ->0.                  (27)

For each nonzero frequency the displacement reconstruction is regular on the
core: for `m=0`, set `xi_r=u_r/(i sigma)` and solve the remaining two algebraic
components from the theta shear and incompressibility equations. The
eigen-equation then gives `eta=curl(xi cross W e_z)`; a solenoidal cutoff of a
vector potential outside the core does not change that curl. Thus these are
not passive exterior-vorticity modes. Equations (23) and (27) prove that the
full translation complement has no neighborhood of zero with a bounded
inverse, even though the isolated `m=1` core complement does.

## 5. Actual leading two-coordinate compression

The singular local part of the kinetic energy of any smooth concentrated tube
is

    H_log=C_log Length[X],
    C_log=rho Gamma^2 log(R/a)/(4*pi),                      (28)

because the near-core velocity is `Gamma/(2*pi d)`. More explicitly, split the
Biot--Savart double integral into centerline separations
`a<|sigma-sigma'|<r_0` and its complement. Taylor expansion of the curve and
unit tangent makes the first region
`rho Gamma^2/(4*pi) integral d sigma integral_a^r0 dq/q`; two shape
variations may be passed through this fixed-separation integral for a fixed
Fourier number. The core and far regions have integrable kernels after the
circulation and centroid are fixed, hence contribute a bounded twice-varied
finite part. This proves the logarithmic coefficient; it does not evaluate the
finite part.
For

    X(theta)=(R+r(theta))e_r(theta)+z(theta)e_z,

expansion of length and axial impulse in the translating frame gives

    (H-cP)_log^(2)=C_log/(2R)
       integral(r_theta^2+z_theta^2-r^2)dtheta.             (29)

The physical filament limit of the KKS form is

    Omega_log(q_1,q_2)=rho Gamma R integral
       (z_1 r_2-r_1 z_2)dtheta.                            (30)

Hence for Fourier number `l` the actual leading centerline compression in
coordinates `(z_l,r_l)` is

    M_l^log=Gamma log(R/a)/(4*pi R^2)
       [[0,l^2-1],[-l^2,0]],                               (31)

with eigenvalues

    +/- i Gamma log(R/a) l sqrt(l^2-1)/(4*pi R^2).          (32)

The factor `l^2-1` follows from the physical impulse subtraction and makes the
rigid `l=1` displacement neutral. This derives the universal leading matrix
from the smooth-tube local energy, rather than importing a hollow-core finite
constant.

Equation (31) is an actual leading compression, but it is not yet a Riesz
reduction of the full ring. The finite `O(Gamma/R^2)` matrix needs the `U_2`
jet in (22), and the Schur complement contains the nonzero `m=0` channel (23).

## 6. Full graph-domain Riesz verdict

Let `tau_l(delta)` denote the magnitude in (32), with
`k=l delta/a`. For fixed `delta`, define the physical axisymmetric separation

    d_0(delta)=dist(tau_l(delta),
      {|k|sqrt(mu_n(k)):n>=1}).                             (33)

On an energy metric in which the uncoupled `m=0` block is skew-adjoint, its
graph resolvent on a contour around `tau_l` is bounded by `1/d_0`. The usual
Feshbach/Riesz argument therefore requires

    ||delta A_1+delta^2 A_2+...||_(D->X)<d_0(delta)/2.       (34)

Neither (33) nor (34) follows uniformly as `delta->0`: the set in (33)
accumulates at zero, while `tau_l(delta)->0`, and (23) makes the coupling
nonzero. A fixed rank-two `P_1` therefore has no earned uniform contour.

The correct representation change is to enlarge `P_1` by all `m=0` modes
within twice the curvature-coupling scale of `tau_l`. That projection is finite
rank for each fixed nonresonant `delta`, but its rank and separation are not
uniform as the ring thins. Exact resonance requires a coupled Hamiltonian
block rather than a two-coordinate matrix. Constructing its high-mode coupling
decay and a nonresonant Cao parameter (or a tame infinite-block normal form) is
the next concrete dependency.

This is not a nonexistence result. It refutes the **uniform fixed-rank
two-dimensional Feshbach route** and replaces the earlier artificial ambient
critical-layer obstruction by the actual internal `m=0` Kelvin accumulation.

## 7. Route verdicts

- Route A is **established** for the physical compact column: support closure,
  simple translation kernel, sign (5), scalar Kelvin continuation (8), and the
  reduced graph inverse (9).
- Route B establishes the exact Hodge/Leray blocks (18)--(22), the nonzero
  `m=0` coupling (23), and the leading matrix (31), but the uniform fixed-rank
  ring Riesz conclusion is **refuted with the mechanism (27), (33)--(34)**.
- Route C is **blocked by the enlarged coupled Evans determinant**: its scalar
  `m=1` factor is now known, but it must include the near-resonant `m=0` modes
  and the Cao second core jet.
- Route D establishes the positive leading centerline Hamiltonian for `l>=2`.
  A full Hamiltonian eigenpair is **blocked by the same internal-mode coupling
  and missing finite Hessian jet**, not by the sign of (31).

The full Cao rotating wave remains open. The strongest next construction is a
finite-`delta` enlarged KKS/Feshbach block with explicit high-`n` coupling
decay and a nonresonance estimate, followed by the solid-torus nonlinear
reconstruction. No quantum or particle identification follows here.
