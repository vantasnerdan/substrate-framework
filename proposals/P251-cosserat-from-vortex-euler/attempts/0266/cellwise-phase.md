# Cellwise phase removes the leading mixed row; the full transverse gain

This executes candidate B after schema activation. Sources and fields are
those of the frozen 0250/0260 ring construction. The periodic application
assumes the actual compact assembly from 0265; the isolated fixed-ring
application uses its actual whole-space projector. The parent geometry
input remains separate.

## 1. An exact admissible change of preparation

Let xi_0 be one of the compactly supported isovortical generators on the
chosen annulus in a reference cell centered at X_c. Define its periodic
Bloch envelope by

    xi_K(x)=exp(-iK dot(x-X_c)) xi_0(x)                 (1)

on that cell and extend periodically. Its support has a positive gap
from the cell boundary, so (1) is smooth as a periodic envelope. Directly,

    div_K xi_K=exp(-iK dot(x-X_c)) div xi_0=0.           (2)

Thus this is an exact permitted divergence-free Euler/Lin initial
generator. A compact vector potential is not required to justify (2).
The physical generator in cell m is exactly
exp(iK dot(X_c+Lm))xi_0(x-Lm): its slow phase is constant across each
ring. Translation covariance between cells remains the Bloch covariance.

At initial time the local material preparation therefore has no internal
K dependence. In whole space P_K=e^(-iKx)P e^(iKx), so its physical
velocity and subsequent physical Euler/Lin evolution are the same fixed
single-ring solution up to that constant phase. In the compact periodic
assembly, the image-pressure kernel between separated cells is smooth.
The local transported parametrix is the same; all image/mean differences
in a retained finite-window observation have arbitrary-order oscillatory
decay on the high-action carrier. The full projector and its ray-wise
mean are retained before this estimate. Choosing the carrier after the
finite K and normalized source costs makes these remainders o(K^2),
including the required time derivatives. This does not assert exact
decoupling of the periodic pressure.

Consequently the leading local configuration and moment histories for
this family are independent of K. The mixed first moment caused by
exp(iKx) in the old microscopic source has been removed by an exact
source choice, rather than discarded from a Taylor series.

## 2. The acoustic first-K optical row

Use the n=0, m=+1 or -1 acoustic combination whose actual leading
G_z vanishes, and whose covariance tilt is zero by axisymmetry.
Its local spin is generally nonzero and axial; it is not set to zero.
At the reference stationary tag center X=0, the material mean velocity
is zero because the positive action tag is invariant and bounded.

For a collapsed material observable O, differentiation of its Fourier
weight gives

    delta[e^(-iK X)O]=delta O-i(K dot delta X)O_0.       (3)

Equation (1) removes any other leading K derivative of delta O. For
the literal covariance-angle chart the unperturbed variation angle is
zero. Alternatively, applying its projector to the scalar phase times
Q_0 gives P_n Q_0 n=0. Thus (3) has no first-K angle row. The linear
displacement current G has zero unperturbed value, so it also has no
centroid-phase term.

The background mechanical spin is S_0=s_0 n, and the acoustic centroid
variation is delta X=x_0 a_n n. Here a_n is the scalar source amplitude
in that realization. Its remaining term in (3) is proportional to

    n (K dot n) a_n.                                    (4)

Choose a_n proportional to D dot n for the desired polar acoustic
vector D. For proper whole-state rotations the average of (4) is zero:
integral_S2 n_i n_j n_l dn=0. Reflected frames preserve this vanishing.
The zero-K local spin response, being an axial output of a polar input,
cancels under the reflected law as already established in 0250.

Hence the actual leading acoustic-to-(theta,G,S) block has c1=0 for
this source family. The current improvement involving curl U begins
at order three for this acoustic correction, since its hybrid U starts
at order two; it does not restore an order-one row. The optical-to-
hybrid odd block may remain nonzero. Its Schur product with the
vanishing c1 is zero at order two. Lower carrier/pressure terms are
retained in the finite-window error budget before taking this limit.

## 3. Evaluate the full oblique material tensor

It remains to compute the acoustic diagonal using its full tensor.
Let r=R+s cos(theta), z=-s sin(theta), x=r e_r+z e_z,
e_theta=-sin(theta)e_r-cos(theta)e_z, and

    u=V e_theta+W e_varphi,
    xi=C[V e_theta+W_a e_varphi] f,
    W_a=W-gamma m r,  C=N/(lambda nu),
    f=A(s)e^(iNI+im theta),  nu=m V/s.                   (5)

The exact G ratio gamma is used in the actual preparation. Its circular
leading value is W/(2mR), so W_a and its scaled derivatives stay bounded
as R is chosen large. Terms omitted from (5) are one carrier order
lower. The leading material velocity variation is Du xi, not v. At a
fixed s, after dividing out f,

    Du xi=C[(V^2/s)(-cos(theta)e_r+sin(theta)e_z)
       +(W_a/r)(-V sin(theta)e_varphi-W e_r)].           (6)

Substitute (5)-(6) into the literal material tensor

    delta B_ijl=integral rho chi [
       (Du xi)_i x_j x_l+u_i(xi_j x_l+x_j xi_l)].        (7)

The source is physically K independent, so (7) is now the actual
second-moment contribution in 0241; the removed first-moment source jet
is not secretly added back. Centroid corrections in the axis row have
the previously checked zero baseline first-momentum contraction. The
whole tensor and remaining centroid phase are retained in the symmetry
calculation below.

With one common full theta/varphi integration, its even uniaxial entries
are

    b_parallel=delta B_zzz=i pi^2 m C V^2 s/2,
    b_perp=delta B_zxx=delta B_zyy
                          =i pi^2 m C V^2(R^2/s+3s/4),
    b_mixed=delta B_xzx=delta B_yzy=-b_parallel/2.        (8)

For example the mixed entry's complete varphi integral is
C pi V^2 s sin(theta)(3cos(theta)^2-1), exactly minus half
the zzz integrand after its full varphi integral. The two W W_a terms
cancel in that entry. This is why the oblique output cannot be inferred
from the axial scalar alone. The chiral mixed entry B_yzx is O(R),
whereas b_perp is O(R^2). It does not change the leading large-R
transverse coefficient; its actual parity and lower terms remain in
the finite-R source map.

Writing K_n=K dot n, the even uniaxial contraction is

    B:KK=n[b_parallel K_n^2+b_perp|K_perp|^2]
                                      +2b_mixed K_n K_perp.

After the transverse projection, (8) gives

    P_K(B:KK)=P_K n [b_perp|K_perp|^2
                                      +2b_parallel K_n^2].       (9)

The factor two on the parallel term differs from the axis-output-only
formula. The full law, rather than a realization-by-realization scalar
inverse, gives an especially simple robust gain. For D transverse to K,

    E[(D dot n)P_K(B:KK)]
                  =|K|^2 (4b_perp+2b_parallel)D/15.     (10)

This follows from E[n_i n_j]=delta_ij/3 and the fourth sphere moment;
equivalently E[n_x^2 n_z^2]=1/15. Thus no division by a directionwise
gain near n parallel K is needed.

The common action weight, density, sparse-tag coefficient and envelope
in (8) are restored before normalizing. For m=+/-1,
m C V^2=N s V/lambda has a fixed positive real sign on the small
Bessel annulus, after one common complex phase. Choose a nonnegative
envelope there. Equation (10) then has a nonzero coefficient. The
hybrid acceleration multiplies this by the nonzero clock factor and
the fixed -1/2 momentum/positive density factor from 0241. Both real
parities are supplied by the conjugate preparations.

The exact finite-R source has smooth core convergence. After dividing
the complete tensor by R^2, (8)'s b_perp coefficient has a nonzero
limit, while b_parallel,b_mixed and the chiral entries divided by R^2
vanish. Consequently the actual averaged transverse coefficient in
(10) remains nonzero for one sufficiently large **fixed** R. This
uses a rescaled full determinant with a strict margin, not a claim
that every small unscaled oblique coefficient keeps its sign.

## Result and next checks

The cellwise phase route makes the leading lower-left acoustic mixed
row vanish and gives a nonzero full averaged transverse quadratic
gain by (8)-(10). On the common acoustic/optical band the normalized
joint leading block is consequently triangular with invertible
diagonal blocks. Actual lower carrier and periodic-image pressure
terms are then controlled after choosing K; their effect on the
Schur coefficient is retained rather than called automatically small.

The exposing symbolic check must differentiate (7) in Cartesian
components with one common angular measure and evaluate the sphere
moment in (10). Independent review of this new source-family and
centroid-row transfer remains to be performed. The compact stationary
field and final same-field normalization remain separate parent inputs.
