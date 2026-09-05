# Exact toroidal pressure and Kelvin pullback

This calculation uses the actual 0211 ring, after central activation of
0220. The ambient pressure remains the complete-space Euler pressure.

## 1. Exact Piola variables and action measure

Use coordinates (x,y,z), y=R phi, r=R+x and h=1+x/R. Their physical
map is X=(r cos(y/R),r sin(y/R),z); its positive volume Jacobian is h.
For a physical vector v define its divergence-preserving flux vector

    q_v=h(DX)^-1 v=(h v_r,v_phi,h v_z).                    (1)

Then div q_v=h div_X v. In particular the exact source imported from
0213 is pulled back by (1) BEFORE taking its straight poloidal-zero
component. This is the precise divergence-preserving interpretation
of that import; naive component averaging would omit a curvature row.

For a physical force F its covector is alpha=(DX)^T F. Define

    N=h g^-1=diag(h,1/h,h),   g=(DX)^T DX.

The exact physical Leray response is

    q_v=N(alpha−grad p),
    div(N grad p)=div(N alpha).                           (2)

The pressure domain is x>−R, y periodic of period 2pi R, z real,
with regularity at the physical axis and the full decaying harmonic
condition at infinity. No artificial source-disk boundary is used.
Its kinetic metric and force pairing are

    integral |v|² dX=integral q_v dot N^-1 q_v dx dy dz,
    integral v dot F dX=integral q_v dot alpha dx dy dz.    (3)

Let a=q_u, w=q_omega and q=q_xi. Cross products pull back as

    alpha_(xi cross omega)=(q cross w)/h.                 (4)

Thus the complete KKS pairing has integrand
rho w dot(q_1 cross q_2)/h, with its overall sign fixed by the chosen
orbit convention. The volume and metric factors are not optional.

## 2. Full Kelvin generator in these same variables

Write Pi_N alpha for the divergence-free response in (2). The Lie
bracket of physical fields pulls back to h[a/h,q/h]. Therefore

    A_R q=Pi_N[(q cross w)/h]−h[a/h,q/h]
         =Pi_N[(q cross w)/h]−[a,q]/h
             +(a_x q−q_x a)/(R h²).                      (5)

This is the exact fixed-Kelvin reconstruction generator and contains
the actual pressure. Here a_x and q_x denote components, not
derivatives. For the actual streamfunction phi_R=psi/R,

    a=(−phi_R,z, G(phi_R)/h, phi_R,x),
    w=G'(phi_R) a.                                       (6)

The first and last entries mean minus the z derivative and the x
derivative, respectively. These identities hold through the smooth
outer taper and retain its exact zero factors.

On any fixed-domain differentiable metric family, N=1+epsilon D+
epsilon²E, differentiating (2) gives the complete pressure responses

    Pi_0=P,
    Pi_1=P D P,
    Pi_2=P E P−P D (1−P) D P.                            (7)

P is the flat Leray operator with the chosen axial Fourier frequency.
This is derived by solving each pressure equation, not by dropping
the pressure-gradient pieces. For the raw toroidal chart,
D=x diag(1,−1,1) and E=x² diag(0,1,0).

The chart's far axis/infinity move with R, so (7) alone is NOT a
global remainder theorem for the ring. The actual scalar Green
kernel must supply that domain comparison. Keeping this distinction
prevents a local metric Taylor expansion from silently bounding a
complete-space inverse.

If a=a0+epsilon a1 and w=w0+epsilon w1, its formal first coefficient
on that fixed pressure comparison is the explicit operator

    A1 q=P D P(q cross w0)+P(q cross w1−x q cross w0)
          −[a1,q]+x[a0,q]+a0_x q−q_x a0.                 (8)

Every displayed contribution belongs to the same physical generator.
The last two metric-bracket terms and the inner P in P D P are
precisely the terms lost by treating this as a flat force-only
perturbation.

## 3. Actual material action-angle coordinates

The positive active cross-section has smooth nested stream levels
phi_R=H, including the finite smooth outer vorticity edge. The actual
poloidal velocity does not vanish on that edge; the central elliptic
limit is nondegenerate. On each regular level define

    T_R(H)=R^-1 integral r dl/|grad phi_R|,
    Omega_R(H)=2pi/T_R(H).

Choose theta as elapsed poloidal time times Omega_R, and orient H
outward with decreasing values. The exact cross-section volume
identity is

    dX=(1/Omega_R) |dH| dtheta dy_old,
    dI/dH=−1/Omega_R.                                    (9)

The sign of I is chosen outward so the physical measure is positive.
Remove the nonuniform toroidal advance by y=y_old−R chi(H,theta),
where the periodic zero-mean chi solves

    Omega_R chi_theta=dot phi−average_time(dot phi).

This shear leaves (9) unchanged. In the resulting positive
volume coordinates (I,theta,y),

    dX=dI dtheta dy,
    u=(0,Omega_R(I),W_R(I)),
    W_R=(R G(H)/T_R) integral dl/[r|grad phi_R|].          (10)

All period averages here are time averages on actual material
streamlines. They do not impose constant angular speed on a painted
tag. Smoothness at the elliptic center is handled in Cartesian
action-angle charts; the tail estimate is on a regular outer annulus.

Let g_R be the physical metric in these coordinates, and let
ell=(g_R u)_theta and b=(g_R u)_y be the covariant velocity rows.
The toroidal coordinate derivative is unchanged by the shear, hence

    b=G(H) EXACTLY,                                      (11)

pointwise, not only after averaging. Since volume is one, the curl
equation omega=G'(H)u gives

    −partial_I b=G'(H) Omega_R,
    partial_I ell−partial_theta(g_R u)_I=G'(H) W_R.

In particular

    partial_I average_theta(ell)=G'(H)W_R.              (12)

Equations (9), (11) and (12) are the exact circulation/force-free
identities for this actual curved metric. They are the replacement
for an unjustified claim that the period-averaged profile is itself
an ordinary Euclidean straight column.

## 4. First-order normalization and the real second-order problem

At order 1/R, the centered radial Green solution has only the
cos(theta) core deformation; the geometric Jacobian and toroidal
advance corrections have the same first poloidal order. On each
fixed level their complete first period, advance and shell-volume
averages vanish. Equation (9) and the periodic shear therefore
normalize the background contravariant u and omega to their radial
limits THROUGH FIRST ORDER. This statement is about the coefficients
after the actual coordinate change, not the original Cartesian
profile at a shifted exponential edge.

More explicitly, expansion of the actual scalar equation gives
(-Delta−f_delta'(phi0))phi1=−partial_x phi0. Its first-order source
and outer matching are poloidal order one. The translation-speed
and center-border rows fix the corresponding homogeneous
translation; the mass border leaves no forced radial component.
The outer matching includes its true growing/logarithmic dipole,
so this is not an incorrectly imposed decaying m=1 inverse. On a
fixed regular flux contour the resulting normal displacement and
each first integrand variation are a coefficient times cos(theta).
Their full period integrals therefore vanish. This directly
licenses the first-mean cancellation used above; the second-order
mean is retained, not discarded by that parity argument.

At second order the averaged periods and metrics genuinely change.
Retain them with the exact factors (11)-(12): comparing translated
flat functions pointwise by their relative ratio would fail near
their exponentially flat zeros. Those exact averaged covariant
identities are the input to the next weighted form calculation.
The full nonzero-poloidal pressure return and its domain remainder
are still calculated separately; no first-order selection rule is
declared to have summed that return.
