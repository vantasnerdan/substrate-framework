# Exact radial localization of a positive-action helical Euler pole

Fix `c,C,rho>0`, integer `m=2`, `k0=-m/c`, and an axial period on
which this Fourier sector is single valued. Set `d=c^2+r^2`,
`h=r e_theta+c e_z`, and `f=C/d`. All norms below are on the SAME
fixed `(m,k0)` solenoidal radial velocity space with physical measure
`r dr`, including the full whole-transverse-space Leray projection.
No radial boundary condition other than regularity and finite energy
is imposed. The positive simple ground pole `i sigma`, its actual
fixed-Kelvin action, and the three interior marker controls are the
0185/0189 inputs. This construction strengthens their finite-time
radial localization to an exact localized-background eigenmode.

## 1. A stationary finite-energy background and an operator-norm estimate

Choose a fixed smooth `chi_cut`, equal to one for `r/R<=1`, zero for
`r/R>=2`, and between zero and one. Define

    f_R=chi_cut(r/R) f, u_R=f_R h,
    p_R'=r f_R^2.

This is exact smooth steady incompressible Euler. Its velocity has
compact radial support, so its FULL background energy per axial period
is finite, not merely a renormalized perturbation energy. Its outer
vorticity and pressure are the ones derived from this profile. In
particular it is not assumed force-free or constant-helical-momentum
in the return annulus.

For fixed azimuthal/axial indices the complete generator is

    L_f v=-P_mk[i(m+ck)f v+2f Jv+f' v_r h],
    J=e_z cross.

There is no radial advection derivative. The Leray projection is the
same norm-one operator for both backgrounds. For `R>=c`, direct
differentiation, not a residual estimate on one selected mode, gives

    ||f_R-f||_infty <= C/R^2,
    ||sqrt(d)(f_R'-f')||_infty
                  <= C(2+||chi_cut'||_infty)/R^2,

and consequently

    ||L_R-L|| <= C[|m+ck|+4+||chi_cut'||_infty]/R^2.       (1)

Thus even the nonlocal pressure difference is bounded in operator norm.
There is no inverse period, pressure truncation, or omitted outer tau
block in (1). The helical-momentum equation at `k=k0` is exactly

    (v.h)_t=-(d f_R)' v_r.                              (2)

The right side is generally nonzero in the annulus.

## 2. The exact pole, not a quasimode-to-spectrum assertion

At the uncut background, the orthogonal decomposition into the
solenoidal perpendicular velocity and `tau h/d` is legitimate at
`m+ck0=0`: each piece is separately solenoidal, and its physical norm
uses `tau/sqrt(d)`, not the unbounded raw tau coordinate. The tau block
is stationary, its coupling into the perpendicular block is bounded,
and the latter is the compact skew operator of 0185's positive radial
problem. Hence the simple positive ground pole `i sigma` is isolated
in the FULL fixed-sector velocity generator, not only in its tau-zero
restriction.

Fix a spectral contour around this one pole, separated from zero and
the remaining fixed-sector spectrum; let `M_Gamma` be the finite maximum
of the uncut resolvent norm there. Choose finite R so that the right
side of (1) times `M_Gamma` is smaller than one. The resolvent Neumann
series and contour integral then give an actual rank-one Riesz
projection and eigenvalue `lambda_R`, with

    |lambda_R-i sigma|+||v_R-v||_2 <= C_Gamma/R^2,        (3)

after fixing the continuous mode normalization. This argument compares
the entire operator. It does not use a small residual as evidence of
an eigenvalue of a nonnormal generator.

The same projection is obtained along `L+s(L_R-L)`, `0<=s<=1`;
these operators are the actual intermediate radial Euler profiles.
This also fixes its phase convention and rank without requiring a
cutoff parameter to be analytic at infinity.

## 3. Exact Kelvin reconstruction, reality and the new action

Initially allow the spectral parameter to be complex. Write
`lambda_R=i sigma_R`, with `sigma_R` close to the positive sigma.
The actual material displacement is

    xi_R=v_R/(i sigma_R)+f_R' v_R,r h/(i sigma_R)^2.     (4)

Here `v_R,r` denotes the RADIAL COMPONENT, not a derivative. The
commutator is `[u_R,xi_R]=-f_R' xi_R,r h`, so (4) satisfies the full
Lin relation. Its divergence vanishes because `h.grad` annihilates
the scalar coefficient in this exact helical sector.

For a general radial f, put `T'=2rf+d f'`. Its vorticity is
`omega=(0,-c f',2f+r f')` in cylindrical components. The radial part
of the difference between `xi cross omega` and
`v+grad pi/(i sigma)` is

    f' tau/(i sigma)+f' v_r T'/(i sigma)^2=0,

using the ACTUAL equation `i sigma tau=-T'v_r`. The other two
components follow directly from the momentum equations. Therefore

    xi_R cross omega_R=v_R+grad pi_R/(i sigma_R),
    P(xi_R cross omega_R)=v_R.                         (5)

This proves fixed-Kelvin compatibility even in the annulus with nonzero
tau; it does not borrow the force-free simplification of the original
profile. All exterior pressure and velocity tails remain present.

The maps in (4) are bounded on physical L2. Equations (1)--(3) and
`||omega_R-omega||_infty=O(R^-2)` imply continuity of the full actual
KKS integral `rho integral omega_R.(xi_a cross xi_b)`. It converges to
the strictly positive coefficient `beta_0=h_m/sigma` of 0185.
Conservation of the exact KKS pairing of an eigenmode and its complex
conjugate excludes `Re lambda_R !=0`, since otherwise their nonzero
pairing would scale by `exp(2 Re lambda_R T)`. Thus `sigma_R` is real
and positive. For the real coordinates with

    a_t=sigma_R b, b_t=-sigma_R a,

write the positive full KKS coefficient as `beta_R`. The actual Euler
Hamiltonian on this phase plane is then

    H_R=h_R I, h_R=sigma_R beta_R>0.                   (6)

This is the localized field's own action, not the old kinetic metric
assigned to a changed mode. The base is an actual equilibrium of the
finite-energy Euler Hamiltonian, and (5) places these perturbations on
its Kelvin leaf. The identity `H=-Omega A` fixes its physical clock
and second variation. No rotating-frame frequency shift is selected.

## 4. The inner and outer radial equations

The exact first-order pressure/material-radial system is useful both
for regularity and observation continuity. Put `F=xi_r=v_r/(i sigma_R)`.
For ANY radial f at the exact helical carrier,

    v_theta=-m pi/(sigma_R r)-(2f+r f')F,
    v_z=m pi/(sigma_R c)-c f' F,
    F'=(-1/r+2m f/(sigma_R r))F
                     +m^2(1/r^2+1/c^2)pi/sigma_R^2,
    pi'=[sigma_R^2-2f(2f+r f')]F
                                  -2m f pi/(sigma_R r). (7)

Only the nonzero frequency occurs in the denominators. In particular
there is no hidden artificial boundary at a local inertial turning
radius. Smoothness at the axis selects the regular solution. For
`r>2R`, f=0, so

    pi''+pi'/r-m^2(1/r^2+1/c^2)pi=0.

The actual finite-energy pressure is the decaying `K_m(mr/c)` branch,
not zero pressure on a cutoff wall. It gives exponential exterior
velocity and displacement tails and finite perturbation action.

For `r<R`, equation (2) gives `tau_R=0` exactly, and hence

    v_r=i m phi_R/r, v_theta=-c^2 phi_R'/d,
    v_z=c r phi_R'/d,
    sigma_R H_m phi_R=m beta(r) phi_R,
    beta(r)=8Cc^2/d^3.                                 (8)

The regular inner solution is the actual one at the shifted eigenvalue
sigma_R. It is not the old global ground function inserted into the
new field. In a fixed interior annulus its pressure and all needed
radial derivatives converge at rate O(R^-2). One way to see this is to
use the nonzero axial Fourier number to control the pressure L2 norm
by its gradient in (5), then use (7) on the fixed annulus. The regular
inner solution and its chosen normalization depend smoothly on sigma_R.
These are local ODE estimates with fixed coefficients and a fixed
observation region, not a claimed global high-Sobolev operator bound.

## 5. One stationary positive material tag and BOTH physical currents

Fix the three smooth interior radial controls whose nonsingular moment
matrix was proved in 0189, within a fixed annulus where `chi_tag>0`
and `I_obs!=0`. For sufficiently large R this support lies entirely in
the UNCHANGED core. Use the actual new phi_R and sigma_R in the rows

    B(r)=r^(m+1), R_row(r)=r phi_R/d,
    J_row(r)=-c^2 r^2 phi_R'/d
                           +2mCc^2 r phi_R/(sigma_R d^2),
    I_R=integral chi_tag r^(m-1)(m phi_R+r phi_R')dr.

Their three-bump matrix converges to 0189's invertible one. For a
predeclared eta>0 and fixed `B_* !=0`, solve the actual linear equations

    integral chi_tag b_R B=B_*,
    integral chi_tag b_R R_row=0,
    integral chi_tag b_R J_row
                           =-eta h_R B_*/(rho pi Lz I_R). (9)

All coefficients, including h_R, come from the localized eigenmode and
its full exterior reaction. Matrix continuity supplies a finite exact
solution `b_R` converging to the original marker control. A sufficiently
small positive epsilon keeps

    w_R=chi_tag(r)[1+epsilon b_R(r)cos(m(theta-z/c))]

between zero and one, with nonzero reference quadrupole. The base label
fraction and the scalar `(x+i y)^m exp(-im z/c)` are stationary on its
support. We select the actual particular Lin displacement (4), not an
arbitrary initial tag perturbation with an unrecorded constant offset.

Direct material integration inside this unchanged core gives the same
FORM of the 0189 equations, with the new mode in every entry:

    theta=c_R b, c_R=I_R/(sigma_R epsilon B_*),
    M_R=h_R/(sigma_R^2 c_R^2)>0,
    G_b=rho pi Lz epsilon integral chi_tag b_R
             [c^2r^2 phi_R'/(sigma_R d)-m f' r^2 phi_R/sigma_R^2],
    S_a=rho pi Lz epsilon integral chi_tag b_R J_row,
    S-G_t=(2mC rho pi Lz epsilon/sigma_R) a
                                     integral chi_tag b_R R_row.

Equation (9) therefore gives the EXACT all-time linear-mode identities

    G_z=eta M_R theta, S_z=eta M_R theta_t,
    G_z(0)+integral_0^T S_z=eta M_R theta(T).             (10)

The spin is the actual `rho integral w[r v_theta+(r^2 f)'xi_r]`;
the moving-position term is retained. For m=2 the reference planar
centroid and momentum vanish by angular symmetry. The linked actual
axial momentum remains `S_z=-c delta P_z`. Thus neither the global
action nor the local physical spin is a fictitious independent rotor.

## Scope and next continuation

The route is established: an exact smooth radially compact-velocity
steady Euler field, finite background energy per axial period, a genuine
fixed-sector optical pole with positive Euler action, and one positive
stationary tag with the complete all-time physical current match.
Its perturbation velocity is not compactly supported; its pressure
return and induced exterior motion are part of the theorem.

The isolating contour is only for the selected `(m,k0)` sector. Other
azimuthal/axial sectors can have transport bands through this frequency.
Neither a bent tube nor a transverse array nor an EPS exterior preserves
that fixed-sector decomposition automatically. No full-space spectral
gap, carrier curvature, finite-bandwidth exact clock, or coupled
translation/rotation continuum is inferred. Those remain distinct
positive construction targets; this is not parent completion.
