# A positive-action noncritical Euler mode on the actual annular supplier

This constructs an actual mode of the constrained straight-column limit
in 0198, including its axial shear and full radial pressure. It is not
an 0185 constant-pitch mode transferred by a local approximation. The
closed thin-ring and array operators remain different downstream objects.

## 1. Admissible profile and the exposing first candidate

Choose positive finite `r0,a,W0`, with `r0>2a`. Let W be a smooth
nonnegative bump supported in `(r0-a,r0+a)`, with a unique nondegenerate
maximum W0 at r0. For example its interior is

    W(r)=W0 exp[1-1/(1-((r-r0)/a)^2)].

Set

    Omega(r)=sqrt(2) W(r)/r,
    chi(q)=W(sqrt(q))/(sqrt(2)sqrt(q)),
    u=r Omega e_theta+W e_z,
    Z=2Omega+r Omega', Phi=2Omega Z, p'=r Omega^2.       (1)

The flat zero extension is smooth, since the support is separated from
the axis. This is exactly the assigned pressure-localizable profile
`chi(r^2)[2r e_theta+sqrt(2)r e_z]`. Its absolute energy per axial
period is finite. At r0, `W'=0` implies `Z=Omega`, so

    Phi(r0)=2Omega(r0)^2>0.                            (2)

The first comparison considered a planar annular-vorticity patch with
inner/outer radii a0,b0 and constant vorticity `2Omega0` between them.
Writing `q=a0/b0`, the exact two-edge displacement matrix is

    omega/Omega0 = [[1,-q^(m-1)],
                    [q^(m+1),m(1-q^2)-1]].            (3)

For m=2 its discriminant is identically zero; it is not a simple
positive pole. For m>=3 and sufficiently small q its inner branch
is a useful positive-action PLANAR edge wave, but `omega=m Omega(r)`
crosses the active annulus. The compulsory axial equation then gives
`v_z=-W' F`, with `F=v_r/[-i(omega-m Omega)]`. On that constant-
vorticity annulus, `W'=Omega0(1+a0^2/r^2)/sqrt(2)` is nonzero.
The planar branch cannot therefore be promoted to a smooth full 3D
eigenvelocity by omitting its axial response. This is a route-scoped
critical-layer defect, not a no-go for the admitted profile class.

The following full-pressure construction avoids that mechanism by
putting the actual frequency strictly above the whole fixed-sector
advection range, rather than pinning a radial boundary.

## 2. Full radial Euler and its real quadratic pencil

Fix integer m=2 and a positive axial carrier k. Use the physical phase
`exp(i m theta+i k z-i omega T)` and define

    q_k=m Omega+k W, s=omega-q_k,
    K_r^2=k^2+m^2/r^2, g=2m Omega/r.

The four exact Euler component equations are those of 0198 with its
time sign converted to this convention. Set the actual material radial
displacement `F=v_r/(-i s)`. Where s is nonzero, they give

    v_r=-i s F,
    v_theta=m P/(r s)-Z F, v_z=k P/s-W' F,
    F'=-(1/r+g/s)F+K_r^2 P/s^2,
    P'=(s^2-Phi)F+(g/s)P.                             (4)

Both Z and W' have been retained; their cancellation in the first
equation for F uses the actual derivative `s'=-m Omega'-k W'`.
For `f=rF`, equation (4) is equivalent to

    P=s(s f'+g f)/(r K_r^2),
    Q_omega[f]=integral_0^infinity
       [(s f'+g f)^2/(r K_r^2)+(s^2-Phi)f^2/r]dr=0.    (5)

The differential Euler equation is the first variation of (5), not
the single integral identity alone. Indeed its first-order expression
is `P'-(g/s)P=(s^2-Phi)f/r`; integration by parts gives (5) with
the regular/decaying endpoint conditions. This is the complete real
selfadjoint radial pressure pencil at each real noncritical omega.

Use the radial Hilbert space `L2(dr/r)` and form domain with

    integral [|f'|^2/(r K_r^2)+|f|^2/r]dr<infinity.     (6)

Near zero its derivative coefficient behaves as `r/m^2`; the finite
energy regular solution has `f~r^m`. Outside the compact velocity
support, (4) gives the actual pressure equation

    P''+P'/r-(m^2/r^2+k^2)P=0.

Thus the physical inner pressure is `A I_m(k r)` and the outer pressure
is a multiple of `K_m(k r)`. A regular solution and a decaying solution
are matched through the flowing annulus; neither cutoff edge is a wall.
The induced perturbation velocity is not compactly supported.

## 3. An exact noncritical variational crossing

For sufficiently large finite k, q_k has its unique global maximum
at `r_k=r0+O(k^-1)` and `Phi(r_k)->Phi(r0)`. This follows from the
nondegenerate W maximum and compactness away from it. Write
`qmax(k)=q_k(r_k)`. Choose fixed

    delta0=Omega(r0)/2,
    delta1>sqrt(max_r max(Phi(r),0))+delta0.

For every omega in `[qmax+delta0,qmax+delta1]`, s is bounded below
by delta0 on the WHOLE radial domain. In particular there is no
particle critical layer in either moving edge or pressure tail.

The powers below use fixed core length and velocity units; dimensionally
the trial width is `ell_ref (k ell_ref)^(-3/4)`.
At the lower endpoint choose a smooth test supported in width
`ell=k^(-3/4)` about r_k and normalize it in `L2(dr/r)`. Fixed profile
constants give

    qmax-q_k=O(k ell^2)=O(k^-1/2) on its support,
    integral |f'|^2/(r K_r^2)=O(k^-1/2),
    Q_(qmax+delta0)[f]
                   =delta0^2-Phi(r0)+O(k^-1/2)<0.      (7)

At the upper endpoint `s^2>Phi` everywhere, so the complete form is
strictly positive. The frequency monotonicity needed between these
endpoints is also quantitative. Direct differentiation gives

    partial_omega Q[f]
      =2 integral s[|f'|^2/(r K_r^2)+|f|^2/r]dr
                         -integral (g/(r K_r^2))' f^2 dr.

Since g is smooth and supported in a fixed annulus separated from
zero, `r |(g/(r K_r^2))'|<=C_m/k^2`. Therefore, after increasing
the fixed carrier lower bound if necessary,

    partial_omega Q[f]>=delta0 integral
                        [|f'|^2/(r K_r^2)+|f|^2/r]dr>0. (8)

The selfadjoint Sturm operator associated to (5) has positive leading
coefficient. Outside a compact interval it is the positive free
radial operator with spectral bottom omega^2; localized coefficient
perturbations have compact resolvent difference by one-dimensional
elliptic regularity and local Rellich compactness. Equivalently, the
regular/decaying radial matching gives the same discrete sub-threshold
spectrum. Hence a negative lowest Rayleigh value in (7) is a discrete
ground eigenvalue. Continuity and (8) make it cross zero exactly once
before the positive upper endpoint. The zero lies strictly below the
essential threshold and has a strictly positive, simple ground radial
function f. Its frequency derivative is nonzero by (8).

This constructs an ACTUAL noncritical simple Euler mode, with

    qmax(k)+delta0<omega(k)<qmax(k)+delta1,             (9)

not a numerical or quasimode spectrum. The full fixed-sector operator
has no radial advection; outside its advection range, elimination of
the nonradial components reduces its Fredholm/matching problem to
(4). Thus the simple matching root is an isolated pole of that full
fixed-(m,k) Euler operator. No assertion about other poloidal sectors
is involved.

## 4. Actual Kelvin displacement, full KKS and positive laboratory action

The complete Lin displacement associated to (4) is

    xi_r=F,
    xi_theta=i[m P/(r s^2)-2Omega F/s],
    xi_z=i k P/s^2.                                   (10)

Equivalently `xi=v/(-i s)+(rOmega' e_theta+W'e_z)v_r/(-i s)^2`.
Its divergence is zero: the derivative of s cancels the shear term.
The full Euler equations give
`xi cross curl u=v+grad[P/(-i s)]`, with the derivative of s retained
inside the gradient. Hence `P_mk(xi cross curl u)=v`. This is the
actual Kelvin orbit, not a second independent circulation oscillator.

Take the real and imaginary mode columns in this fixed lab clock.
Their exact KKS coefficient is

    beta=2pi rho Lz integral r[
          (m Z/r+k W')F P/s^2-Phi F^2/s]dr.            (11)

This expression includes the axial vorticity/shear contribution. To
reduce it, use `m Z/r+k W'=g+q_k'`,
`P'=(s^2-Phi)F+gP/s`, and `s'=-q_k'`. Integration by parts cancels
the q_k' pressure term, leaving

    beta=-2pi rho Lz integral
          [s f^2/r+(s f'+g f)f'/(r K_r^2)]dr
         =-pi rho Lz partial_omega Q[f]<0.             (12)

The endpoint terms vanish by the actual regular/decaying pressure
solutions. This proves nonzero wave action and its sign without a
rotating-frame substitution. In real coordinates
`a_T=-omega b,b_T=omega a`, the actual Euler Hessian is

    H_phase=-beta omega I>0.                          (13)

The frequency omega is positive by (9). Positivity is not borrowed
from a rotating observer, an ideal MHD energy or the old helical field.

## 5. Actual carrier branch and limits of transfer

The coefficients of the full radial problem are analytic in `(omega,k)`
in this noncritical region. The simple ground root and its normalized
mode therefore form an actual analytic branch near every selected
sufficiently large finite carrier. Its pressure, displacement, KKS and
fixed interior observations have all finite parameter derivatives needed
for a second carrier jet. The positive signature and strict separation
persist locally. The mode action is finite per fixed reference axial
period, or as a phase-averaged Fourier-fiber density; a monochromatic
mode is not called a finite-total-energy R3 packet.

Moreover `omega(k)=k W0+O(1)` by (9). The globally simple large-k
ground branch cannot be constant; there are arbitrarily large finite
carriers with `omega'(k)>0`. Analyticity makes zeros of its nontrivial
derivative isolated locally. Thus a nonzero group slope can be selected
by an exact existence argument rather than numerical sampling. No
sign of omega'' or of a coherently averaged curvature is inferred here.

The following file constructs a stationary physical material angle
and its complete displacement/spin current through carrier two. This
activates a real optical supplier rather than stopping at a radial norm.
The mode is noncritical only in its specified m sector. Across ALL
poloidal harmonics, the bands `m Omega+k W` can cover its frequency.
Neither a closed-ring pole nor an array/Bloch interaction follows
automatically from local coefficient convergence in the 0198 geometry.
Actual finite-time transfer and its conditioned physical observations
remain available next routes; no parent completion is claimed.
