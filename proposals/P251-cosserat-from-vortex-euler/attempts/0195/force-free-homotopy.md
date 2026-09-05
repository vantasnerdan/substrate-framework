# A global smooth force-free endpoint with the same meridional core

This executes Candidate B's geometry after the small-swirl construction.
The missing optical-mode/current continuation is retained explicitly.

## 1. Earn the square-root regularity before defining the field

The formal primitive G(phi)=sqrt(2 integral_0^phi f(s)ds) is positive
for phi>0. Smoothness of an arbitrary square root at a zero is not
assumed. The 0186 construction leaves the smooth flat boundary taper
free, and the finite-radius border, isolated optical contour and tag
Jacobian all have open nonzero margins after its fixed profile is chosen.
Use the following explicit positive boundary-tail repair if needed.

For a sufficiently small positive delta, set

    f_mod(t)=A_delta t^(-2) exp(-2/t), t>0,
    f_delta=(1-chi(t/delta)) f_mod+chi(t/delta) f(t),

where chi=0 for t<=delta and chi=1 for t>=2delta, smooth with values
in [0,1], and A_delta=exp(-1/delta)>0. Extend by zero for t<=0.
The function f_delta is nonnegative, flat at zero, positive for t>0,
and unchanged on the entire fixed inner plateau. For every fixed k,
its C^k difference from f tends to zero as delta tends to zero: all
derivatives of f vanish faster than each fixed power at zero, and
the explicit exponential dominates every cutoff derivative delta^(-j).
It need not be monotone in the tiny transition interval; none of the
open fixed-profile continuation estimates relies on newly asserting
monotonicity there.

The radial mass-bordered logarithmic equation of 0186 continues by its
invertible radial derivative to this small law change. Its radial source,
smooth streamfunction and finite-radius border remain close. The full
optical contour and tag Jacobian, if used at the no-swirl starting point,
continue by their fixed-profile graph/finite-dimensional margins. Thus
this is a selected actual nearby seed, not an edit of frozen 0186.
Fix delta before the ring radius. Below call this repaired law f again.

Define G=0 for phi<=0 and

    G(phi)^2=2 integral_0^phi f(s)ds for phi>0.       (1)

On 0<phi<delta the integral is exactly
A_delta exp(-2/phi)/2, and hence G=sqrt(A_delta)exp(-1/phi).
Both G and G' are C-infinity and flat at zero. On phi>0 ordinary
positive square-root calculus applies; the full G is C-infinity.
In particular G G'=f globally, with no boundary sheet or singular factor.

## 2. Uniform finite-parameter global Green construction

For eta in [0,1] prescribe the ACTUAL laws

    F_eta(psi)=R eta G(psi/R),
    B_eta'(psi)=-(1-eta²)f(psi/R)/R.

The exact stationary Euler equation is

    -Delta*phi=[(1-eta²)r²/R²+eta²]f(phi).           (2)

In the full 0186 stream kernel its source is

    q_R,eta(r,phi)=[1-eta²+eta² R²/r²]f(phi).        (3)

Use the same equation and physical mass/center borders (3) of
swirling-ring.md, with (3) replacing q_epsilon. On the fixed meridional
disk r=R+x, the bracket in (3) is positive and differs from one by
O(1/R), UNIFORMLY for eta in [0,1], together with every required fixed
spatial derivative. The renormalized full stream kernel differs from
the planar logarithmic kernel by O(log R/R). The entire bordered map
therefore converges, uniformly in eta, to the SAME invertible planar
mass/translation border, not to a different unknown problem at eta=1.

The parameter-uniform small-perturbation IFT produces one solution family
for all eta in [0,1] at sufficiently large finite R. The streamfunction,
core boundary and meridional velocity converge to the selected radial
seed, uniformly in eta. This construction does not require iterating a
local eta-continuation through unproved intermediate invertibility.

The source remains positive on its core and zero outside. Its actual
mass and center are fixed. The translation speed and outer threshold obey

    U_R,eta=Gamma log R/(4pi R)+O(1/R),
    mu_R,eta=3 Gamma log R/(8pi)+O(1),               (4)

uniformly in eta. To check the first coefficient, expand the TRUE kernel
and source. The kernel contributes log R times (x+3y)/(4pi R), while
the extra source factor contributes -2eta² y/R. Its new log term has
only a source first moment, fixed to zero by the center border. The
target x coefficient and its translation-cokernel projection are unchanged.
The nonlinear errors are O(log² R/R²), as in 0186.

The global extension and maximum-principle argument now give actual
global steady fields, with no outside positive island. On the axis and
throughout the exterior, phi is negative and F_eta is exactly zero;
the azimuthal velocity F_eta/r has no axial filament. The flat G and f
make velocity and vorticity globally smooth across their compact toroidal
support. The far field is the actual uniform -U_R,eta e_z. A regular
streamfunction gauge is used as in 0186.

At eta=1, the displayed field is genuinely GENERALIZED FORCE-FREE:

    curl u=G'(phi)u, B'=0.                          (5)

This is an actual smooth global ring, not just a local toroidal Dirichlet
solution or an approximation by a field with a different Bernoulli law.
It improves the source-only geometry interface in force-free-comparison.md.

## 3. Nonzero periodic core and actual endpoint twist

The inner plateau has f=f_c=2Omega_c>0, and the selected planar seed has
V_theta=Omega_c s there. Its streamfunction decreases outward with
phi_s=-Omega_c s. At eta=1 the leading toroidal core speed is
W(s)=G(phi(s)), with

    W'(s)=-2Omega_c² s/W(s),
    W(s)^2=W(0)^2-2Omega_c² s².                    (6)

Choose a sufficiently small fixed inner disk, where W stays positive.
These are exact planar-limit relations following from G G'=f; a large
constant axial offset is NOT an independent parameter in (1).

For the actual thin ring, let I be the physical meridional volume action.
Writing y=I/R, its limit is s²/2. Uniform smooth core convergence gives

    Omega_R(Ry) -> -Omega_c,
    R V_R(Ry) -> W(s),
    rho_R(Ry)/R=Omega_R(Ry)/(R V_R(Ry)) -> -Omega_c/W(s),

including one y derivative on this inner disk. Consequently

    d rho_R/dI -> -2Omega_c³/W(s)^3 !=0.            (7)

The flux-action derivative is rho_R'/V_R and is also nonzero. Since
vorticity in (5) is pointwise parallel and nonzero there, its trajectory
return map is the SAME map: both optical geometric tube fields share
the nonzero twist, not merely two independently chosen tori.

The core toroidal speed tends to W(0)>0, giving a genuine closed
nonzero streamline/vortex line. To obtain nonresonance without any
integer-arithmetic assumption, retain the continuous radius parameter.
Differentiating the full kernel/border gives profile errors with
R derivatives O(log R/R²); this follows either by differentiating the
displayed logarithmic kernel coefficients or the same near/far estimates.
Thus its actual core rotation number obeys

    rho_R(0)=-(Omega_c/W(0))R+O(log R),
    d rho_R(0)/dR=-Omega_c/W(0)+O(log R/R).          (8)

For sufficiently large R this derivative is nonzero. Select a radius
with Diophantine core rotation inside the open large-R existence region,
and select an inner boundary with Diophantine rotation using (7).
Near this inner core f is constant and G is analytic on phi>0, so
the return map is analytic there. The same flux-form invariant-circle
persistence argument already licensed in swirling-ring.md applies.
The core is an unknotted circle, not an arbitrary prescribed knot.

## 4. The precise optical remainder, not a hidden small perturbation

The meridional background stays close throughout the homotopy, but its
toroidal velocity grows from zero to W(s)=O(1) in fixed core units.
For n/R near a fixed nonzero optical carrier p, the transport term
i n V becomes i p W, also O(1). The vorticity's poloidal derivative
response changes by O(1) as well. Thus meridional closeness does NOT
make the full endpoint Euler resolvent close to the no-swirl one.
No positive optical pole has been continued all the way to eta=1 here.

The exact endpoint retains the full active-support operator with
u=Omega partial_theta+V partial_varphi and
omega=G' u. All bands m Omega+n V and the full inverse-curl response
belong in its spectral problem. A separated positive endpoint pole or
an actual continuous branch with controlled collisions is the next
operator construction, together with its physical material observable.

Also, G(phi)>0 throughout the positive core. The inner-supported-swirl
stationary tag repair is unavailable at this endpoint. Its finite-arc
tags obey the actual transport in operator-and-current.md. The current
cannot be imported from a stationary tag on the eta=0 field. A halo-
supported flux profile that increases W(0) while controlling the complete
outer Euler response is a plausible next candidate, not an assumed
resolvent estimate or an independently supplied axial offset.
