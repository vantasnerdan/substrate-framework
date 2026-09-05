# Exact affine inner current and its complete radial mass border

This selects candidate A. The new seed is a smooth truncated Lundquist
seed, not a claimed small C1 change of the old constant-vorticity plateau.
It is built in the exact current-law class of 0195 and carries its own
full-plane inverse. No numerical profile or spectrum is used.

## 1. Smooth current with an exactly affine inner interval

Fix a physical inverse length lambda>0 and a positive stream amplitude
A. Choose a small positive stream-value thickness delta; its dimensionless
small parameter is delta/A, with A and lambda fixed in every asymptotic
estimate below. Let

    e(t)=exp(-1/t) for t>0, e(t)=0 for t<=0,
    h(t)=e(t)/(e(t)+e(1-t)) for 0<t<1,

with h=0 below zero and h=1 above one. This is smooth, nondecreasing,
flat at both ends, and h(1-t)=1-h(t). In particular its integral on
[0,1] is exactly one half. Set

    H_delta(phi)=integral_0^phi h(s/delta)ds for phi>0,
    H_delta(phi)=0 for phi<=0,
    G_sigma(phi)=sigma lambda H_delta(phi), sigma=+1 or -1,
    f_delta(phi)=G_sigma G_sigma'
                =lambda^2 H_delta(phi)h(phi/delta).          (1)

Every function in (1) is C-infinity, flat at phi=0. For phi>=delta,

    H_delta=phi-delta/2,
    G_sigma=sigma lambda(phi-delta/2),
    f_delta=lambda^2(phi-delta/2).                          (2)

Thus G is literally affine on an open positive-flux interval, not just
approximately linear in radius. Its primitive satisfies exactly
`G_sigma^2=2 integral_0^phi f_delta`. There is no independently appended
axial offset. Moreover

    Q_delta(phi)=f_delta'(phi)
      =lambda^2[h(phi/delta)^2
                    +H_delta(phi)h'(phi/delta)/delta]>=0.  (3)

On the taper, 0<=H_delta<=delta, so Q_delta is bounded by
lambda^2(1+||h'||_infinity), uniformly as delta decreases. It equals
lambda^2 on the affine interval and zero on the exterior. Higher
smooth norms may depend on the fixed delta; they are not used as
uniform bounds when subsequently choosing the large ring radius.

## 2. Actual radial seed and its thin flux taper

Solve the regular radial equation

    phi''+phi'/s+f_delta(phi)=0,
    phi(0)=A+delta/2, phi'(0)=0.                            (4)

As long as phi>=delta its EXACT solution is

    phi(s)=delta/2+A J0(lambda s),
    V(s)=-phi'(s)=A lambda J1(lambda s),
    Z(s)=f_delta(phi)=A lambda^2 J0(lambda s),
    W_sigma(s)=G_sigma(phi)=sigma A lambda J0(lambda s).     (5)

Let j=j_(0,1), the first positive zero of J0, and s_*=j/lambda.
On (0,j), J0>0 and J1>0: the identity (x J0')'=-x J0 gives
J0'<0 there. J1(j)>0 also follows, since a zero derivative at a
zero of J0 would violate uniqueness of the Bessel equation.

The inner endpoint a_delta is uniquely determined by
J0(lambda a_delta)=delta/(2A). The inverse-function theorem at j gives

    a_delta=s_*-delta/[2A lambda J1(j)]+O(delta^2).          (6)

At that point phi=delta and V is bounded away from zero. On the
transition 0<phi<delta, the source is bounded by lambda^2 delta and
the ODE preserves a strictly negative radial derivative. It reaches
phi=0 at a radius c_delta=a_delta+O(delta), with
c_delta=s_*+O(delta), V(c_delta)>0. These statements follow by integrating
(s phi')'=-s f_delta(phi): over a distance O(delta), its bounded source
changes s phi' by O(delta^2), while the term phi'/s is uniformly bounded.
Continuation past the crossing is the actual logarithmic exterior

    phi(s)=-M_delta log(s/c_delta),
    V(s)=M_delta/s, M_delta=c_delta V(c_delta)>0.             (7)

The flat source makes the crossing C-infinity. The entire radial
solution obeys -Delta phi=f_delta(phi), with no imposed radial wall.
Its positive source is supported on the disk s<c_delta, not an annular
quiet-cavity construction. Set Gamma_delta=2pi M_delta; it tends to
2pi A j J1(j)>0. Select this actual circulation before the ring border.

The mirrored choice sigma changes only the swirl W and the eventual
local curl sign. The meridional seed, circulation and all elliptic
border coefficients are exactly the same.

## 3. Full radial mass border, including all angular kernels

Use the fixed-delta radial potential Q(s)=f_delta'(phi(s)). On the
inner disk Q=lambda^2. On the O(delta)-thin transition it is uniformly
bounded by (3), then it is zero. Hence Q converges in radial L1 on
bounded intervals to lambda^2 1_(s<s_*).

Normalize the regular radial homogeneous solution by Y(0)=1,Y'(0)=0:

    (s Y')'=-s Q Y.                                        (8)

It is exactly J0(lambda s) until a_delta. Integrating (8) across the
thin transition and using the ODE's uniform bounds proves that its
exterior logarithmic coefficient is

    B_delta=c_delta Y'(c_delta)
               =-j J1(j)+O(delta) !=0                     (9)

for a fixed sufficiently small delta. This is an explicit strict
mass-border margin. The radial zero-source-mass condition is precisely
zero logarithmic coefficient; it therefore removes the sole regular
radial homogeneous direction. No near-zero numerical eigenvalue is
being used to decide the border.

For angular order m>=1, V=-phi'>0 is the regular m=1 solution and
has the actual exterior V=M_delta/s. Its ground-state identity gives

    integral s[|g'|^2+(m^2/s^2-Q)|g|^2]
      =integral s V^2 |(g/V)'|^2
                        +(m^2-1) integral |g|^2/s.         (10)

Regularity and exterior decaying matching remove the boundary term.
Thus only the two translations occur at m=1, and there is no m>=2
decaying kernel. In the z-even subspace only the x translation remains.
The exact translation/cokernel pairing is

    integral x Q partial_x phi
           =integral x partial_x f_delta(phi)=-Gamma_delta. (11)

This is also the nonzero center-row detection of the translation.
On a disk D_b containing the support with phi(b)<0, let K0 be the
full-plane logarithmic Green map. The complete derivative

    (Y,a,u) -> (Y-K0 QY-a+u x,
                         integral QY, integral x QY)       (12)

is Fredholm of index zero and has zero kernel by (9)-(11). The
constant is fixed by the remaining radial Green equation. Hence
(12) is invertible. This establishes precisely the full-plane mass/
center border needed by 0195 for this new affine-current seed.

## 4. Exact compatibility with the low-carrier profile interface

On the WHOLE positive disk and taper, put Omega=V/s and Z=f_delta(phi).
Then Omega>0, Z>=0, and

    kappa_ep^2=2Omega Z>=0,
    W_sigma W_sigma'=G_sigma G_sigma' phi'
                       =-Z V=-s kappa_ep^2/2.              (13)

W is smooth and flat zero outside the positive disk; V retains its
physical exterior circulation. Thus the changed seed satisfies the
structural profile identities requested by 0206's reflected one-sided
coercivity route. That route's actual pole, tag and continuum conclusions
remain its own calculation. Constant-plateau coefficients have not been
copied to the new Bessel profile.
