# Positive variance repair and a different stationary-shape candidate

## 1. A third actual annulus removes the first amplitude drift

This implements the failure-derived moment-repair-contract after root
acknowledgement. In the simple-node variables y=(x-z)/h choose supports
(-1,e,1), 0<e<1, and POSITIVE material weights (A,1,tau). Write
M_j=sum w_i y_i^j. The signed response weights are w_i y_i/M1, so

    mu_y=M2/M1,
    Var_response=M3/M1-(M2/M1)².

Direct polynomial elimination gives

    M1 M3-M2²=tau e(1-e)²-A[e(1+e)²+4tau].

Hence the strictly positive weight

    A=tau e(1-e)²/[e(1+e)²+4tau]                     (1)

sets the ACTUAL signed phase variance to zero. Both noncentral annuli
still contain positive physical material. The phase sign is supplied
by the nodal Euler displacement, not by negative density.

Put tau=b e². As e->0 at fixed b>0,

    M0 M2/M1² ->1+2b, mu_y=O(e), M1=e+O(e³).

For this three-ring physical marker the initial band equation is

    E+z+h mu_y-(3z/2)[M0 M2/M1²-1]=0.                (2)

The limit has the positive solution b0=(E+z)/(3z) and b derivative
-3z. Together with the nonzero A derivative in(1), this is an
invertible two-equation Jacobian. Fix e small and then h small. The
IFT supplies positive weights solving variance zero and initial band
edge simultaneously. The second quotient derivative has leading term

    mu_PP=(9z²/(2h))(M0/M1)[M0 M2/M1²-1]+O(1)>0.

It preserves the positive squared-clock curvature of band-edge.md.
The characteristic angle amplitude is now h M1, not h; this changes
the actual mass and moment targets and is retained in their finite
normalization. Smooth narrow annuli replace the rings by continuity;
widths are chosen after the finite positive margins of A,e,tau.
The nonlinear real mode F, full pressure/Kelvin correction and packet
filter perturb this same finite Jacobian, so their actual initial
equations are solved with the full data, not a copied nodal ratio.

The spin/G/reference moment controls of0164 remain available on the
open annular supports. Their condition numbers, the amplitude h M1,
and the two added unmarked equations are fixed before the carrier,
packet width and transfer accuracy. The local nonlinear reference
normalization is included as in the existing moment IFT.

## 2. The next exact cumulant and the actual current

The third CENTRAL moment of the signed response after(1) is

    k3=2 e tau(e-1)²(e+1)(e³+2e²+e+4tau)
                              /(e²+e+2tau)³>0.        (3)

For physical x this is h³ k3. It is not removed by the variance
equation. The exact small-time characteristic-function expansion gives

    log[A(t)/A(0)]
      =-i c P mu t/2-(cP)² Var t²/8
                         +i(cP)³ k3_phys t³/48+...,
    gamma=gamma(0)+(cP)³ k3_phys t²/16+...,
    ell=-(cP)² Var t/4+O(t³).                         (4)

Thus zero variance suppresses the leading ell term but leaves a
nonzero gamma_t/gamma current connection. If desired, replace the
second target by the ACTUAL current derivative at t=0:

    Var=c P k3_phys/[2 gamma(0)]

in the leading model. This is an O(delta) perturbation of(1) for
the fast painted clock, so the same invertible Jacobian supplies the
actual first-time current match as well as the initial band edge.
It is a finite-time-jet repair, not a proof that all subsequent
current coefficients vanish. All higher actual F/G derivatives in
0172 remain present.

There is a useful exact rigidity behind this limitation. In a
neighborhood of one simple node, F(x)=(x-z)g(x), and g has a fixed
nonzero sign. Absorb |g| into the positive material measure dnu.
If both the signed central variance and the signed third central
moment about mu vanish, then

    0=integral (x-mu)³(x-z)dnu
       +(mu-z)integral (x-mu)²(x-z)dnu
     =integral (x-mu)²(x-z)² dnu.                    (5)

Positivity forces support only at x=z and x=mu. The node contributes
no angle amplitude. A smooth positive finite-width marker therefore
cannot set both moments to zero. This is a one-node marker-class
identity, not a universal Euler or physical-angle obstruction.

## 3. Full-current rigidity for a finite annular phase sum

Even allowing gamma and amplitude to vary inversely does not make
the two-ring construction's exact current automatically coordinate-only.
For a finite annular response write c_obs(t)=sum b_j exp(-i w_j t),
where all w_j lie in the same nonzero fast-frequency band. Set
R=|c_obs|, J=Im(conjugate(c_obs)c_obs,t)=gamma R².
The current condition gamma_t/gamma+R_t/R=0 is gamma R=constant,
hence

    J²=L² |c_obs|².                                  (6)

Let w_max-w_min=d>0 and retain the nonzero endpoint amplitudes.
The outer Fourier coefficient of J at d is a nonzero multiple of
(w_max+w_min)b_max conjugate(b_min); the sum is nonzero because the
whole band has one sign. Thus J² has a nonzero coefficient at2d,
whereas |c_obs|² has no frequency beyond d. Equation(6) is impossible.
Degenerate equal-frequency groups are first combined exactly; this
argument does not rule out exact cancellation between identical
frequencies or a single surviving frequency.

The same support argument applies to the leading smooth-annulus
integral when its real amplitude has fixed sign near the two endpoints
of its nonzero compact frequency support. Near the outer difference
frequency the convolution integrand has fixed sign and w+w' stays
away from zero; no endpoint cancellation removes that support. Its
square reaches twice the support width. Near a simple real node with
positive material density these endpoint hypotheses hold. This is an
elementary signed-convolution support calculation, not an assumption
that an approximate mode is an exact Euler eigenmode.

Accordingly finite positive annuli give controlled initial and finite
jet repairs, but not the claimed exact all-window current within this
painted fast-band representation. The full Euler packet still has
its actual remainder; the conclusion above concerns the explicit
response family, and does not close other stationary-observation or
normal-mode routes.

## 4. Different representation: positive stationary radial shape

The already registered candidate B uses the stationary elliptic
material-domain covariance of0155, not fast painted lobes. Its
ground-state radial response is F(x)=exp(-x/2), E=2 and the leading
physical slow clock at a fixed material marker is

    gamma(P,0)/c=-E P^-1/2+P mu(P)/2,
    mu(P)=integral x exp(-x P^1.5/2)dnu
                         /integral exp(-x P^1.5/2)dnu.

Here the response is a positive probability measure. Write its
reference mean, variance and third central moment as m,V,k3. Exact
exponential-family differentiation gives

    mu_P=-3V/4, mu_PP=-3V/8+9k3/16,
    gamma_P/c=E/2+m/2-3V/8,
    gamma_PP/c=-3E/4-15V/16+9k3/32.                  (7)

Choose the POSITIVE reference response with weights32/35 and3/35
at x=5 and x=50/3. It has

    m=6, V=32/3, k3=928/9.

Equation(7) then gives

    gamma=c>0, gamma_P=0, gamma_PP=35c/2,
    partial_P² gamma²=35c²>0.                        (8)

These are actual physical moment values, not a supplied frequency.
The physical radial material weights are obtained by undoing the
positive factor exp(-x/2), followed by a common scaling to keep
0<=chi<=1. Positive smooth narrow annuli and a small adjustment of
one moment preserve the simple initial band equation and both strict
signs. In particular eta=A_F²/(C_tag N) is strictly positive, although
it obeys the0155 upper bound1/3 and is not asserted equal to1.

The all-time observation is still its actual characteristic function.
For the sign convention in(7), gamma_t(0)=ell(0)=0 but

    partial_t(gamma_t/gamma+ell)(0)=-140c²/9.           (9)

Its dephasing is order one on the NATURAL slow optical scale t~1/c.
It cannot be declared negligible there just because c/Omega is small.
Thus(8) is a positive stationary-shape INITIAL band-edge construction;
an exact full-window mode/current would need a genuinely stationary
observation of an actual Euler normal mode or an additional measured
boundary/current reconstruction. The source's smaller-overlap option
has not been discarded by the old eta=1/2 mismatch.

## 5. Why positive overlap is a useful next route, not an eta=1 gate

Retain the actual overlap eta in a prospective stationary mode, instead
of identifying its spin with canonical momentum by name. If its
constant action coefficients are j,kappa>0 and its measured integrated
spin is eta j q, the leading physical displacement/angle map is

    U=X-eta j h q/(2rho), Phi=q+hX/2,

where h is the real curl-helicity multiplier. Pulling back the SAME
diagonal phase action gives, at order h,

    M_U,Phi=j(eta-1)h/2,
    K_U,Phi=-kappa h/2.

At the optical gap omega²=kappa/j the actual translational forcing
row is -eta kappa h/2. Equivalently with M_U,Phi=-b h and
K_U,Phi=-g h, the physical invariant is

    g-kappa b/j=eta kappa/2>0.                        (10)

So a positive nonunit overlap does not erase physical coupling.
The remaining current and gradient inertia must remain in the actual
map. Equation(10) is an exact interface implication, NOT a claim that
the nonstationary marker in(8) already has that constant full-time row.
It gives the next normal-mode/domain construction a concrete positive
target without adding rigid mass or silently imposing eta=1.

## Route ledger and continuation

The fixed positive two-annulus and repaired three-annulus INITIAL
band-edge constructions are established. The inference of exact
all-window autonomy from finitely matched painted-marker moments is
refuted for the explicit nondegenerate fast-band response, with its
interference/current mechanism named. The stationary-domain ground
mode supplies a different positive initial band edge and positive
nonunit overlap; its full optical-window normal-mode/current realization
remains that route's missing construction, not an obligation no-go.
An exact stationary observation, with its boundary flux and actual
normal-mode current, is the next implementation target. The full
parent objective and other registered stationary-law routes stay active.
