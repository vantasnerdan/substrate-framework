# A fixed physical two-annulus marker with positive curvature and small current connection

## 1. The fixed-radius derivative refutes the single-annulus shortcut

Keep the ACTUAL n=8 Lundquist Euler packet of0147. Only its nonnegative
radial material fraction is redesigned. In its leading circular core let
x=g_* r², P=p/p_*, cD_*=sqrt(2)Omega delta_*, and

    F(x)=exp(-x/2)L_n^1(x),
    sigma(P)=-2Omega+(2n+2)cD_* P^(-1/2)+O(delta_*²).

For a fixed physical radial marker, its exact leading angle numerator,
after the intrinsic phase, is proportional to

    A(P,t)=integral dnu(x) F(x P^(3/2))
                              exp(-i cD_* P x t/2).       (1)

Here dnu is a nonnegative FIXED physical radial measure (including the
known r³ Jacobian); the mode F may change sign. Common real carrier
normalizations cancel from its phase but not from its mass. Equation(1)
comes from the transported material quadrupole, not a renamed eigenphase.

A single narrow annulus at x=x0 makes its phase
sigma(P)-cD_* P x0/2. The second term is linear in P. Therefore its
leading squared-clock curvature is NEGATIVE:

    partial_P² gamma² = -3Omega(2n+2)cD_*+O(delta_*²).

Choosing x0 so that gamma(P=1) itself has a favorable shift does not
repair this sign. Holding the scaled x fixed during differentiation
would move the physical tag with the carrier and give the wrong answer.

## 2. Two positive material annuli exploit an actual mode node

Choose any simple positive zero z of L_8^1. Fix small h>0 before the
high-carrier hierarchy. Take two positive smooth annuli centered at
x_+=z+h and x_-=z-h, with their radial measures normalized to weights
2 and1. Their widths are O(h²); all supports and weights remain fixed
during carrier differentiation. The delta-ring calculation below is the
leading smooth-annulus asymptotic, not the physical material definition.

At P=1, F(z+h)=F'(z)h+O(h²), F(z-h)=-F'(z)h+O(h²).
The total angle numerator is F'(z)h+O(h²), nonzero for small h.
Its normalized response weights approach 2,-1. These signs belong to
the actual nodal Euler displacement. Both material densities are positive.

Put

    mu(P)=[2 x_+ F(x_+ P^(3/2))+x_- F(x_- P^(3/2))]
                    /[2 F(x_+ P^(3/2))+F(x_- P^(3/2))].   (2)

The exact quotient derivatives, using F(z)=0,F'(z)!=0, give

    mu(1)=z+3h+O(h²),
    mu'(1)=-12z+O(h),
    mu''(1)=108 z²/h+O(1).                               (3)

For transparency, write d=z(P^(3/2)-1). The leading nodal quotient is
z+h(3h+d)/(h+3d); differentiating it gives the same singular coefficient
in(3). The full smooth F changes only the displayed bounded remainder.
Symmetric smooth bumps of width O(h²) preserve these asymptotics by
Taylor's theorem, including two P derivatives. Choose a carrier interval
|P-1|<c h/z; its numerator stays away from zero by a fixed multiple of h.

The actual initial physical clock is

    gamma(P,0)=sigma(P)-cD_* P mu(P)/2+O(delta_*² C_h).

Consequently

    partial_P² gamma(P,0)²
        =216 Omega z² cD_*/h+O(Omega cD_*)
                                      +O(delta_*² Omega² C_h)>0. (4)

Select h first and delta_* sufficiently small afterwards. No eigenvalue
or frequency was supplied; (4) differentiates the actual material phase.
The finite packet's Gaussian carrier observation gives the same positive
a²=1/(1+c²)² prefactor from0147 and controlled narrow-band corrections.
Choose p_*L large relative to the now fixed annular derivative constants.

## 3. The physical clock variance, not its curvature, controls the connection

At P=1 the signed second central moment of the two-ring response is

    Var(x)=(-8 h²)+O(h³).                                (5)

A negative signed response variance is allowed: it describes the
interference of a nodal displacement under a positive material measure,
not the variance of a positive probability distribution. For the smooth
annuli the total variation of the normalized response remains bounded.
The two phase frequencies differ by O(cD_* h). Direct Taylor expansion
of(1), or its two exponentials before the smooth limit, gives on
|Omega t|<=T0

    ell=(partial_t |A|)/|A|=O_T(Omega delta_*² h²),
    gamma_t/gamma=O_T(Omega delta_*³ h³),                 (6)

at leading order. At higher pressure order use the actual simple node
of the real radial angle amplitude, shifted by O(delta), and retune the
two positive annular weights by their smooth O(delta/h) corrections.
The fixed physical radial supports are then frozen at p_*. The true
material frequency width remains O(Omega delta h), with relative
O(delta C_h) corrections, so retained physical pressure terms do not
create an unexplained order-delta amplitude drift. The same fixed-radius
derivative of the actual node preserves(3)-(4) for sufficiently small
delta after h.

Actual quasimode residuals, finite-band, nonlinear reference and geometry
errors are retained as an additional differentiable packet-own error
epsilon_dyn. Arbitrary fixed finite orders of the imported construction
control APPROXIMATION residuals, not physical pressure coefficients.
Those coefficients were kept in the preceding node/frequency-width
argument. All accuracies are selected after h and its finite conditioning
constants. There is no uniform h->0 assertion.

For the actual scalar action j=-beta/(gamma |c|²), Pi=j(theta_t-ell theta),
the exact identity is

    j_t/j=-gamma_t/gamma-2ell,
    (-j ell)-j_t=j(gamma_t/gamma+ell).                    (7)

Preserve the ACTUAL initial G0=j(0)theta(0) with the added0158 in-tag
control, and match the actual spin S=Pi+e at the declared finite order.
Then the integrated physical angular current obeys

    H=G0+integral_0^t S
     =j(t)theta(t)+integral_0^t j(gamma_t/gamma+ell)theta
                                                    +integral_0^t e. (8)

Thus its non-coordinate part is bounded, in the physical optical phase
norm, by C_T j[delta_*² h²+epsilon_dyn] times the initial angle/rate scale.
The corresponding Poisson row has bound C_T[delta_*² h²+epsilon_dyn]/Omega.
This is a controlled ACTUAL current remainder, not an identity H=j theta.

## 4. Existing spin and initial G controls survive as actual material controls

The selected chi is a nonnegative sum of two smooth annuli, positive on
open subannuli. The twenty radial reference/pressure/carrier functionals
of0147 and the additional G row of0158 remain independent there: multiplying
their analytic functions by a fixed positive radial weight does not create
a relation on an open interval. Smooth bump controls supported there give
a finite invertible moment matrix. Its conditioning may grow as h shrinks;
all those constants are fixed BEFORE delta, packet length and transfer
accuracy are selected.

The angle numerator in(1) has size h, so the normalization is not copied
from the old r² marker. With the actual density factor mu_tag included,
the reference quadrupole and the target spin/G equations are scaled by
their actual angle numerator. In the leading normalized system a choice
Q_* of order mu_tag delta_* h gives bounded forcing; the reference
quadrupole remains strictly nonzero. The same enlarged finite moment
IFT solves the old forty matching/reference equations and G0=j0 theta0,
using their newly computed targets. The actual nonlinear denominator and
its Jacobian are retained, as in0147; its reference dephasing is controlled
by the preserved three exact background moments, not a fixed clock.

Reduce the common positive tag fraction after the finite annular shape
is fixed so mu_tag ||chi||_infinity(1+|eps b|)<1. A common marker scaling
maintains a strict margin. This affects the actual spin/action matching
normalization and is included there. It neither signs physical density
nor requires j to stay nonzero in a zero-width or infinite-volume limit.
Every selected finite h,delta and cell has a nonzero tag and positive j.

## 5. Natural scales and the continuum meaning

Let epsilon_c=delta_*² h²+epsilon_dyn. The first-gradient hybrid shift is
U=X-curl H/(2rho), with the symmetric-shape current still explicit. Its
additional noncommuting position row from(8) is bounded by

    |{U,theta}-{X,theta}| <= C_T |K| epsilon_c/(rho Omega).

Relative to the natural first-gradient optical/translation coupling this
is an epsilon_c error, not an all-k configuration identity. At the optical
gap, any second-order dispersion comparison also keeps the actual ratio

    (j/rho) Omega² epsilon_c / B,
    B=partial_p²(gamma²)/2 ~ C Omega² delta_*/(p_*² h)>0.  (9)

It is this ratio, together with the exact matrix pullback, that must be
small when optical curvature is the consumer. For any selected finite
packet the density j=J_packet/Vcell can be chosen small but positive by a
finite cell/population choice; packet-own optical curvature and spin errors
are divided by that SAME volume and are not weakened by doing so. No
uniform Bloch radius under that cell choice is inferred: macro K and the
fixed finite time window must then lie within the actual chosen-cell bounds.

This construction controls the angular current part. The remaining
symmetric-shape current and the actual fixed-cell joint matrix are retained
in the main response proof. The ordinary continuum is licensed only to
its declared spatial order and these explicit relative errors, not as an
exact unrestricted Euler invariant manifold.
