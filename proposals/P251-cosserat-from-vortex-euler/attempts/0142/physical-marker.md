# One transported marker: physical curvature and momentum through time

## 1. Physical observation, fixed before carrier differentiation

Use the actual Euler columns of lundquist-packet.md, m=2,n=2. Choose
a sufficiently small reference delta_* and p_*=lambda/delta_*² first.
All marker functions below are then fixed in physical space as p varies
in a fixed relative neighborhood of p_*. A subscript * on a length
denotes a construction choice, not a variable in partial_p.

Let B_*=(C log(1/delta_*))^(1/2), b_*=B_* ell_*, and choose a smooth
nonnegative chi(r), equal to r²/(4b_*²) for r<b_*, zero for r>2b_*,
and at most1. Choose an even nonnegative axial marker chi_z, with
Z0=integral chi_z and Zc(k)=integral chi_z cos(kz) nonzero on the
carrier interval. Its width can be a fixed small multiple of1/p_*.
The Fourier-fiber action is assigned a fixed axial action measure Lz,
also a fixed multiple of1/p_*. Lz is not changed to2pi/p when taking
carrier derivatives. This is a fiber action normalization, not a claim
that restriction of an infinite column to a finite axial window is
itself an invariant Euler subsystem.

The initial material fraction is

    w0=mu chi(r)chi_z(z)[1+eps b_tag(r)cos(2theta)],
    0<mu<=1/2, |eps b_tag|<1.                                  (1)

It is transported by the actual Euler flow, including displacement
of its boundary/level sets. Signed b_tag specifies opposite local
quadrupole orientations, not negative mass. Its initial central
quadrupole is proportional to the nonzero real number

    Q=integral chi b_tag r³ dr.                                (2)

C2 symmetry makes the planar centroid zero. The observed angle is
half the argument of the central quadrupole, relative to that of the
same marker transported by the unperturbed flow. This definition fixes
its clock and orientation; no Floquet winding is available for sign
selection. The unperturbed quadrupole is exactly proportional to

    integral chi b_tag r³ exp(2i O(r)t)dr,                      (3)

not Q exp(2iOmega t) when O varies.

The complete material angular momentum variation is obtained by
varying both position and velocity under the material map:

    delta S_z=rho integral w0 r[v_theta+Z(r)xi_r] dV,
    Z=2O+rO',
    r[v_theta+Z xi_r]=m P/s                                   (4)

for an exact Fourier mode, by the azimuthal Euler and Lin equations.
For a quasimode the same identity includes its displayed residual.
The rO' term is necessary on the actual Lundquist background. Pressure
and the entire material-position term are kept. The observable is
tagged axial spin, with surrounding-fluid and marker-boundary exchange;
it is not the total angular momentum of the whole infinite column.

## 2. Why this particular excited packet has positive physical curvature

Put x=g r², L(x)=L_2^1(x)=x²/2-3x+3, and

    cD=sqrt(2)Omega delta, sigma=-2Omega+6cD+O(delta²).

At leading order the quadrupole variation uses the radial weight
chi r² f_2 dr. On the Gaussian region it is proportional to
x² exp(-x/2)L(x)dx. The transported axial phase is
exp(-i cD x t/2). Its normalized integral is evaluated exactly:

    a(t,p)=(w²-6w+8)/(3w^5), w=1+i cD t,
    a(0,p)=1, integral x³ e^-x/2 L / integral x² e^-x/2 L=38/3.
                                                                    (5)

The only numerator zeros are w=2,4, away from Re w=1. Thus the leading
angle row never vanishes at any real t. On a fixed optical window the
actual row remains nonzero for small delta. Its phase derivative is

    gamma=sigma+partial_t arg(a)+actual profile/base corrections
         =-2Omega-cD/3+O_T(delta² polylog(1/delta)),
    connection=partial_t log|a|+corrections
              =O_T(Omega delta² polylog(1/delta)).              (6)

The reference-shape correction in (3) is included here: O-Omega is
O(Omega delta³ B_*²), while Q has one small cancellation of order
delta, so its normalized contribution is at most O(delta² polylog).
The exact profile corrections start at relative order delta. Their
effect on the phase shift, which is already order delta, starts at
delta². The Gaussian cutoff error and its first two scaled carrier
derivatives can be made O(delta^N) by increasing C before freezing
the marker. There is no radial boundary in the field equations.

For a function of delta alone,

    p² partial_p²=delta² partial_delta²/4+3delta partial_delta/4.

Consequently, with every marker and Lz held fixed,

    gamma²=4Omega²+(4sqrt(2)/3)Omega² delta+O_T(delta² polylog),
    p² partial_p² gamma²=sqrt(2)Omega² delta+O_T(delta² polylog)>0.
                                                                    (7)

This is curvature of the physical transported-marker clock. The
intrinsic sigma² curvature has the opposite sign and remains recorded
in analytic-route.md. Equation(7) does not identify the laboratory
carrier kU+2Omega with a macroscopic optical gap.

## 3. Full physical-angle action, including time connections

The two actual Euler columns preserve their exact initial KKS beta>0.
Let the actual angle row in solution coordinates have length |c(t,p)|
and continuously chosen phase psi(t,p), with gamma=psi_t. The sign of
c is fixed at t=0. For sufficiently small delta, gamma<0. The exact
time-dependent one-form pullback of0115/0140, followed by eliminating
the complementary phase coordinate, is

    L_theta=M/2[(theta_dot-connection theta)²-gamma² theta²],
    M=-beta/(gamma c²)>0, connection=c_t/c,
    Pi=M(theta_dot-connection theta).                          (8)

For the prepared sine phase theta=c sin psi, Pi=-beta cos psi/c.
This is an exact action on the two-dimensional space of actual linear
Euler histories, with its time-dependent observation map. It is not
obtained by projecting Euler onto a guessed invariant finite span.
The actual quasimode comparison controls its coefficient errors.

The scalar equation retains the connection and mass derivatives. Its
coefficient of theta is

    gamma²+connection²-connection_t+(gamma_t/gamma)connection,   (9)

and its damping/metric-connection coefficient is
-gamma_t/gamma-2connection. The extra terms in (9), including their
first two scaled carrier derivatives, are O_T(delta² polylog). Thus
the positive leading curvature in (7) survives in this physical
equation as well. No autonomous action or unqualified dispersion law
is being inferred by deleting these terms.

## 4. Seven fixed marker controls: the time/carrier mismatch is repaired

At leading order the spin-pressure polynomial is

    P(x)=L(x)-2L'(x)=x²/2-5x+9.                                (10)

After factoring the nonzero dimensional prefactors, the three
carrier-jet spin rows and three first slow-time rows are

    e^-x/2 P, D(e^-x/2 P), D²(e^-x/2 P),
    x e^-x/2 P, D(x e^-x/2 P), D²(x e^-x/2 P),
    D=-1/2+(3/2)x partial_x.                                 (11)

Replacing p²partial_p² by (p partial_p)²-p partial_p is an invertible
triangular row change. These six rows, after extracting e^-x/2, span
all polynomials of degree at most5: their coefficient determinant is
19683/8192. Together with the independent reference row1 this gives
seven independent analytic functions on every open radial annulus.
More explicitly the equivalent family
{1,e^-x/2,x e^-x/2,...,x^5 e^-x/2} is an extended Chebyshev system:
six derivatives of e^(x/2) times a purported combination leave a
nonzero constant times e^(x/2), unless the coefficient of1 is zero;
Rolle then proves the required zero-count bound. Seven distinct
radial points and sufficiently narrow, fixed smooth bumps therefore
give a nonsingular moment matrix. Their supports and widths are
frozen after p_*, not transported to different radii as p changes.

Use these bumps for b_tag. The exact target is the actual canonical
momentum from (8), not an assigned inertia. At t=0 solve simultaneously

    Q=Q_selected !=0,
    partial_p^j[S(0,p)-Pi(0,p)]_(p*)=0,
    partial_p^j[partial_t S(0,p)-partial_t Pi(0,p)]_(p*)=0,
    j=0,1,2,                                                  (12)

where the first line of momentum conditions is taken in the sine-angle
phase and the second in its cosine-angle partner. The unused phase
rows vanish by reflection parity. For the second group subtract the
already matched common sigma contribution and divide by cD. Its
leading rows are precisely the second line of (11).

Equation(12) is linear in the seven marker coefficients. At t=0 the
angle denominator is Q, so Pi=-beta/c0 is linear in Q. Its first time
derivative contains the reference-shape phase rate (3); multiplying
it by Q leaves the linear O-weighted reference row, not a nonlinear
constraint. The exact pressure, profile, KKS and reference-shape
corrections perturb the dimensionless moment matrix by O(delta
polylog). Its inverse therefore persists for sufficiently small delta.

Choose Q_selected of order delta in the normalized reference moment,
so the required bump coefficients are bounded despite the spin row's
leading order-delta cancellation. Scaling Q and all coefficients by
one common small number enforces |eps b_tag|<1 without changing any
matching equation. These are derived geometric marker weights, not
fitted frequencies or fitted elastic constants.

For a coherently prepared standing pair, one can instead prescribe the
same fixed positive weight eta in S-eta Pi in both groups of (12).
The coefficient matrix is unchanged up to the same reference-row
subtractions; eta=1/2 is the standing-pair normalization of0138/0140.
It is chosen before the tag is solved, not inferred by averaging two
independent reaction momenta. The physical label mass stays nonnegative.

The n=0 failure recorded in analytic-route.md is genuine: its
carrier and slow-time rows are dependent with incompatible targets
for positive curvature. The n=1 six-row matrix still has rank5.
The n=2 determinant in (11) is the materially different repair.

## 5. Natural-scale error, not absolute quasimode accuracy

All estimates are for |Omega t|<=T0, j<=2 scaled carrier derivatives,
and the two real prepared phases. Gaussian moment ratios and the fixed
seven-bump inverse are bounded after the explicit delta normalization.
Taylor expansion of the *actual transported exponential*, through
its first slow-time moment, has remainder bounded by
C_T delta² times the second absolute x moment. The same bound holds
after two p derivatives, because differentiating x=g r² only inserts
polynomials in x. Equation(12) cancels the constant and first slow-time
terms of both physical spin and Pi. Hence

    max_(j<=2) |p^j partial_p^j(S-Pi)| / Pi_scale
      <= C_T delta² polylog(1/delta).                          (13)

The bound is a row norm over both phases; it does not divide by a
particular trajectory's vanishing instantaneous Pi. Pi_scale is its
nonzero initial row norm. Time derivatives needed in (8)--(9) satisfy
the corresponding finite-order estimate.

For transparency the actual Euler replacement is made at a much
higher order than (13). After scaling r=ell R and the fixed axial
length by1/p_*, all marker norms and the finite moment inverse have
only polynomial losses in delta and logarithmic factors. Two scaled
carrier derivatives of the exact Euler propagator cost at most p²,
or delta^-4. Two time derivatives cost at most another p², or
delta^-4. Dividing by the small spin/reference moment costs at most
delta^-2; weighted-observation normalizations can be covered by two
further inverse powers. In particular, for the dimensional profile
f=r times a dimensionless Gaussian polynomial,
||xi||_2=O(ell² sqrt(Lz)), ||r w0||_2=O(b_*² sqrt(Z0)), and the
matched spin row is of order rho Omega delta ell^4 Zc/B_*². Their
Cauchy--Schwarz ratio costs only delta^-1 B_*^4, since
sqrt(Lz Z0)/Zc is bounded by construction. The normalized second
group of matching rows costs one additional delta^-1. Thus there
is no inverse power hidden in the shrinking physical marker volume.
The two spare inverse powers cover reference-row divisions and their
finite products. The deliberately conservative bound
for these normalized rows is

    C_(N,T) delta^(N-12) polylog(1/delta).                      (14)

The arbitrary finite-order construction in lundquist-packet.md can
be taken to N>=16 before selecting delta. This makes (14) smaller
than the O(delta²) physical remainder. Equivalently, if one uses a
different fixed family of marker norms, first compute its finite
algebraic loss exponent d and choose N>d+2. Nothing requires a
uniform-in-N estimate or an assumed exact pole. The original Kelvin
defect is removed at initialization, not hidden in (14).

The positive curvature signal in (7) is order Omega² delta/p²;
its full coefficient error is O_T(Omega² delta² polylog/p²).
Their ratio tends to zero. Choose delta so the latter is below, for
example, one quarter of the leading signal. Then choose the actual
nonlinear perturbation amplitude last so its finite-time quadratic
remainder and first two preparation derivatives lie below the same
margin. This is a finite-time actual-Euler physical observation
theorem, not numerical evidence below an unmeasured floor.

## 6. Exact scope earned

The construction supplies a smooth constant-lambda Euler background,
pressure-resolved optical packets with exact Kelvin initial data,
positive KKS and scalar material-angle mass, a fixed nonnegative
transported marker, positive physical second-carrier curvature, and
spin/action agreement through two carrier jets with controlled
order-delta² error on a finite optical window. Every physical time
and parameter connection is retained. Its axial action is explicitly
Fourier-fiber normalized; full ambient spin exchange is not deleted.

It does not prove an isolated Lundquist eigenvalue, a nonlinear
periodic rotor, an exact autonomous continuum, or a knotted/global
periodic EPS transfer. The parent's0145 and acoustic constructions
remain distinct obligations. No accepted claim is modified here.
