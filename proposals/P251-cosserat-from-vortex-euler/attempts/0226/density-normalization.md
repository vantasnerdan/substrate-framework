# Nonzero density scaling and a separate positive-action normalization

This is the activated candidate C's first exact construction. It proves
a finite-moment and full quadratic-action normalization in the actual
straight literal-curl supplier of 0222, and derives its ring density
scales. It does not replace the curved whole-space response estimate by
an unproved higher-order estimate.

## 1. The actual positive forms

In the straight inner region,

    Omega(s)=A lambda J1(lambda s)/s,
    Z'(s)=-lambda^2 s Omega(s).

For the m=-1 actual Kelvin lift of 0222, and its real phase pair, the
full inherited phase and energy MATRIX are

    beta=rho pi L lambda^2 integral s Omega(s) A0(s)^2 ds,
    H=rho pi L lambda^2 integral s Omega(s)^2 A0(s)^2 ds.
                                                               (1)

They are positive on every nonzero compact radial preparation in the
positive inner region. Thus H=nu beta is a quadratic mean condition,
not a signed linear-response condition.

Fix a smooth positive tag chi(s/a), with chi' nonzero on an open annulus
strictly inside the literal-curl core. Set A0(s)=a F(s/a), normalize the
tag response functional L(F)=integral b(x)F(x)dx to one, and include the
compact inverse-return condition integral x^2 F dx=0. The latter may be
completed by a separate off-tag bump as in 0222. Every extra condition
below is on an actual preparation, not on tag probability.

Let epsilon=(lambda a)^2 and

    Omega_0=A lambda^2/2,
    y_epsilon(x)=[Omega_0-Omega(ax)]/(Omega_0 epsilon).

On a fixed bounded annulus y_epsilon=x^2/8+O(epsilon), with all fixed
derivatives, and y_epsilon is strictly increasing for sufficiently small
epsilon. Choose c strictly inside its annular range with a fixed margin,
and nu=Omega_0(1-c epsilon)>0. For a fixed integer q impose

    L(F)=1,
    L[(y_epsilon-c)^j F]=0,  1<=j<=q.                  (2)

The compact return remains an additional homogeneous linear row. Further
finite rows can be retained when their compatibility with L(F)=1 has
been independently established; the construction does not presume that
an arbitrary additional row is independent of L.

## 2. Exact nonlinear action repair, with positive probability

The new action condition is

    Q_epsilon(F)=integral x Omega(ax)
                         [Omega(ax)-nu]F(x)^2 dx=0.      (3)

This condition can be solved simultaneously with (2) and the compact
return. Here is an explicit function-space construction, rather than
an assumed sign cancellation of energies.

First select q+1 small disjoint bumps in the region where b is nonzero,
with distinct y-values. Their normalized moment matrix tends, as the
bump widths shrink, to a Vandermonde matrix. It is therefore invertible
for fixed sufficiently small widths. Solve for a signed F0 satisfying
(2), and add an off-tag return to cancel integral x^2 F0. Keep two
further disjoint open intervals, one on each side of y=c, unused.

In either unused interval the space of compact smooth functions killed
by all the finitely many homogeneous linear rows is infinite-dimensional.
Choose nonzero Z+ and Z- there. They have disjoint support from F0 and
each other. The quadratic form Q is strictly positive on Z+ (y<c) and
strictly negative on Z- (y>c). If Q(F0) is nonzero, choose the opposite
sign Z and set

    F=F0+t Z,   t=sqrt(-Q(F0)/Q(Z)).                    (4)

All linear rows stay fixed, while Q(F)=0 exactly. If Q(F0)=0 use F0.
No signed ensemble probability has entered: F is a signed coherent
initial displacement in one positively weighted realization. Both
beta and H remain positive integrals of F^2. Their ratio is exactly nu.

Divide Q by Omega_0^2 epsilon before taking epsilon to zero. The
resulting form converges smoothly to integral x(c-x^2/8)F^2. Choose
the supports with fixed margins from y=c. The moment matrices and the
nonzero denominator in (4), in these normalized units, have finite
nonzero limits. Consequently the construction costs only fixed
q-dependent profile constants, not extra inverse powers of epsilon.
One may choose the initial F0 so Q(F0) is separated from zero, yielding
a nondegenerate root and a smooth local continuation. A zero Q(F0)
already solves the exact condition; no division by it is used.

For any fixed finite time and derivative order the physical response is

    L[F exp(-i Omega(ax)t)]
       =exp(-i nu t)[1+O_q((Omega_0 epsilon T)^(q+1))].    (5)

Taylor's theorem and (2), with the fixed L1 norm of bF, prove this
estimate and its differentiated versions. Initial time derivatives
through q match exactly. The actual action (1) simultaneously has
H=nu beta, so a signed linear flattening has not been mistaken for
quadratic-action flattening. The unmeasured radial complement remains
an actual continuum of Euler histories; (5) is an observed finite-time
statement, not an invariant single eigenmode.

The leading ring tag-spin row is the same first frequency moment of
the response. Its initial ratio to angle rate is positive in this
construction. The positive tag fraction may then be selected by the
actual 0222 phase/spin matching formula. That fraction scales the
mechanical tag current, not the phase or the velocity field. Curved
spin, G and shape corrections are not removed by (2)-(4).

## 3. Actual density scales

Lengths below are expressed in one fixed reference length ell_0; lambda
denotes its dimensionless curl parameter. Keep Omega_0 fixed, take
A=2Omega_0/lambda^2, a=lambda^(-1/2), and set R=bar_R/lambda,
bar_R tending to infinity. The preparation F and its signed returns
above have fixed scaled support; thus a is much smaller than the full
core radius of order lambda^-1. The derived 0217/0222 integrals give

    Ic proportional to a^2,  Ia proportional to a^4,
    B proportional to a^2,  c_theta proportional to R^-1,
    beta proportional to rho R lambda^2 Omega_0 a^4,
    M_phase=beta/(nu c_theta^2)
                         proportional to rho R^3 lambda^2 a^4.

For the proposed number density P^-3, P=d R with d fixed and sufficiently
large to contain the vorticity ring, this yields the bookkeeping identity

    M_phase/P^3 proportional to rho lambda^2 a^4
                              =positive finite constant.        (6)

Equation (6) establishes a nonvanishing density scaling available to a
future actual stationary array or ensemble. It does not construct that
array by superposing rings: the ring velocity and pressure extend beyond
its compact vorticity, and the nonlinear cross interactions remain a
separate stationary-field construction.

The positive tag fraction is proportional to lambda^2 a^2=lambda,
so it is admissible for small lambda. Its tagged mass per volume tends
to zero like the inertia density divided by R^2: the finite inertia
comes from the real lever arm R. This distinction is essential. It
does not turn a dilute tagged mass into all-ambient translational mass
rho. Infinitesimal material amplitudes are taken before this geometry
limit; a fixed finite angle would eventually leave the thin core.

The background velocity scale is Omega_0/lambda, the tube cross-section
is of order lambda^-2, and its length is of order R. Its kinetic energy
per volume scales, up to fixed profile and exterior factors, as
rho Omega_0^2/(lambda^2 bar_R^2). This is another independently retained
quantity, not the optical inertia density (6).

## 4. What the density repair does and does not solve

The old fixed-a dilution is repaired at the level of the actual derived
phase and tag scales. The nonlinear condition (3), rather than a signed
ensemble energy average, repairs the reference action normalization.
These are positive established constructions.

The existing curved estimate is only a fixed-parameter estimate of order
log(bar_R)/bar_R, with its actual parameter-dependent constants. Merely
declaring K R small cannot make that error o(K^2). For example with
bar_R=lambda^-b and K=o(lambda^(b+1)), the displayed geometric error
scales as lambda^b log(1/lambda), larger than K^2. No higher-order
curvature theorem follows from (6).

The next analytic construction is a full-pressure compact multipole
recursion or another actual curved response improvement. At each finite
order it must retain pressure sources, their exterior moments, physical
tag rows, and the independent quadratic normalization. The elementary
fact that finite radial moment constraints leave an infinite-dimensional
preparation space supplies possible controls, not a proof that all of
that recursion's forcing rows are compatible. A finite-radius Richardson
mixture likewise cancels only its explicitly solved linear rows. Its
positive quadratic forms would need their own actual solution, such as
(3), on the curved geometry. The shared-field acoustic/ambient-current
construction remains active.
