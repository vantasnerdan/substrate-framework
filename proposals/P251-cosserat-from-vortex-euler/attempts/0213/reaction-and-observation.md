# First fast reaction and actual Euclidean observation

## 1. A computed reaction at the bending action scale

Let epsilon=1/R, k=n/R with fixed n>=2, and take a positive spectral
contour with omega comparable to epsilon² log R and
c=omega/k in the positive sector. Fix the harmonic amplitudes and
profile; the actual residual from `curved-residual.md` has size

    ||eta_b||_{H^j(core)} <= C epsilon² log R.

Project this computed source onto its straight axisymmetric channel.
The exact pressure/coadjoint lift (7)-(10) of `weighted-graph.md`
defines its response Xi_0. The full ring volume measure contributes
one factor R. The absolute inherited KKS reaction therefore satisfies

    |rho integral omega dot (eta_b cross Xi_0)|
       <= C rho R ||eta_b|| ||Xi_0 cross omega||
       <= C rho R k^-1 log^N R ||eta_b||²
       <= C rho R^-2 log^(N+2) R.                         (1)

This estimate uses the actual computed source and the proved response
norm, not a postulated matrix inverse. It includes the particular
pressure terms in the response. The sign convention for the KKS form
changes its overall sign but not the bound. For a fixed profile the
positive bending Hessian from 0206, with THIS source circulation, has
size rho Gamma² log R/R. The ratio of (1) to that positive scale
tends to zero as R^-1 times a finite logarithmic power.

Thus the first response of the formerly dangerous threshold channel
is genuinely lower order in the inherited action. This does not by
itself sum its feedback through the curved operator.

The exact next Schur object is obtained by fixing the full physical
KKS projection onto the bending coordinates and its complementary
coadjoint space. On n=1 the separate exact Euclidean projection is
used, including its translating-frame Jordan row. On n>=2 the
centerline translations are retained as bending coordinates, not
discarded as global symmetries. For blocks of the ACTUAL generator,

    z−A_bb−A_bf (z−A_ff)^−1 A_fb                          (2)

is the full reaction. Formula (1) controls its first straight m=0
term with the actual source, and supplies an appropriate norm for
its continuation. The full (z−A_ff)^−1 is not replaced by that first
term without an estimate.

## 2. The exact remaining curved feedback estimate

The poloidal-zero projection of a first-order cosine/sine curvature
coefficient vanishes. A zero-average first-order core deformation
also has this selection rule. A uniformly invertible nonzero-poloidal
complement would then produce a two-step m=0 correction of order
epsilon²; paired with k^-1 polylog(1/c), this is small at k=n epsilon.

The quantity to compute is this two-step pressure/coadjoint operator
in the weighted source norm used above, including the actual
flux-coordinate pullback. A naive unweighted operator estimate does
not decide it: Xi_theta in (7) can be large where w is flat, although
its induced force and velocity have the proved bounds. Treating that
large gauge component as an independent forcing, or dividing by the
vanishing vorticity before retaining the pressure cancellation, would
lose the meaningful estimate. This is the precise remaining
construction for a full curved pole on the optical-time window.

## 3. An actual finite-time history already follows

Let xi_b(t) be a bounded smooth bending history with time derivatives
of order log R/R², and v_b=V xi_b. Let v_exact solve the complete
linearized Euler equation about the actual global 0211 ring with
v_exact(0)=v_b(0). No cell or exterior wall is introduced. For each
fixed Sobolev order s, the actual residual has

    ||(partial_t−L)v_b||_{H^s(R³)}
       <= C_s R^(1/2) log R/R².                           (3)

Indeed its force is eta_b cross omega, supported in volume O(R),
with the derived core derivatives; the full-space Leray projection
is bounded in unweighted H^s. Standard linearized Euler energy
estimates, obtained by integration by parts of transport and the
bounded Leray multiplier, give for 0<=t<=T

    ||v_exact−v_b||_{H^s}
       <= C_s t exp(C_s t) R^(1/2) log R/R².              (4)

The fixed-profile global velocity derivative bounds make C_s
independent of sufficiently large R. The constant far velocity is
transport only and does not contribute to that growth bound.

The corresponding ACTUAL material displacement solves
xi_t+[u,xi]=v_exact with the same initial xi. The prepared lift obeys
that equation with residual eta_b. The same transport estimate gives

    ||xi_exact−xi_b||_{H^s(D)}
       <= C_s(t+t²) exp(C_s t) R^(1/2) log R/R²,          (5)

on the actual invariant material tube D. It follows directly that
the transported material tag, not a frozen Eulerian section, has the
same controlled displacement and covariance errors on every fixed T.

For n=2, the unnormalized Euclidean quadrupole response has leading
size rho R² times the bending amplitude. Its relative error from
(5) is at most C_T log R/R². Its center-of-mass displacement and
global spin vanish exactly by the azimuthal rank selection rule:
the axisymmetric full Euler operator preserves n=2, whereas global
vector moments have toroidal orders at most one. The actual spin
formula used for this statement is

    delta S_D=rho integral_D [r cross xi_t+2 xi cross u],

with the centered position and all material-boundary terms included
by the exact invariant-domain identity. A nonzero section angular
flux is not substituted for this zero global material-spin row.

The positive inherited quadratic energy of the chosen bending initial
phase is conserved by the full linearized Euler flow. Its positive
initial KKS pair is likewise preserved. These are real same-background
action statements, not an assertion that the observable executes the
restricted two-coordinate oscillator exactly.

In particular (4)-(5) are NOT uniform at T proportional to
R²/log R. At fixed T their absolute quadrupole error can be as large
as the small predicted change of the bending phase. They prove an
actual positive-energy, nearly stationary quadrupole history and a
controlled physical observation map; they do not determine the tiny
optical frequency by comparing two almost unchanged curves. The
weighted full-feedback construction in (2) is the next route for
that stronger conclusion.
