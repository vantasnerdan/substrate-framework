# Actual acoustic rotation and its independent physical spin repair

Use the SAME accepted C016 field and preparation, not a fitted rotational
row. Let x=(X,Y,Z), u=(psi,v), v=J grad psi, curl u=-u, psi=B cos Y+A cos Z,
A=1/100,B=1. The stationary tag chi(psi) has compact support in its regular
elliptic island and the full axial circle. Its normal centroid is zero,
Q=int chi(Y^2-Z^2) is nonzero, and theta=int YZ delta_chi/Q. Normalize all
densities by the same axial length and cell volume. The reference tag is
unchanged from0216/0218; j0 denotes its actually measured optical G/theta.

## 1. Rigid and rate germs are actual Euler/Lin solutions

Let A_X x=e_X cross x and R_X=A_X x. The infinitesimal rotated stationary
field is w_R=A_X u-Du R_X=[u,R_X]. Differentiating the actual rotated Euler
family proves L w_R=0 with its rotated pressure. For the affine velocity
R_X, direct full-pressure linearization gives

    L R_X-w_R=-2P(A_X u),    A_X u=(0,-grad psi).

The sign in the second identity uses J^2=-I. Thus the right side is zero
on the full Euler pressure problem. Consequently

    xi=(d+t v0) R_X,
    w=v0 R_X+(d+t v0) w_R

solve Euler and Lin exactly as first-affine germs. These are not separately
periodic rotations of a fixed lattice: their affine terms are the first
laboratory Bloch derivatives of the common translation/velocity columns.
Only that first-gradient statement is used below.

The accepted preparations realize these germs, rather than merely admitting
some other symmetry family. In the actual C015 API take the antisymmetric
contraction

    prepared(kappa=e_Y,D=e_Z)-prepared(kappa=e_Z,D=e_Y).

Its d=kappa_Y D_Y-kappa_Z D_Z is zero, its symmetric current is zero, its
negative-helicity return is zero, and its material-rate correction is zero.
Its first velocity lift equals A_X u. These are exact identities of the
complete finite Fourier fields, including the higher-harmonic corrector.

For C016 the full same-lift velocity source is
    b_V=-P[(kappa.u)V+q_V].
The same antisymmetric contraction gives -2P(A_X u)=0 and d=0. Each actual
passive velocity/configuration return carries the symmetric coefficient
kappa_Y V_X+kappa_X V_Y and hence contributes zero to this contraction.
The zero initial remainder therefore remains zero under the full Euler
evolution. No pressure row or subsequent history is dropped.

For xi=Omega R_X, integration by parts gives delta Q_YZ=Omega Q, hence
the literal theta=Omega. Rigid rate gives theta_dot=Omega_dot. This is a
normal covariance measurement, not an assigned phase coordinate.

## 2. The measured ensemble Gram removes the apparent factor three

For whole-field rotations R and axial scalar detectors, use the already
declared physical estimator

    Phi=M^-1 E[n theta],    M=E[n n^T]=I/3.

Here the parity of n theta is axial; improper rotations transform the
complete detector/frame law, not n alone. This estimator is the unique
least-squares reconstruction of a common vector from its actual scalar
projections. It returns a rigid vector Omega exactly because theta=n.Omega.
It does not multiply the actual spin by three.

More generally the complete body first-gradient response is a linear row
T_ab(t) B_ab, B=grad U. Haar averaging of its whole-frame push-forward uses
E[R_i1 R_ja R_kb]=epsilon_ijk epsilon_1ab/6. Thus only the antisymmetric
body contraction T_ZY-T_YZ survives. Its value is exactly one for D and
t for V by section1. Symmetric strain may produce a body tag response but
has zero averaged axial vector. Therefore the actual prepared acoustic
observation is

    Phi_ac=curl(D+t V)/2+O(|K|^2)

on every fixed time window. Under the complete inversion-symmetric law,
the polar-to-axial observable is odd in K, so its next analytic spatial
term is cubic, not quadratic. For the finite smooth preparations the
ordinary full-pressure remainder bounds of C016 apply. This is not an
acoustic-time limit or a statement about an unprepared Euler ensemble.

## 3. Rigid inertia is not optical inertia: construct the missing spin

The literal rigid displacement dipole and mechanical spin are

    G_rigid=I_tag Omega,  S_rigid=I_tag Omega_dot,
    I_tag=rho int chi(Y^2+Z^2).

The background axial tag spin is unchanged by an axial rigid rotation.
There is no general identity I_tag=j0. Matching the optical measured mass
by simply renaming this rigid inertia would change the physical claim.

Use instead the actual smooth mean-free stationary vector

    b=(0,v),   [u,b]=0,   L b=0.

Both facts follow directly: T psi=0, and v.grad v is a gradient by the
actual stationary two-dimensional Euler equation. Equivalently
(psi,(1+epsilon)v) is a stationary Euler family with its pressure scaled
by (1+epsilon)^2. The background axial velocity remains fixed. This b,
rather than the full u, also has zero chi-weighted normal mean and no
axial material drift.

For any affine Omega(t),

    xi_b=c Omega(t)b,    w_b=c Omega_dot(t)b

are exact Euler/Lin solutions. They change the velocity for the rate
column; this is an allowed actual acoustic initial velocity, not claimed
to be a fixed-circulation Kelvin tangent. C016 itself allows such initial
velocity/configuration data. Their tag variation is identically zero
since b.grad chi=0, so their centroid, covariance and angle remain zero.
Their actual moments are nonetheless

    G_b=c S0 Omega,  S_b=c S0 Omega_dot,
    S0=rho int chi r cross v=rho int chi r.grad psi.

S0 is nonzero without a numerical assumption. In this elliptic island
psi has a strict maximum; choose F'(psi)=chi with F=0 below the tag's
outer level. Then F is nonnegative, positive on a set of nonzero area,
and compact in the island. Integration by parts gives

    S0=-2rho int F<0,     int chi v=int J grad F=0.

Thus the explicit finite coefficient c=(j0-I_tag)/S0 constructs

    theta_ac=Omega,    G_ac=j0 Omega,    S_ac=j0 Omega_dot.

The same positive whole-field law yields j=j0/3 and the projected internal
rows E[n S_X]=j Phi_dot, E[n G_X]=j Phi.
The added fields need no signed probability or changed tag fraction.
Time reversal changes S0 and b together, preserving the displacement
normalization; the complete data, rather than one amplitude alone, are
pushed forward under reflection and reversal.

## 4. What changes downstream and what is still separate

These projected rows are not silently equated with a full three-dimensional
tag moment. This tag wraps the axial circle: normal cross-section moments
are literal, while a full axial lever arm requires its cell-boundary flux.
Transverse spin/dipole components under generic affine input must be
combined with those flux and ambient rows before applying the full hybrid
current identity. The scalar detector construction alone does not remove
them by symmetry. The boundary-current continuation0232 treats that exact
localization distinction.

The added b-return itself has zero transverse full-spin components in
a centered axial cell: they contain the factor int X dX=0, and b has
no axial velocity or axial boundary motion. Its full added spin is
therefore genuinely n c S0 Omega_dot before averaging. This proves a
usable full-vector adjustment direction, but does not determine the
unrepaired baseline's transverse response. Averaging that baseline
requires the antisymmetric contractions of all three body spin rows,
not only the one body-X angle detector. The missing row is thus explicit.

This new acoustic preparation leaves the complete Euler velocity mean
unchanged through first K: b is mean free and stationary. It leaves the
tag density and its symmetric shape exactly unchanged at that order.
Its full material/ambient dipole and the hybrid-current reconstruction
are retained as separate actual rows. Since the added G itself has order K,
its curl first enters
physical U at order K^2. That is a second-gradient chart correction, not
a change of the leading acoustic acceleration: its leading time history
is affine, and its second time derivative vanishes. Normalize the actual
initial U and U_dot by that near-identity second-order physical map;
the acoustic squared-frequency coefficient a is unchanged. Corrections
to curl U in the angle row then first enter at order K^3.

The added actual phase, Jacobi energy and acoustic/optical cross forms
are finite physical forms and are not set to zero. The reviewed0228
normalization applies to their resulting joint two-jet only when its
actual source-cost/output hypotheses are met; here b and c are fixed
smooth finite-cost sources on the same field. No normalizer is used to
manufacture their angle or spin. Existing optical second-jet construction
and transfer to stationary EPS geometry remain separate active tasks.

This result supplies the actual first-gradient acoustic angle and projected
G/S row, including the previously missing measured axial spin repair.
The complete hybrid angular-current join is still required by0227.
It does not claim a completed same-field optical/EPS continuum merely
from a symmetry tangent or the conditional branch matrix.
