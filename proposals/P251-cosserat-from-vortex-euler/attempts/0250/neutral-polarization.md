# Neutral constant-curl polarization and its finite-N Kelvin boundary

This continuation checks the proposed leading Euler polarization on the
fixed-background Route C and identifies exactly what it does and does not
construct.  Notation here follows the proposal: `a` is the leading velocity
amplitude and `b=k cross a` its leading vorticity amplitude.  This is the
opposite letter assignment from equation (10) of `fixed-ring-response.md`.

## 1. Exact sign and normalization

Let `A=Du` along a trajectory of a steady incompressible Euler field, so
`tr A=0`, and use the physical convention

    A-A^T=[omega cross].                                  (1)

Let an invariant phase `I` satisfy `D_t I=0`, put `k=grad I`, and take a
transverse leading velocity amplitude

    k_t=-A^T k,
    a_t=-A a+2k(k.Aa)/|k|^2,
    k.a=0.                                                (2)

The last term in (2) is the full-pressure principal symbol.  Define
`b=k cross a`.  Direct differentiation gives

    b_t=-(A^T k) cross a-k cross (A a).                   (3)

For `tr A=0`, the cofactor identity is

    -(A^T k) cross a-k cross (A^T a)=A(k cross a).        (4)

Using `A a=A^T a+omega cross a` and
`k cross (omega cross a)=omega(k.a)-a(k.omega)` in (3) yields

    b_t=A b+(k.omega)a.                                   (5)

Thus the proposed sign is correct.  If `curl u=omega=lambda u`, with
nonzero constant `lambda`, and `k.omega=0`, then `k.u=0`.  The choice

    b=u,
    a=-k cross u/|k|^2                                   (6)

has exactly `k cross a=u`; no factor of `lambda` is missing.  Along a
steady trajectory, `u_t` in material notation is `D_t u=A u`, so (5) is
satisfied.  Since cross product with nonzero `k` is invertible on `k`'s
orthogonal plane, equations (5)--(6) imply (2) as well.  If instead one
normalizes the vorticity amplitude to `b=omega`, the velocity amplitude is
`a=-lambda k cross u/|k|^2`.

Multiplication by any transported scalar `q`, `D_t q=0`, preserves the
construction:

    b_q=q u,
    a_q=-q k cross u/|k|^2.                              (7)

On an exact action-angle region one may take
`q=exp[i ell.(theta-Omega(I)t)]`.  Whenever
`nu(I)=ell.Omega(I)` is nonconstant, (7) supplies the proposed real scalar
transport-clock channel at leading Euler-symbol order.  Real and imaginary
parts give its two phase quadratures.  No frequency approaching zero is
used.

## 2. The leading Kelvin symbol is singular, but an exact finite-N lift exists

The polarization (6) is not a regular leading-order Kelvin displacement at
the exactly normal covector `k`.  For an oscillatory divergence-free Lin
displacement with leading amplitude `xi`, the Kelvin velocity symbol is

    P_k(xi cross omega),
    k.xi=0.                                               (8)

When `k.omega=0`, the vector identity

    k cross (xi cross omega)=(k.omega)xi-(k.xi)omega      (9)

shows that `xi cross omega` is parallel to `k`; (8) therefore vanishes.
Equivalently, the leading vorticity variation
`k cross (xi cross omega)` is zero.  Hence treating (6) as an ordinary
order-one Kelvin column would give a false nonzero action normalization.

The total-phase construction below displays the apparent one-power loss, but
it is not the minimum exact repair.  The subsequent
`isovortical-channel.md` constructs a divergence-free generator directly:
its normal part maps to (6) exactly and a tangential `beta u` part solves the
cohomological divergence equation without changing `xi cross omega`.  The
present calculation remains useful because it identifies the same source-
norm scaling in principal-symbol variables.  Use the exact transported phase

    Psi_N=N I+ell.theta-nu(I)t,
    D_t Psi_N=0,                                         (10)

whose total covector is

    K_N=N dI+ell.dtheta-t dnu,
    K_N.u=nu(I),
    s_N=K_N.omega=lambda nu(I).                          (11)

On a subband separated from `nu=0`, `s_N` is nonzero even though
`N^-1 K_N` tends to `k=dI` and `k.omega=0`.  The principal Kelvin inversion
for a transverse velocity amplitude `a_N` is then

    xi_N=K_N cross a_N/s_N.                              (12)

For (6), equations (11)--(12) give

    xi_N=N u/[lambda nu(I)]+O(1).                        (13)

The sign in (13) is positive:
`k cross[-k cross u/|k|^2]=u`.  Thus the total-phase inverse loses one power
of `N`.  The exact construction in `isovortical-channel.md` realizes the same
leading velocity with an `O(N)` cohomological tangential correction, and
closes the finite-`N` isovortical initial-data question on an open subannulus
of the exact 0211 ring.  Observation normalization and the common preparation
diagonal must still include that loss.  The limit polarization itself lies
at a singular boundary of the normalized Kelvin symbol; it does not by
itself supply a nondegenerate KKS pair or its Jacobi energy.

## 3. What remains on the final periodic field

The calculation removes the need to assume a stable two-dimensional Floquet
band for this scalar controller: one neutral Cauchy polarization is enough at
leading order.

It applies literally on the exact 0211 axisymmetric ring before any periodic
approximation.  On a compact regular inner annulus take `I` to be its physical
section-flux action.  The exact flow is tangent to `I`'s levels, `grad I` is
nonzero there, and 0211's strict flux twist makes an integer combination of
the poloidal/toroidal transport frequencies nonconstant on a smaller annulus.
After avoiding a zero of that combination, equations (6)--(13) therefore give
an exact leading scalar Euler-symbol band on one already fixed compact ring.
This is stronger than a transported tracer: its leading velocity and
vorticity amplitudes obey the full-pressure Euler symbol.  The later
`isovortical-channel.md` solves the finite-`N` divergence-free Kelvin
preparation and all-order fixed-time pressure expansion on this exact ring;
its remaining burden is the literal measured gain and normalizer matrix.

The normalized principal Kelvin symbol still vanishes and the source carries
the one-order loss displayed by (13), but the exact cohomological generator
in `isovortical-channel.md` closes the initial isovortical leaf and permits
0112's arbitrary-order pressure recurrence to start from the actual Leray
field.  That does not prove that the physical gain matrix or full KKS matrix
has the desired rank.

The construction is harder again on the final 0145/0147 periodic approximant.
It does not yet produce an actual all-order quasimode on that field.

The exact axisymmetric precursor has a smooth invariant action and makes
(10) literal.  A generic rational-direction periodic Beltrami approximant is
only close to that integrable field.  KAM persistence gives selected invariant
tori, not automatically one smooth open action chart on which compact smooth
Euler preparations and the band synthesis can be varied.  A finite-order
normal form `I_q` has an eikonal defect multiplied by the large carrier:

    D_t exp(i N I_q)=iN(D_t I_q)exp(i N I_q).             (14)

Therefore finite-order action persistence is sufficient only with an
explicit joint choice of normal-form order, packet width and `N` making (14),
the full pressure hierarchy and the singular Kelvin loss (13) smaller than
the physical gain divided by the finite band-synthesis coefficient cost.
One unit Floquet multiplier at leading order does not prove this hierarchy,
a spectral band or stability.

The fixed positive tag ladder in `fixed-ring-response.md` can absorb an
additional polynomial power of `N`: increase its WKB order before selecting
the sparse carrier.  But its nonzero leading coefficient must now be computed
from the actual finite-`N` Kelvin column (12), not from (6) alone.  The same
calculation must expose, on one common preparation:

1. both real scalar parities and the literal acoustic point-to-hybrid
   acceleration gain;
2. tag/ambient separation and all initial observation rows;
3. the full KKS/Jacobi-energy cross matrix after the one-power displacement
   loss (the two transport quadratures have exactly zero mutual KKS entry);
   and
4. the physical `G`, spin/current, pressure and current-memory rows needed by
   C-CST-017's joint normalizer.

For a stationary material observation the fixed ladder should depend only on
the invariant action, for example
`chi(I)=chi_0(I)[1+epsilon sum_j a_j cos(N_j sigma(I))]`.
It is then transported trivially and stays positive.  The laboratory position
or covariance moment, rather than a co-moving angular tag phase, must supply
the `ell` harmonic in (10); putting the same transported angular phase in both
the tag and the perturbation would cancel the desired clock.  The leading
gain is consequently the actual Fourier coefficient of the physical centroid
(`ell=1`) or covariance/spin (`ell=2` and companions) weight against (6).
Nonvanishing and joint rank of those coefficients are not consequences of
the polarization identity and remain to be computed.

## 4. Route statement

`route_verdict: established at exact leading-polarization scope; the later
exact isovortical continuation repairs finite-N preparation on the 0211
ring, while the fixed-field supplier remains open at gain/normalizer and
final-periodic normal-form closure`

`evidence_scope: exact full-pressure Euler symbol identity, correct
constant-curl normalization, transported scalar channel, and explicit
one-power singular Kelvin lift`

The exact 0211-ring continuation and its remaining construction are recorded
in `isovortical-channel.md`.  Transfer to one fixed final periodic background
still needs a controlled invariant-action normal form, followed by the
complete measured gain and action/current calculation.  This is a positive
new controller polarization and a strictly smaller missing mechanism; it is
not a parent no-go and does not alter the separate 0248/0253 geometry
transaction.
