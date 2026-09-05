# Independent review of 0266 cellwise phase and acoustic Schur gain

Reviewer: Herdr geometry-review pane `w3:p4`, fresh non-author.  I reviewed
the active 0266 source and receipt at hashes
`6cf3dea30789a39c02f910fcb60ba01ab7fc85e911b5fc62e6d22ce332308794`
and `6ba77b965e3d60db294f8c48494803d0e0e379e5e2ad3ac4497b3e5d1bcf6c59`.
The corrected exact oracle is
`180e49da24de7f2ecc238a4778e8a6fba50c809079cfec59c75d531ac772779a`;
its recorded exit is zero.  The initial structural-equality stop is an
implementation failure, not contrary equation evidence: the active script
replaces Python structural equality by simplification of the symbolic
difference without changing the expected formulas.

## Verdict and strongest supported statement

`route_verdict: established at the conditional whole-law finite-window
Schur-gain scope, with one non-load-bearing full-tensor wording correction`

Assume the actual compact periodic background and separated-cell chart, and
import the 0260/0250 axisymmetric `n=0,m=+/-1` local controls at exactly their
recorded scope.  For a fixed finite observation/time-derivative inventory
and a fixed sufficiently small nonzero macroscopic wavevector (or retained
ray), the cellwise rephasing is an exact admissible Bloch preparation.  The
carrier and pressure order may then be chosen after that wavevector and the
finite normalized source costs.  Under this hierarchy, the periodic
image/mean contribution is `o(|K|^2)` on the retained window, the averaged
first-`K` acoustic-to-`(theta,G,S)` row is zero, and the complete parity-even
oblique material tensor gives

    E[(D dot n) P_K(delta B:K K)]
       =|K|^2 (4 b_perp+2 b_parallel)D/15
       =|K|^2 4 i pi^2 C V^2(R^2+s^2)D/(15s)             (1)

for `D dot K=0` in the circular leading model.  After the common physical
weights and phase are restored, this coefficient is nonzero on a suitably
chosen small core annulus.  Smooth fixed-core convergence preserves it for
one sufficiently large fixed finite `R`.  Thus the first odd lower-left row
`c_1` in 0250's Schur matrix vanishes for this whole source law and
`a_Sch=a_2` is nonzero.  No directionwise scalar division near `n parallel
K` is used.

This statement is conditional on the actual compact background and local
controls.  It is a finite-window, ordered high-carrier pressure estimate,
not exact pressure decoupling, global pressure localization, an infinite
lattice limit uniform through `K=0`, or parent completion.

## Exact source and phase audit

With `div_K=div+iK dot`, the proposed envelope obeys the exact identity

    div_K[e^(-iK dot(x-X_c)) xi_0]
       =e^(-iK dot(x-X_c)) div xi_0=0.                    (2)

The positive zero collar makes this a smooth periodic envelope even when
`K` is not reciprocal-lattice valued.  Multiplication by the physical
Bloch factor `e^(iK dot x)` leaves the source in cell `m` equal to the fixed
local source times the constant phase `e^(iK dot(X_c+Lm))`.  Hence, in whole
space, conjugacy of the Helmholtz projector and of the linear Euler/Lin
equations gives the fixed single-ring response up to that constant phase.
For the periodic application the response is not exactly cell-decoupled:
the local singular kernel is the whole-space one, while the image kernel and
the retained ray-wise zero mode are smooth tests against the high-action
carrier.  The 0262/0265 integration-by-parts license gives arbitrary fixed
inverse-carrier order, including any predeclared finite number of time and
observation derivatives.  This supports precisely the ordering stated
above and no stronger uniform-in-`K` pressure claim.

For any collapsed central observable `O`, direct differentiation gives

    delta[e^(-iK dot X)O]
       =delta O-i(K dot delta X)O_0.                      (3)

There is no additional internal first-`K` source jet after (2).  The
axisymmetric acoustic column has `delta X=x_0 a_n n`, zero covariance tilt,
and exact zero `G` by the actual finite-carrier `G`-canceling ratio.  The
stationary invariant tag has zero baseline centroid velocity and baseline
`G_0=0`.  For covariance, the phase term is a scalar multiple of the
uniaxial reference covariance `Q_0`; the simple-eigenaxis angle differential
annihilates that scaling direction.  For spin, `S_0=s_0 n` and the only new
centroid-phase row is proportional to

    a_n n (K dot n),                                     (4)

where the polar acoustic amplitude is `a_n=(D dot n)` times an even scalar
function of `K dot n`.  Its sphere average is zero under `n -> -n`.  The
intrinsic zero-`K` polar-to-axial spin row separately cancels in the paired
reflected law.  These facts establish the whole-law `c_1=0`; they do not say
that every local realization has zero spin.

## Full oblique tensor audit

The active Cartesian calculation correctly derives the parity-even entries

    b_parallel=delta B_zzz=i pi^2 m C V^2 s/2,
    b_perp=delta B_zxx=delta B_zyy
                 =i pi^2 m C V^2(R^2/s+3s/4),
    delta B_xzx=delta B_yzy=-b_parallel/2.                (5)

Consequently

    delta B:K K
      =n[b_parallel K_n^2+b_perp|K_perp|^2]
         -b_parallel K_n K_perp,

and projection removes the component parallel to `K`, giving

    P_K(delta B:K K)
      =P_K n[b_perp|K_perp|^2+2b_parallel K_n^2].         (6)

The independent sphere identities
`E[n_x^2]=1/3` and `E[n_x^2 n_z^2]=1/15` yield (1).  The mixed entries in
(5) are load bearing: deleting them changes the coefficient of
`b_parallel`.  Finally, `mC V^2=NsV/lambda` has one fixed real sign on the
chosen Bessel-core annulus.  A common complex phase and conjugate real
quadratures therefore turn (1) into a real nonzero physical row.  The
`R^2/s` term has a strict nonzero scaled limit, while the parallel, mixed,
chiral, lower-carrier, and smooth finite-core corrections are lower after
division by `R^2`.

## Exact unsupported extension and minimum correction

Equation (7) of 0266 and the receipt call the oracle's expression the
literal central material tensor with all 18 symmetric entries.  That is too
strong.  If

    A_ij=integral rho chi u_i x_j

at the centered reference tag, differentiation of the defining central
moment gives

    delta B^cen_ijl=delta B^(7)_ijl
                    -delta X_j A_il-delta X_l A_ij.       (7)

The stationary axisymmetric tag has `A+A^T=0`, and axial symmetry leaves
only its axial-spin pair `A_xy=-A_yx`.  Since this acoustic source has
`delta X` parallel to the ring axis, (7) does not alter any entry in (5),
(6), or (1).  It does alter the chiral entries with one `z` and one
transverse moment index, including the displayed local `B_yzx`; the change
is of the same `O(R)` order as that chiral entry.  Those terms are
reflection-odd and their whole-law contribution vanishes.  Thus this is a
real false local-tensor extension but not a missing load-bearing step for
the averaged acoustic Schur gain.

The minimum correction is to insert (7) before describing the full tensor,
and to describe the current oracle as computing the uncentered leading
expression whose parity-even sector is also the central sector.  To restore
the stronger all-18 claim, extend the oracle with `delta X` and the baseline
`A` tensor, show the changed chiral entries explicitly, and verify their
reflected-law cancellation.  No new source family, gain parameter, or
pressure construction is required.

## Supplier boundary

The precise 0266 supplier license is:

> Conditional on one actual compact stationary periodic background with the
> recorded invariant core/source gap, and on the exact 0260/0250 local
> axisymmetric controls, cellwise Bloch rephasing produces a physically
> constant-on-each-ring source.  For each retained finite window and finite
> derivative inventory, with the macroscopic wavevector chosen before the
> high-action carrier and pressure order, the full periodic pressure-image
> error is below quadratic order.  The isotropic/reflected whole law has zero
> first-`K` acoustic-to-`(theta,G,S)` row and a nonzero full-oblique
> transverse quadratic hybrid gain on one sufficiently large fixed ring.

This review does not re-review or establish 0263 compact-field existence,
extend 0260/0250 beyond their consumed controls, control modes or moments
beyond the registered finite inventory, supply the separate KKS/Jacobi and
current normalization, or conclude issue200.

## Post-review bounded correction check

The 0267 minimum scientific correction is now closed.  I checked the amended
0266 source at
`e5c3e641c9f9614875e57861981182bdc80d3f4bb0d5e1b8b617115a8ff88765`
and the active centered oracle at
`22002752fcee1e333f2a23b4966bf45058e8ec22b725cbd37ea384fc46b37e0e`.
An independent execution of the active oracle returned exit zero.

Amended equation (7) is the exact central variation (7) above.  The oracle
first derives the shell diagnostic

    A_0,xy=-2 pi^2 R W,       A_0,yx=2 pi^2 R W,

and then, correctly, does **not** identify that shell diagnostic with the
radially integrated physical values.  It retains the latter as the
independent symbol `a_spin` and retains the independently integrated
centroid variation as `delta_X_z`.  It subtracts both center terms from all
tensor entries.  In particular,

    delta B_yzx=-pi^2 C R V W-a_spin delta_X_z,            (8)

with the opposite signed correction in its paired chiral entry.  Reversing
the background swirl also reverses `a_spin`; hence these additions disappear
from the parity-even tensor.  The rerun reproduces, unchanged,

    b_parallel=i pi^2 C V^2 s/2,
    b_perp=i pi^2 C V^2(R^2/s+3s/4),
    b_mixed=-b_parallel/2,

and the whole-law coefficient
`4 i pi^2 C V^2(R^2+s^2)/(15s)`.  The original uncentered implementation and
output are preserved separately.  This is exactly the minimum repair named
by the independent review; no new physical hypothesis or tuning parameter
was introduced.

`post_correction_verdict: 0266 is established as stated at the precise
conditional whole-law finite-window supplier scope above`

One evidence-record repair remains outside 0267's write authority.  The
current 0266 `receipt.md` still records the pre-correction hashes
`6cf3dea...` and `180e49d...`, rather than the active hashes above, and its
prose does not distinguish the newly preserved uncentered oracle from the
active centered one.  Updating that receipt to identify the active source,
oracle, stdout (`eb80ace5...`), zero exit, and preserved predecessors is the
minimum provenance repair.  It does not reopen the scientific verdict.  The
same-field finite-`K` history/action/current theorem being assembled in 0265
is a separate joint-review boundary and is not assessed here.
