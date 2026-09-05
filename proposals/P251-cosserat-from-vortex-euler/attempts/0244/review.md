# Independent review — actual hybrid acoustic supplier

Reviewer: `herdr optical-review pane w3:p3`, a separate fresh Codex process
from the `/root` coordinator, 2026-09-05. I authored no part of 0241 or 0243
and made no implementation or source contribution to this transaction. My
earlier work was the independent 0242 optical review; its supplier license is
used here as an already completed input boundary and is not under renewed
review. No expected conclusion was supplied.

## Decision and strongest supported statement

**Established as stated at the prepared hybrid-acoustic supplier scope.** On
the accepted fixed C016 cell, the source constructs both parity quadratures of
an actual passive Euler response and uses them to correct the complete
point-to-hybrid acceleration row. On every fixed compact time window and in
every preselected finite `C^r` list, its finite smooth-band preparations can
normalize the actual hybrid initial observations to `D,V` and achieve

    U_H,tt = -a |K|^2 (D+tV) + o(|K|^2),
    U_H,D  = D-t^2 a |K|^2 D/2+o(|K|^2),
    U_H,V  = tV-t^3 a |K|^2 V/6+o(|K|^2).

The same construction preserves the first-spatial-order 0231 tag angle and
spin rows. Its finite second-order tagged pressure responses are retained,
not asserted zero. After the actual initial observation map and every new
finite phase, energy and cross jet are fixed, reviewed 0228 supplies the full
action normalization at its conditional same-cell, polynomial-cost,
finite-window scope. This is a prepared linear Euler/Lin sequence, not a new
stationary medium or a parent completion result.

No false or missing load-bearing construction was found in the frozen scope.
No correction is requested; the minimum repair is `none`.

## Frozen transaction and evidence boundary

The criteria are exactly `0244/README.md`. The review was performed on branch
`research/pr199-completion`; the integration head observed at final evidence
capture was `b69839e642fe1dff494c25d7e5387cb694041625`. The source and relevant
captured evidence were:

| Artifact | SHA256 |
| --- | --- |
| `0243/README.md` | `68a6f3eb170fadab5933d86d12376e33089edfaa06bf9573dfb14e2958fdcd33` |
| `0243/hybrid-acoustic-repair.md` | `f6287306c9ea1ffbd445c0cb61e3220316b36d842532b6ab3c83d21d1a10a99a` |
| `0243/verify_controls.py` | `2c8bf1b374d707cac2b24f12c5491d6b2f96891add2f7339a197beb2eb389c49` |
| `0243/first.stdout` | `6a27bab2534953e037781c4a0b8b7dca802fbd0433de70f84c12fd2ba4fd9c7c` |
| `0241/joint-residual.md` | `554e42179f9c24a1b30df3df83919b9e21b1036f5c08cf5686203af386d9a3e0` |
| `src/substrate_framework/euler_joint.py` | `cf4a34d3dfa6bc36f9f7d6cde51d24182243bff881404465985a75c8b5249e0d` |
| `tests/test_euler_joint.py` | `90b1f291eb689c8d1ca306ac3d0a013074b85ac148a4868a0ad8e95d482025e8` |
| `0241/second-pytest.stdout` | `6a05181c1a71c9ada89ff9df19bdb54bacfc244d81fd319ff8c82991e7b5ccee` |

The 0241 material-moment identity and its focused API/tests are used only to
audit equation (1), centroid phase, absolute-velocity dependence, both parity
moment inverses and the initial derivative map. The unrelated 0241 joint
current and geometry residual is not evidence for this supplier. C015/C016,
0234/0242 and 0228/0230 are unchanged dependencies and were not re-reviewed.

## Mean-to-hybrid order audit

The defining observation is not the point mean alone. With `X` the material
centroid, `r=x-X`, and the complete tagged first and second momentum tensors
`A` and `B`, the resolved-minus-centroid momentum has the exact expansion

    Delta J = -i A K - B:(K tensor K)/2
              -i (K.X)(-i A K) + O(|K|^3),

with the appropriate time and input variations of `A`, `B` and `X`. Hence

    U_H,tt = X_pt,tt-rho^(-1) partial_t Delta J

is the correct coefficient to target. This retains shape rate, absolute
velocities, centroid phase and any `K` jets of the moment data. The 0241 exact
point-mass regression specifically exposes the incorrect spin-only reduction,
and its uniform-boost fixture exposes deletion of the second moment. These
tests corroborate the analytic identity; they do not substitute for the
continuous material-tag argument in 0243.

For a fixed baseline, covariance makes the transverse degree-two remainder a
scalar coefficient, while whole-state time reversal makes `b_D` even and
`b_V` odd. Defining those coefficients from the full Euler/Lin and material
observation before choosing controls is therefore noncircular. The required
control targets `-a-b_D` and `-at-b_V` have exactly the two parities supplied
below.

## Actual passive responses and finite-band inverse

For the C016 wrapped-streamline generator `T=omega(c) partial_theta`, the
axial field `exp(-tT)g e_X` is an actual three-component, pressure-free passive
Euler solution at zero Bloch wave number. Direct substitution gives

    g=G(c) sin(theta)  ->  <sin Z exp(-tT)g> proportional to cos(omega t),
    g=G(c) cos(theta)  ->  <sin Z exp(-tT)g> proportional to sin(omega t).

Thus the even response is a genuine second quadrature, not a renamed odd
oscillator. The full whole-state reflection/time-reversal law transports the
background, vorticity and signed initial field together. Positive band
densities and positive whole-field probabilities remain distinct from the
signed coherent weights. Both quadratures are mean free, smooth and have the
same finite quadratic norm displayed in equation (7).

For parity `p`, the derivative map is the matrix with entries

    M[j,l] = integral eta_l(c) omega(c)^(2j+p) dc.

On ordered disjoint positive-frequency supports, multilinearity expresses its
determinant as an integral of a strictly signed Vandermonde in `omega^2`.
Therefore the actual finite-width matrix, not merely its point-band limit, is
invertible. Matching a parity polynomial's initial derivatives is exact.
Parity-polynomial density in `C^r`, the point-band Taylor remainder, and later
band narrowing give the stated finite-window approximation. The weights may
grow without bound as accuracy increases, but every selected preparation has
finite explicit norm; uniform control cost was not claimed.

## Tag separation and initial observation map

The zero-wave sources live on wrapped streamlines below the separatrix, while
the material elliptic tag lies above it. Stationary transport preserves this
separation, and zero initial configuration makes the leading Lin response
vanish on the tag. The complete Bloch pressure correction can first reach the
tag one spatial order later. Consequently:

- the first-spatial-order tag angle, spin, shape and centroid rows are zero for
  both parity families;
- a second-order tagged pressure row is allowed and is retained;
- in `Delta J`, a tagged first-momentum variation carries another explicit
  `K`, while a second-momentum variation carries two. The added control thus
  changes `Delta J` only at order three or higher and cannot alter the
  degree-two hybrid acceleration coefficient.

This is the load-bearing distinction between preserving the first tag rows
and claiming that all higher tag observations vanish. The source makes the
former claim only. Its second-order tag rows remain available to the finite
optical/current target and later joint ledger.

The new source has zero lower mean/tag observation and therefore perturbs the
finite initial `(U_H,U_H,t)` map by `O(|K|^2 C_N)`. After a finite controller is
fixed, that map is invertible for small `K`. Precomposing by its inverse makes
the actual initial observations exactly `D,V`; because the desired
acceleration starts at `|K|^2`, this changes it only by
`O(|K|^4 C_N)` with the fixed finite costs included. Curl adds one more `K`, so
the first acoustic angle changes only at order three. Initial material/current
rows are retained rather than identified with this normalization.

## Precise use of 0228 and common ordering

The new passive controls generally add nonzero action norms and cross forms.
0243 does not infer their cancellation from disjoint tag support. It first
fixes the source, the actual physical observation map, and the complete finite
phase/energy/cross jet, then presents their difference from the desired joint
form to the reviewed 0228 triangular normalizer. That inherited theorem
supplies actual same-cell Kelvin preparations; it includes intrinsic control
jets and all finite cross constraints through degree two. It does not supply
the passive history, the hybrid observation or a missing physical clock.

The imported physical-output conclusion is equally limited: it applies on a
fixed finite window, after all source/preparation norms and the finite
derivative list have been included, and assumes the reviewed polynomial
source-cost/off-tag estimates. For each requested acoustic/optical history
accuracy, the low-frequency controllers and their possibly large costs are
fixed first. They are therefore constants in the later `h` diagonal. A finite
`D_N` can then bound the enlarged list, followed by moment order
`m>2D_N+2`, remote integration order `q`, `|K|=h^(D_N+1)`, and finally small
`h`. This makes the cubic, remote-normalizer, initial-map and clock errors
`o(|K|^2)` without assuming a polynomial norm bound in the acoustic band
width. The ordering is nested, not circular.

## Oracle assessment and findings

The primary oracle is the raw analytic source. The captured 0243 verifier
records 14/14 exact checks for the zero-pressure transport, both angular
quadratures and norms, both finite derivative inverses, the failure of the old
odd-only family on a nonzero even target, and the extra-`K` tag-moment order.
The 0241 focused pytest receipt records five passing tests for the independent
moment and parity APIs. The negative cases expose precisely the two plausible
false greens: silently reusing only the odd family and retaining an order-one
first tag momentum in `Delta J`.

No additional exposing check was warranted. The source identities, the
finite-band determinant proof and the captured negative controls decide the
concrete uncertainties in the frozen criteria. I found no counterexample,
absent actual Euler source, circular target definition, invalid 0228 import or
initial-map singularity within the stated asymptotic scope.

## Precise supplier license and parent status

The current/geometry join may import 0243 as the following supplier and no
stronger one:

> Given the accepted fixed C016 stationary cell, common transverse `D,V`
> inputs, the complete material hybrid observation, the measured positive
> C015 acoustic coefficient, a fixed compact time window and finite derivative
> list, and the reviewed 0228 same-cell form-normalization/error hypotheses,
> there is a sequence of smooth finite-energy actual linear Euler/Lin
> preparations with positive whole-field law whose actual hybrid initial data
> are `D,V`, whose hybrid displacement histories satisfy equations (9)-(10)
> through `o(|K|^2)`, and whose first-spatial-order tag angle/spin rows are
> preserved. Both even-`D` and odd-`V` acceleration corrections are supplied by
> actual passive wrapped-streamline quadratures. All ambient current, material
> moments, finite second-order tagged pressure rows, and complete action/cross
> forms are retained; the latter are normalized through degree two by the
> reviewed 0228 construction under the stated nested cost ordering.

This license does not give uniform control norms, an acoustic-time
`1/|K|` limit, nonlinear finite-amplitude stability, vanishing of every
second-order tag observation, a local boundary couple law, the full optical
current/acoustic-spin/pressure-torque compatibility diagram, stationary
Euclidean EPS geometry or density, an accepted claim, or P251 completion.

- Verification: `symbolic_verified` for the exact transport, parity, moment
  inverse and order algebra, with analytic finite-window approximation and
  diagonal estimates.
- Review: independently audited and supported at the prepared hybrid supplier
  scope.
- Compatibility: compatible extension on the accepted C016 fixed cell; no
  accepted statement is changed.
- Epistemic: established route evidence inside active P251, not canon.
- Route verdict: established as stated.
- Minimum repair: none.
- Correction check: not needed.
- Parent verdict: active; the current/action/geometry join remains open.

Signed: `herdr optical-review pane w3:p3`, separate fresh Codex process and
independent non-author reviewer of 0243, 2026-09-05.
