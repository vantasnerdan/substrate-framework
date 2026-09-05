# 0252 — localizable axisymmetric center compatibility

Parent issue200 remains unchanged at843fcb2/v0.182.0. Candidate0248B
asks a localizable compact carrier to retain0211's nonzero elliptic closed
core. This attempt freezes a precise compatibility theorem before a new
symbolic verifier: a smooth axisymmetric steady incompressible Euler field
that is localizable near a nondegenerate poloidal stagnation point at r>0
has zero toroidal velocity there. Hence it cannot carry the specific nonzero
elliptic core required by that candidate. This is a route-level conclusion,
not a no-go for smooth compact Euler, finite invariant tubes, or issue200.

Target_kind: fixed_theorem. One proof route differentiates the two actual
Euler first integrals, angular momentum r*u_theta and speed squared, at the
poloidal center. Their gradients vanish by invertibility of the actual
poloidal velocity derivative. The cylindrical metric factor then decides
the core speed. All fields and pressure are C2 locally, r0>0, and localizable
means u·grad p=0; no analyticity or globally regular Bernoulli boundary is
assumed. The accepted geometry inputs remain unchanged.

Exact analytic calculus and symbolic first-jet identities are the strongest
oracle. Mutation cases retain either a degenerate center or remove
localizability to expose the relevant hypothesis. No sampled eigenvalue,
Morse-index numerical claim, tolerance or empirical comparator is planned.
The script will derive the actual angular-momentum and speed jets, with raw
output/exit captured in Herdr script pane w3:p2.

Pass licenses: compatibility verdict for0248B and candidate generation away
from that matching mechanism. Requires: axisymmetric steady Euler, local
pressure first integral and invertible poloidal derivative. Does not license:
refutation of a compact nonlocalizable carrier, a degenerate core, a shell
with no core, or the full stationary-density objective. The next candidate
keeps the nonzero core while changing the localization mechanism, or uses
an already exact fixed periodic compact-tube field with0250 controls.

## Executed result and independent review

Both exact derivations and all eight symbolic checks pass on first execution,
with recorded exit0. Independent review in review.md establishes the stated
local theorem without correction. Route verdict: established as stated.
Evidence scope: exact local C2 axisymmetric Euler compatibility. Its application
refutes0248B only when the nonzero elliptic core lies inside the pressure-
localizable region. It does not decide outer-only matching, nonlocalizable
compact constructions or the fixed exact Beltrami background response route.
The latter remains active in0250; new outer-only/nonlocalizable matching is
reserved0253. This successful compatibility check activates those constructions
and supplies no campaign exhaustion or terminal PR license.
