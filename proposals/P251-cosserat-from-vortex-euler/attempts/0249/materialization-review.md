# Materialization review — proposed C-CST-017

Reviewer: `herdr optical-review pane w3:p3`, separate from the `/root`
coordinator, 2026-09-05. This is the bounded materialization/correction check
requested after the completed 0245 and 0247 reviews. I did not author the
candidate payload, the additive API, its tests, or either reviewed source. The
scientific supplier decisions are unchanged and were not substantively
reviewed again.

## Decision

**The proposed C-CST-017 statement, dependencies, evidence mapping and
additive `prepared_joint_symbol` implementation match the completed 0245 and
0247 support boundaries.** I found no scope extension, dependency omission,
formula mismatch, or implementation defect inside this bounded transaction.
The minimum correction is `none`.

This file checks a candidate payload only. `review: accepted`, prospective
`accepted_in: v0.182.0`, and the proposed campaign provenance do not acquire
registry or release authority until the later central claim transaction
materializes and validates them. The compact Euclidean same-field
geometry/density objective remains open.

## Frozen artifacts and hashes

The materialization boundary inspected here is:

| Artifact | SHA256 |
| --- | --- |
| `0249/README.md` | `ae5a8d95cd79c6f6488774eb412b21dd6fea1123b9e1b7045c792a4ddbcdb6db` |
| `0249/claim-payload.yaml` | `c57754f13c46f84539132ba92ee00dc749e152486bc8888798b07682dd82aa6b` |
| `src/substrate_framework/euler_joint.py` | `637fb63102968f7413c46a4a5afe1af2efd5f47888686ba9cd3e4879854eab33` |
| `tests/test_euler_joint.py` | `cd528b8e85e16dbbaf625e67d4446498b2b61271e4a56d3ead24d871d76c2b19` |
| `0249/repaired-pytest.stdout` | `4f3cbdb56c8f11f41c6ee9b54fe30da73500945601b4ec57e01e84ace366f02a` |
| `0249/repaired-pytest.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0249/test_initial_snapshot.py` | `830ee8566b24b5c24b6567b589def78b88204f9ccd4db58885641744b3702660` |
| `0249/first-pytest.stdout` | `e65210151273fe840b1d75eeed51053a149f440e66ae8ddd866c94593a5fd49f` |
| `0249/first-pytest.exit` | `4355a46b19d348dc2f57c046f8ef63d4538ebb936000f3c9ee954a27460dd865` |
| `0249/oracle-repair.md` | `7d304a333d9dc3406f9b03b1a8aab021abefeb31612002e79ee1184f8ee781c4` |
| `0245/README.md` | `db11a9faf813a135650beabc4d767d51a4f11f6c288bdc4a7e822733771d55d2` |
| `0236/periodic-core-and-density.md` | `3ec7e27792e64bdab873c783e680750275c13b7769c6832931caf2d242a9792a` |
| `0247/review.md` at checkpoint `1edefee` | `43d98ae26648083d893d810b5d9bebbf9bebae97402f3298f42cead7973102d0` |

The working branch was `research/pr199-completion` at committed checkpoint
`1edefee`. No source, payload, registry, release or central record was changed
by this review.

## Statement mapping

The first five substantive paragraphs of the proposed statement are the
reviewed 0247 periodic supplier, with the same quantifiers and exclusions:

- one positive whole-field O(3)/time-reversal law and one summed actual
  Euler/Lin preparation per realization, selected along a common ordered
  long-wave sequence;
- the actual physical map
  `T=[[I,-j C/(2rho)],[C/2,I]]`, the acoustic/optical branch equations and
  `e,e_t,e_tt,r_A,r_B=o(|K|^2)` on each fixed compact window and finite
  derivative inventory;
- `mu=rho a`, `alpha=j nu^2/4`,
  `gammaT=j(cT-alpha/rho)`, `gammaL=j cL`, the exact cubic residual and the
  bounded physical-state inverse;
- both inherited action forms pulled through the same physical map, rather
  than an energy-only match;
- literal material-plus-ambient momentum/angular balance, including pressure,
  convective transport, `Q`, initial angular charge and `q_t` memory, followed
  by periodic bulk virtual-work equivalence rather than pointwise traction
  equality.

The payload correctly makes `cT,cL` prepared output targets, not universal
unprepared moduli, and retains `cT>j nu^2/(4rho)`, `cL>0`. Within the explicit
constant-density Euler assumption, its phrase “ambient density” denotes the
same `rho` carried by the material tags and their continuous complement; it is
not a tube filling fraction or an ambient-only mass. The final exclusions
prevent either reading.

The geometry paragraph is exactly the independently reviewed 0245 periodic
supplier. It retains the actual C016 background, nondegenerate elliptic core,
physical flux action, nonzero flux twist

    r_J(0)=-235025/2060602,

normal torsion `4 pi^2 r_J(J_b)`, the local analytic divergence-free/Moser
persistence hypothesis, positive periodic tube fraction and the precise
whole-law measured-coefficient density. It also states the decisive topology:
the core is noncontractible in `T^3` and lifts to an unbounded line. It does
not promote that object to a compact Euclidean ring or claim a nearby
stationary Euler family. `0245/README.md` is the signed geometry review artifact
itself, not merely an unexecuted preregistration, so the evidence edge in the
payload points to the completed decision.

The last paragraph preserves the shared boundary of both reviews: no
all-`K` invariant manifold, acoustic-time uniformity, nonlinear finite-
amplitude result, universal modulus, globally fixed switching-scale
preparation, compact Euclidean geometry/density, or parent completion.

## Dependency and evidence mapping

The declared accepted dependencies have the correct roles:

- `C-CST-009` supplies the accepted conditional second-gradient physical
  action and periodic bulk balance conventions used for the canonical
  representative;
- `C-CST-015` supplies the exact stationary C016 field and positive
  displacement coefficient/preparation;
- `C-CST-016` supplies the same-cell common-`D,V` acoustic action and
  fixed-window diagonal, and already depends on C-CST-015.

The later acoustic, optical, normalizer, current and periodic-geometry
constructions are evidence proving the new claim rather than already accepted
claim dependencies. The payload maps their roles without treating proposal
chronology as authority: raw 0241/0243/0236/0246 sources establish the new
identities and constructions; the 0242, 0244, 0245 and 0247 records supply the
independent scope decisions; the symbolic source/output pairs corroborate the
finite algebra; and the 0249 API/tests materialize the reusable residual.
References inside the signed reviews retain the deeper 0228/0230 action and
0232/0235 current evidence at their exact conditional scopes. No compact
Euclidean evidence is cited or implied.

The proposed evidence list includes the geometry source and its signed review,
the joint source and review, both new supplier proof/oracle groups, and the
additive API with its repaired direct test receipt. That is sufficient for the
statement actually proposed. The preserved initial test failure remains
attempt provenance and need not be mislabeled as positive claim evidence.

The four status axes are mutually consistent for a proposed acceptance
transaction: exact analytic derivations with symbolic corroboration support
`symbolic_verified`; the completed independent decisions support `accepted`;
the result extends rather than changes C-CST-015/016, supporting
`compatible_extension`; and `active` preserves the open Euclidean parent
frontier. These statuses remain prospective until registry materialization.

## Additive API implementation check

`prepared_joint_symbol` is a faithful reusable extraction of the already
reviewed 0241/0247 algebra for the reference axis `K=(0,0,k)` and retained
coordinates `(Ux,Uy,Phix,Phiy,Phiz)`:

1. It sets `alpha=j nu^2/4` and selects the legitimate periodic curvature
   representative
   `c_tr=j cL/2`, `c_s=0`, `c_a=j(cT-alpha/rho)`. Therefore the canonical API
   receives exactly `mu=rho a`, `gammaT=j(cT-alpha/rho)` and `gammaL=j cL`.
2. It constructs `C=i[K cross]` with the reviewed Fourier sign and derives the
   complete physical observation matrix, mass matrix and branch-frequency
   matrix in the correct five-coordinate order.
3. It obtains stiffness from the unchanged canonical
   `micropolar_fourier_stiffness` API and computes
   `defect=K2*T-M*T*D` directly. The expected cubic matrix is not inserted as
   an input, and the exact nonzero cubic remainder remains exposed to callers.
4. It returns immutable matrices and documents that actual positivity,
   supplier existence and analytic error bounds are caller licenses. Permitting
   symbolic signs is appropriate for counterexample probes; the function does
   not claim arbitrary parameters construct an Euler ensemble. It rejects the
   algebraically singular `rho=0` case.

The direct tests independently reconstruct the real cosine/sine-quadrature
energy and compare it with the Hermitian Fourier stiffness, check signed-`k`
reality, the zero-wave limit and observation determinant, substitute arbitrary
history/error functions into the full residual, require every defect entry to
be zero or exactly cubic, retain a nonzero generic defect, check the special
zero-defect parameter relation, and expose the wrong optical observation
factor. These predicates cover the implementation delta rather than repeating
the supplier existence proof.

The final API hash includes the later overview-docstring clarification that
removes the transient word “Unpromoted” and names the physical branch residual
among the module's tools. Comparison with the previously inspected source
shows no function, signature, return value or mathematical statement change.
It therefore does not stale the banked direct-test receipt.

The captured repaired execution reports `8 passed` and exit `0`. The first
execution is correctly preserved with seven passes, one failure and exit `1`.
Its failed line formed a SymPy matrix of Python booleans and compared it to an
integer matrix of ones even though every computed degree predicate was true.
The repaired test replaces only that invalid container comparison with Python
`all(...)`; `oracle-repair.md`, the initial snapshot and both outputs preserve
the diagnosis. The API source did not change in response, so the initial error
is a test-representation failure, not contrary residual evidence.

No rerun was needed: the frozen captured outputs, exit receipts, exact source
and direct algebra decide this materialization check. No full suite was run.

## Materialization license and remaining dependency

The later central transaction may materialize this exact payload and API with
no scientific scope correction, subject to its ordinary registry/release and
changed-surface validation. It may not use this review to assert that the
candidate payload is already registry authority.

- Materialization verdict: support boundaries match.
- API verdict: faithful additive implementation of the reviewed joint symbol.
- Evidence scope: exact prepared periodic history/action/current and periodic
  tube/flux-twist/density conjunction, with symbolic implementation coverage.
- Minimum correction: none.
- Scientific correction check: complete; no second review was performed.
- Remaining dependency: construct and independently review the compact
  Euclidean same-field geometry/density transfer before claiming the full
  parent objective.

Signed: `herdr optical-review pane w3:p3`, independent non-author
materialization reviewer, 2026-09-05.
