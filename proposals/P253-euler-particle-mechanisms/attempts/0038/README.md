# P253/0038: independent joined review of the fixed-carrier Kelvin return and Euler lift

This README preregisters one coherent fixed-claim, non-author review of
attempts `0032` and `0039`. Both units concern the same smooth compactly
supported Gavrilov steady Euler carrier: `0032` derives a source-controlled
returned Kelvin polarization cocycle on a derivative-resonant pressure torus,
and `0039` claims to lift its hyperbolic return to the exact whole-space
linearized Euler evolution on finite-energy dynamically accessible
perturbations. The review will decide the joined implication only after
independently checking both links.

The strongest possible result in this review is an exact **linear**
same-carrier statement: the Baldi/Gavrilov recurrence gives the asserted
hyperbolic `SL(2,R)` Kelvin return, and for every finite circuit count its gain
is realized in the operator or essential norm of the actual whole-space Euler
semigroup on the stated accessible quotient. It is not a nonlinear instability
theorem, a single-particle orbit theorem, or a particle/quantum result.

## Independence, ownership, and activation boundary

`particle-balance-review` did not author or implement `0032` or `0039` and
owns only append-only review artifacts under `attempts/0038`. Root and the
original carrier worker retain ownership of the reviewed attempts, central
proposal state, source/API/tests, validation, memory, and commits. The reviewer
will not edit `0032`, `0039`, any canonical module or test, any central file,
generated documentation, memory, or Git history.

Before this freeze the reviewer opened only `0032/README.md`, SHA-256
`8772d196fa8f86466dcf418f79177ff58ca189f67c2382d3a594ce64aff1475d`,
and `0039/README.md`, SHA-256
`be8db67ef25fc542ef5721f87958386041ad1f70c8fe7ca95157076c54943fe4`.
The reviewer also inspected repository/Git metadata, calculated hashes without
displaying target bodies, and read the three primary sources pinned below.
No `0032` or `0039` derivation, access inventory, result, validation, verifier,
command, stdout, or stderr body was opened.

A required repository-memory search exposed a high-level coordinator summary
asserting a corrected physical-frame coefficient
`tr M=2+22*pi^2 I+O(I^(3/2))`, positive discriminant, and an active global
lift. The user's review prompt itself also names the hyperbolic return. These
are non-authoritative target statements, not evidence and not a preferred
verdict; every coefficient and bridge remains to be derived from the pinned
sources and bodies after activation.

No scientific target body may be opened until root:

1. registers `0038` centrally with this joined scope and independent owner;
2. records the repository-schema invocation, stdout, stderr, and exit under
   `attempts/0038`; and
3. makes `attempts/0038/activation-schema.exit` contain exactly `0`.

After activation, the reviewer will pin this README and the activated target
boundary, perform one substantive review pass, and request at most one bounded
evidence-driven correction pass. The accepted base is release `v0.183.0`; the
preactivation Git head is
`1d378f9aebd10a8ff0dd91f26bf15d371f3e243e`. Shared-worktree changes belong to
their existing owners and will be preserved.

## Frozen target inventory

The exact reviewed artifacts are:

- `0032/derivation.md`, `result.yaml`, and `validation.md`;
- `0032/verify_c2_recurrence.py`,
  `verify_c2_source_recurrence.py`, `verify_pullback_coefficient.py`, and
  `verify_trace_metric.py`, with the three nonactivation command/exit/stdout
  receipt families;
- `0039/access-inventory.md`, `derivation.md`, `result.yaml`, and
  `validation.md`; and
- `0039/verify_accessible_semigroup_lift.py` with its nonactivation
  command/exit/stdout receipt family.

No canonical `src/` module or `tests/` file is named by the two frozen target
READMEs or present as a new same-unit artifact at preregistration, so none is
presumed. If central activation explicitly identifies an already-materialized
source/test artifact as part of this exact frozen transaction, its path and
pre-review hash must be appended before it is opened; otherwise executable
attempt scripts are corroborating oracles only.

The following SHA-256 values were calculated without displaying the bodies.
Activation must preserve them or explicitly identify the replacement frozen
boundary before review:

| Artifact | Preactivation SHA-256 |
|---|---|
| `0032/derivation.md` | `a1c252d729716656d758f8449711463a671b41ba17c3eaf4d15351a180d86a73` |
| `0032/result.yaml` | `a25166e2ff0009313d3ba36ca8d09ac105b86c5e0a806891ac56a3a650865c8f` |
| `0032/validation.md` | `717b03be39c581b7b23b440ed4c9cb5272cfdef76f4c5f7d167e2784d58566a6` |
| `0032/verify_c2_recurrence.py` | `4dae3367b62dbaf7641973ac6edc0b56513531deaadaaa3bb9a6d5c65b4d4d2b` |
| `0032/verify_c2_source_recurrence.py` | `194fa4b83f08689bd77dae3cc0da6f0a26b12ffe72dd570e716367754a5d95c8` |
| `0032/verify_pullback_coefficient.py` | `d2615996405ca5c75c72b8a5498e26dcbf25942729cd0bd6465a2fe313782ab3` |
| `0032/verify_trace_metric.py` | `585c66b7f1eb4e71e56f80bba3ea2fa14866d789c80753ecf5840a04d4ef894f` |
| `0032/first-run.command.txt` | `7dc40971db21c39b24c55562649924df6fc0b1c2b45b0b4ede9aa2a2dae8f867` |
| `0032/first-run.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0032/first-run.stdout.txt` | `075804141bc99d67ff59dc2396a4b0f60280b2428511d5dea1aecb48fe8a79ed` |
| `0032/second-run.command.txt` | `b35ff4f3e56fc2f67412eb5d8f68378572c6111c25a4ac03e19494e49e1e13b2` |
| `0032/second-run.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0032/second-run.stdout.txt` | `d4057a8085d19f17931bfbccbb7efaf118ac6fadbf067dac09279e59f88dd4be` |
| `0032/c2-first-success.command.txt` | `f3829407139d2e944e2a7e936478fb067812e7b7b01c82999cf2baf536c76d9c` |
| `0032/c2-first-success.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0032/c2-first-success.stdout.txt` | `6fb56ad2baf1f7ebdbf2f67e88a000e80bca753c11c0a7859edfd0949c71e878` |
| `0039/access-inventory.md` | `72595f4e2fca7374509da106ecbdce21dcd4ddc26c2cfde6d927413f4649133f` |
| `0039/derivation.md` | `0a47bea9732796586208d31f38490eee16e2d72c8a4a315ffee4e7a7bdea4fbb` |
| `0039/result.yaml` | `9b2e33dfaaac5bf5d85f9b96f9738069bfbfcbe3ad11b1c2e432cb5055a1ff1b` |
| `0039/validation.md` | `b2d932b7470424fc32a433953b15d7d56c3651e047ede5df6081f0aaca9c557f` |
| `0039/verify_accessible_semigroup_lift.py` | `375f4c94c6098098baabdceb6beb3495471e20087f1c720082e798d9d5f762ba` |
| `0039/first-run.command.txt` | `7c48aa1ff2d0cf4f4cdbb9316d55ffc5a3791a3a533c494edb28637cf75711e4` |
| `0039/first-run.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0039/first-run.stdout.txt` | `4599bc20a997562ec7f1909b339e2488a3380b32ce848d2fb48aa84bf5c874b9` |

Attempt `0037`, all later edits outside this activated boundary, and every
stability or particle mechanism are excluded.

## Frozen source versions and applicability boundary

The primary-source inventory is closed before target bodies are opened:

| Source | Exact version and local evidence | Permitted role |
|---|---|---|
| A. Gavrilov, *A steady Euler flow with compact support* | `arXiv:1810.08020v1`; cached PDF `/tmp/primary-source-cache/P253-0010/1810.08020.pdf`; SHA-256 `fcaca85faa77e3876b11d16718037169fc026112dfd3e4248e03e963a0ebc3c9` | Exact smooth compact steady field, pressure localization/cutoff, support and regularity. |
| P. Baldi, *The dynamics of Gavrilov's compactly supported Euler flows* | `arXiv:2302.02982v1`; cached PDF `/tmp/primary-source-cache/P253-0025/2302.02982v1.pdf`; SHA-256 `4fe82be0db15a7f117017e505e6682b3438ed14ebcf56225bd7e74bc759136f3` | Theorem 1.1; action-angle map (3.52)--(3.68); pressure/action relation; source recurrences (4.29)--(4.61), including `K(I)=I+1065 I^3/1024+O(I^4)` and `Omega_2/Omega_1=sqrt(I)(1+7I/4+O(I^2))`. It does not itself state the reviewed Kelvin trace or Euler-semigroup theorem. |
| R. Shvydkoy and M. Vishik, *On Spectrum of the Linearized 3D Euler Equation*, Dynamics of PDE **1** (2004), 49--63 | Cached 15-page journal PDF `/tmp/primary-source-cache/P253-0039/Shvydkoy-Vishik-2004.pdf`; SHA-256 `9307efc6096565f74c5a099af25c834af09ba93b3e56ff34dd076ba289551592` | Equations (1.3)--(1.5), Theorem 1.1, Proposition 2.1, and the weakly-null pseudodifferential packet construction on `L^2_div(T^3)`. The paper states its theorem on the torus; its remark that free space should present no major difficulty is not a whole-space theorem. It does not establish invariance of the repository's dynamically accessible closure, symmetry quotient descent, or the claimed `R^3` graph-domain result. |

No inaccessible source is treated as an imported theorem. Standard smooth
Euler well-posedness for the fixed steady coefficient, maximal transport plus
bounded-shear semigroup theory, Calderon--Zygmund/Leray facts, elementary
finite-dimensional Floquet theory, and compact-operator/essential-norm facts
may be used only with their hypotheses and exact role stated. Any additional
published Egorov, essential-spectrum, or microlocal theorem appearing in an
activated body must be treated as motivation unless its exact version,
statement, domain, and hypothesis map are pinned in the review before use.

Reviewed attempt `0019` through `0026` may supply its exact fixed-time
whole-space linear Euler realization and normalized accessible-packet result.
Attempt `0025` may supply only its explicitly stated moving-fibre, pressure,
collar, metric, and finite-time quantifier ledger. Neither supplier proves the
new recurrence coefficient, repeated-circuit lower bound, essential norm, or
generator statement by inheritance.

## Unit A: source recurrence, accessibility, and hyperbolic return

1. Reconstruct Baldi's action-angle chart and all small-`I` jets actually used
   from the pinned recurrence. Track the cutoff factor, pressure/action
   reparametrization, square-root parameter, remainder order, and the interval
   on which derivative resonance is selected. A symbolic script must derive,
   rather than encode, every source coefficient it claims to verify.
2. Pull the physical steady velocity gradient and Kelvin system into the
   source chart. In cylindrical orthonormal components, retain the rotating
   basis/connection terms as well as coordinate derivatives. Check the
   cotangent lift, physical metric, density, returned frame, and axial-rotation
   identification; a Euclidean chart matrix is not the physical gradient.
3. Derive the exact dynamically accessible principal perturbation from a
   compact divergence-free generator. Check the factor and sign in

       delta omega_pr = i[(k dot omega)d-(k dot d)omega],

   the transverse constraint, the `mu=omega dot k` term in the amplitude
   equations, and whether the chosen growing polarization lies in the actual
   image. Accessibility may restrict amplitudes and cannot be inferred merely
   from `k dot A=0`.
4. Recompute the one-circuit physical monodromy through the order needed for
   its trace. Audit the claimed coefficients

       tr M = 2 + 22*pi^2 I + O(I^(3/2)),
       det M = 1,
       Delta = (tr M)^2-4 = 88*pi^2 I + O(I^(3/2)).

   The determinant must follow in the returned oriented physical frame, not
   from an assumed scalar model. Verify that the remainder is controlled
   uniformly on the selected sufficiently small derivative-resonant shells so
   the positive discriminant and reciprocal real multipliers genuinely
   follow. Distinguish a full flow period from a meridional return modulo axial
   symmetry and derive the frame rotation used in either case.
5. Check that coprime resonance integers and a shell wholly inside one flat
   cutoff plateau actually exist with the stated parameter quantifiers. A
   small-shell asymptotic family must still refer to one fixed carrier when the
   later semigroup statement does.

Exposing mutations include deleting cylindrical frame rotation, changing the
`1065/1024` or `7/4` source jets, omitting the accessibility factor or `mu`
term, using a chart metric in place of the physical metric, reversing the
cotangent sign, or testing only `det M=1` without the trace remainder.

## Unit B: exact whole-space Euler microlocal lift

6. Identify the actual closed operator on
   `X=L^2_sigma(R^3)`, its maximal transport/graph domain, and the bounded
   shear plus whole-space Leray realization. Confirm that the full pressure
   projection and nonlocal Hodge/Biot--Savart tails are retained. Compact
   carrier vorticity or packet generator support does not imply compact
   velocity support.
7. Prove that the kinetic closure of compact coadjoint tangents is invariant
   under the exact linearized flow, and state the exact velocity/vorticity map.
   Verify graph-domain membership of the high-frequency packets and of their
   propagated images. Energy normalization may coexist with graph norms that
   diverge as frequency increases; it cannot silently become a bounded-graph
   family.
8. Define the Euclidean symmetry subspace and quotient norm. Check that it is
   finite-dimensional and invariant, that the semigroup descends, and that the
   weakly oscillatory packet sequence remains separated from translations and
   rotations. A lower bound before quotienting is not automatically a lower
   bound after quotienting.
9. Audit the whole-space finite-time parametrix/Egorov step directly. Derive
   the principal Leray symbol that gives the full Kelvin amplitude, all
   subprincipal commutators, cutoff and chart errors, transported density,
   collar separation, and exterior tails in the physical `L^2` norm. If the
   Shvydkoy--Vishik theorem is cited, distinguish its torus/full-space-
   perturbation hypotheses from the new `R^3` accessible quotient and prove
   the missing extension rather than importing a remark.
10. Preserve the exact quantifier order. For each integer `j>=1` and
    `t_j=jT_*`, the construction may choose a new tube width and then a new
    frequency, with constants depending on `j`; it need not produce one datum
    or one asymptotic expansion uniform as `j` tends to infinity. Conversely,
    a statement proved only for one fixed `j`, or with a carrier/frequency
    changed with `j`, does not establish the claimed one-semigroup growth.
11. Separate operator norm, essential norm, quotient norm, approximate point
    spectrum, and generator spectrum. An operator-norm lower bound needs one
    normalized packet per `j`; an essential-norm bound needs, for each fixed
    `j`, a weakly-null normalized sequence whose images retain the gain modulo
    every compact operator. A Weyl sequence additionally needs the exact graph
    domain and residual. None of these implications may be substituted for
    another by nomenclature.
12. Check the deduction from the per-circuit estimate to

       limsup_j (jT_*)^-1 log ||S(jT_*)|| >= log(lambda_+)/T_*.

    Track every multiplicative loss and ensure it is subexponential in `j`.
    The allowed `j`-dependent localization and frequency must not introduce a
    loss whose logarithm is linear with a coefficient that changes the claimed
    exponent.

Exposing mutations include dropping the Leray order-zero term, replacing a
Hodge tail by compact support, reversing the transported covector, removing a
cutoff commutator, swapping `for every j exists h_j` with one uniform packet,
using a single sequence element to claim essential norm, or quotienting by a
noninvariant symmetry space.

## Oracle, verdict, and exclusions

The strongest oracle is independent source-jet reconstruction, physical-frame
Kelvin reduction and perturbative monodromy, followed by a direct functional-
analytic audit of the exact Euler semigroup and a complete microlocal remainder
estimate on each finite interval. The attempt scripts and captured zero exits
may corroborate algebra, coefficients, factor sensitivity, and a finite-time
operator implication only after their predicates are inspected. They cannot
replace the frame, source-remainder, accessibility, whole-space pressure,
weak-null, graph-domain, quotient, or quantifier proofs.

No production numerical observable, soft eigenvalue, stability boundary,
force, or energy splitting is registered, so `small-ratio-numerics` does not
bind this review. Exact small-`I` asymptotics are audited through source
remainders and analytic sign control, not floating-point sampling.

The joined fixed-carrier route receives one route verdict: `established` only
if Unit A's exact accessible hyperbolic return and Unit B's actual whole-space
linear Euler lift both hold at compatible quantifiers; `refuted` only when a
claimed link contradicts the source equations or functional analysis; or
`blocked` with the precise missing recurrence, accessibility, parametrix,
tail, quotient, domain, or essential-norm construction. Component verdicts
will be recorded separately so a true hyperbolic ray result survives any
global-lift gap, and a valid finite-time lift is not overstated as essential
spectrum or nonlinear instability.

After activation, write only `0038/review.md` and `0038/verdicts.yaml`, pin all
inputs and evidence roles, preserve the strongest supported theorem, and use
at most one bounded correction check if concrete evidence requires it. Report
the README-ready state and eventual final artifacts to both
`herdr agent prompt w4:p2 "..."` and
`herdr agent prompt w4:p1 "..."`.

Active `0037`, nonlinear Euler instability, orbital breakup, stability,
particle mechanics, electron/neutrino identification, spin, statistics,
quantization, and every other quantum claim remain outside this review.
