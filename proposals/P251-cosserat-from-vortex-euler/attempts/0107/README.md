# 0107 — importable paired common-mean reaction reduction

Owner `/root/construction_review`; bounded additive API extraction under
P251. Freeze0102's exact shared macro momentum and independent opposite-
circulation reactions, with inputs rho,Hq,N,P,B,C and t=helicity*k/2.
Derive the physical (U,Phi) mass and its spatial-order-two jet, the reduced
locking stiffness and the physical coupling invariant by the COMPLETE
quadratic Schur reduction. The C=j and C=0 cases are outcomes, not separate
hard-coded formulas. Retain the explicit reaction-domain condition
P-B²t²/rho>0 and the paired even reaction N before varying common V.

Conditional unpromoted API only; existing APIs remain unchanged. Current
physics governance/oracle references were read in0100 and reused here.
Read the complete existing euler_orbit module and tests before editing.
rg identifies the test file and existing P251 consumers; the new symbol
has no consumers yet. Replay the full directly affected test_euler_orbit
file because this is an additive edit to that module. Other new Dirac
modules and the unexpected untracked uv.lock are outside this write scope
and remain untouched. No uv invocation, package export, claim or release edit.

Exact oracle: independent direct variation of V,s_plus,s_minus, full
Schur and physical coordinate congruence, both curl helicities, explicit
nonzero N, matched/unmatched moments, and finite-domain/sign mutations.
No numerical approximation, empirical comparator, or all-k Euler closure.

## Completed extraction and receipt

Added `paired_euler_mean_reduction` and immutable
`PairedEulerMeanReduction` to the existing module, without changing any
prior definition. The result includes the full rational physical mass,
its spatial-order-two jet, the locking matrix, kappa, j, leading coupling
invariant and the actual reaction-domain margin. The derivation builds
the exact3x3 reaction Hessian of COMMON V and paired p,r, then computes
the complete Schur complement. N enters the even reaction row and is
eliminated before the physical field congruence; no closed-form answer
is substituted as an input.

The exposing oracle independently differentiates the original action in
V,s_plus,s_minus with nonzero N and both curl helicities. It compares
the complete rational action, not merely a copied jet, and derives the
jet by differentiation. Matched C=j gives the exact negative macro
gradient mass and nonzero physical coupling; C=0 retains its extra
spin-gradient inertia and gives zero leading transfer. The test which
incorrectly ties s_plus=s_minus loses the zero-wave-number inertia.
Finite reaction-domain violations and malformed data are exposed.

First execution: full affected `tests/test_euler_orbit.py`,15 passed in
4.24s, exit0; original captured output in `stdout.txt`. Ruff and scoped
diff checks pass. Explicit rg after implementation finds the new API only
in its definition and appended tests; its attempt records provenance.
No indirect old consumer can call the previously nonexistent symbol,
and the complete pre-existing test file also passes unchanged. No uv or
lock-file operation occurred and the unrelated untracked uv.lock remains
preserved. No sibling0106 file, claim registry or package export changed.

SHA256 at this API boundary:

- src/substrate_framework/euler_orbit.py:
  6cc71e6faae74d22eedaa26151ecf9937cd27546382616a42872a3368c3b997a
- tests/test_euler_orbit.py:
  2a07be0a192f094d9b0bf37fc13216be2ffa9449ab077a33807d591f5cbcf193

route_verdict: established as stated.
evidence_scope: exact conditional paired common-mean/reaction algebra;
no microscopic field existence, physical mean construction, unrestricted
trajectory closure or claim promotion inferred from importability.
