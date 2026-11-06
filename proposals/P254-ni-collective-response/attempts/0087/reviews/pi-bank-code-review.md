# A87 Pi/bank/next-Dyson original code review

**Verdict: incorrect; confidence 0.94.** Original pre-correction review. Root accepted both findings for correction; this report does not review those corrections. The companion JSON preserves the complete two-finding array and exact reviewed source/contract byte counts and SHA256 bindings.

## Complete findings

1. **P2 — Authenticate committed pair bytes before trusting resume coverage.** `0087/produce_joint.py:123-128` (confidence 0.98). `PairStore.validate()` verifies finiteness/causality, not previously committed numeric identity. A finite change to a committed future-time greater/lesser entry leaves zero retarded valid; resume skips the contraction and publishes new hashes for altered data. The same store owns cubic old0/new0/family stages. Authenticate durable committed-pair digests before trusting coverage.
2. **P2 — Authenticate saved dependency bindings before resuming contractions.** `0087/produce_joint.py:209-212` (confidence 0.97); same defect in `cubic_update.py:253-255`. Resume pins only the new producer files, then overwrites saved owning-dependency bindings without checking them. Changed A41/A50 or cubic/reference implementations can mix old committed and new remaining pairs under fresh provenance. Authenticate saved dependencies before reconstruction or commits; reject changes.

## Scoped evidence

Read all targeted implementations/interfaces, the actual joint-block oracle and its source, owning A41/A50/A56 formulas, A51 cubic routing, A75/A76 consumers, A77 COW/binding traversal, A61/A63 solver/provider/readout, A86 M2 reader and admitted A85 verification. Read actual G1/M2 identities, retained all64 cubic and paid-nine manifests, source schemas, and a partial EI checkpoint.

Static paths directly sum IF+FI+FF, evaluate full114 Pi once, retain all64-old0+actualG1q0 with literal paid-nine reuse and only952 old-pair restoration, compose Sigma=EE+rawEI+joint correction and Pi=joint+embedded111 cubic, separate actual/original identities, and reuse unchanged public solvers/readout with committed-slot skipping and retained-staging refusal.

The existing oracle records exact identity, zero A41 comparison errors and direct 7i versus subtractive zero in its cancellation fixture. Its scope is only the new rectangular identity, not complete physical output or recovery qualification.

**Legacy-byte boundary:** no analogous joint/cubic bank refusal established statically. A77 `_bindings` skips historical `metadata` in archive/metadata wrappers; actual directly traversed joint archives/boundary and inspected cubic/reference binding lists use integer counts. This is not an observed bank admission.

## Execution and unresolved risks

No reviewed-code imports, tests, builds, Native runs, numerical-input loads or old-proof replay were performed. Only review reports were written.

Root reports joint001 failed **before module import** because its runpy launcher omitted script-directory `sys.path`; **no joint numeric output exists**. EE002/EI001 runtime admission is not complete scientific ownership. Complete joint/cubic output, four-family bank publication, and an actual M2-driven all128/fullGamma Dyson/readout remain unestablished. No accuracy, conservation, convergence, calibrated loss/rate, heat or closure conclusion is made.
