# Validation

The existing `hj2.*` symbolic receipt remains immutable evidence for the exact
algebraic identities checked by `verify_hj2.py`; the verifier was not rerun for
this wording/domain correction.

## Established checks

The receipt continues to corroborate the exact source rescaling and first and
second local Cao cells, the Piola determinant/Jacobian cancellation, the first
and second point-difference Biot--Savart shape formulas, the toroidal-distance
and phase expansion, and the algebraic inverse/Leray differentiation signs.
It also checks the scalar resonance and remainder scale arithmetic.

## Analytic boundary exposed by review 0057

The finite symbolic tally does not prove that Cao's auxiliary system (3.36) is
an exact augmented steady row, a C3 carrier branch, a common closed DA graph
range, front/side/overlap Hs symbol seminorms, product decay in two independent
mode indices, or a nonnormal all-sector X_* to D_* graph resolvent. Those
claims are reclassified as blocked or conditional in `derivation.md` and
`result.yaml`. In particular, the scalar matched integral and one-Lambda
commutator checks are necessary identities but not the operator estimates
formerly inferred from them.

No formula predicate or implementation changed, so rerunning the unchanged
oracle would not test the corrected functional-analytic statements.
