# Validation and claim boundary

The analytic derivation supplies the circulation convention, physical clock,
sharp shell maximum, vector Bessel orders, Debye contour bound, and both
radiation ceilings.  The exact helper and focused tests are algebraic
regressions for those formulas.  They do not establish a Cao high-`J`
eigenpair, current trace regularity, KKS normalization, outgoing limiting
absorption, radiated power, or gate time.

The first focused run exposed three SymPy representation-only comparisons:
equivalent square-root factorizations, an unencoded `eta<1` assumption, and
the unevaluated identity `acosh(cosh(alpha))=alpha` for positive `alpha`.
The formulas were unchanged; the repaired assertions use algebraic
simplification and an exact interior rational margin.  That first failure is
preserved separately.

The final public helpers reject decidable luminal/superluminal carrier speeds,
an impossible Debye margin, the unused zero-margin boundary, and a negative
support enlargement.  `focused-v3.*` reports seven passing tests, empty
stderr, and exit `0`.  The six-predicate exact verifier independently checks
the clock, sharp shell maximum, vector orders, both ratios, and Debye rate.

A separate repository-wide run reported 2784 passing tests after this work was
visible in the shared filesystem, but it was launched before the 0087
module/tests were created.  Its launch snapshot therefore makes it a
pre-0087 checkpoint receipt, not validation evidence for this attempt.  The
focused and exact receipts above are the applicable executable evidence.

The bounded 0089 agreement repair changed only the public docstring, focused
verifier print label, and focused test name for
`high_index_limiting_bessel_ratio`: they now call it a finite-`J` leading
predictor at the supplied carrier speed.  The formula, API name, assertions,
and scientific predicates are unchanged, so the existing receipts are
preserved and those checks were deliberately not rerun.  Their historical
stdout retains its original label.

No production numerics are used.  The fixed-`J` route earns a finite-ceiling
theorem.  The fixed-`ell` route earns the finite-`J` predictor `q_0(J)`; an
exact Debye margin requires the open same-carrier supplier to provide a
uniform `C_asym/J` remainder and an actual integer above `J_asym` but below
the speed ceilings.  Shell-trace-to-gate loss remains open on the named
flux/LAP and normalization constructions.
