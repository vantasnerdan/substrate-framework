# 0056 — joint coefficient and continuum-operator reconciliation

Main owns this continuation of N4/N5. The same-field shear, spin and gradient
constructions now make the final coefficient map executable; their material
assembly remains a separate obligation (0057). Base release and objective
are unchanged from 0035. Inputs are the complete same-orbit Schur matrices,
the simultaneous-rotation ensemble and the explicit low-gradient observable
map, not the retired free-director tension formula or an added rigid rotor.

Re-reading the actual current N4/N5 consumers exposed a concrete remaining
defect: N4 has div m=(c_s+c_a) lap Phi+(2c_tr-c_a+c_s) grad div Phi,
but N5 subtracts the two longitudinal coefficients instead of adding them.
The longitudinal spin coefficient is 2(c_s+c_tr), not 2(c_a-c_tr).
This is an implementation/convention failure, not a refutation of the
microscopic route. Exact longitudinal projection is the exposing oracle.

Extract the reusable quadratic energy and full Fourier operator into an
unpromoted canonical module; independently derive the latter by variation
in tests. Update mutable N4/N5 consumers to call it and use the actual
Euler-orbit inertia/locking map. Preserve historical attempt outputs. The
old numerical refinement remains an operator regression with its declared
parameters, not evidence that the microscopic assembly is complete.

Frozen tests: full asymmetric strain/rotation/curvature differentiation,
both curl helicities and longitudinal spin, rigid-frame covariance, complete
mass/potential map, and the singular structure-free branch. No empirical
comparator is opened. Exact algebra fixes the sign and coefficients; no
small-ratio numerical inference or tolerance is used for this repair.

Result: established in the declared conditional scope. The initial new test
compared expanded and factored SymPy expressions structurally, falsely
rejecting the correct longitudinal identity; the preserved failure was
repaired by expanding their difference. Four API tests then pass. Mutable
CST004 and CST005 pass 8 and 17 checks respectively, with DOP853 refinement
retained and the old tension/rigid-inertia substitutions replaced by the
same-orbit K_Psi/J_Psi coefficient map. Independent review confirms the
equations and requested only a prose distinction between excluded
longitudinal displacement and retained longitudinal spin; that is corrected.
No accepted claim or release is changed by this repair.
