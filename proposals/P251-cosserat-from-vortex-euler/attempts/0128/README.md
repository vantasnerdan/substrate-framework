# 0128 — localized actual 3Ω correction with measured angle and spin

Owner `/root/construction_review`, this new directory only. Accepted
v0.175.0 conditional claims and frozen 0123 remain unchanged. Parent
owns the slow-sideband/cancellation algebra and registers this attempt.

Frozen positive target: construct smooth finite-norm Kelvin/Euler data
whose finite-time physical core-angle and actual material-spin rows
approximate the parent's explicit local 3Ω polynomial mode. The
correction amplitude is O(K²) at fixed microscopic localization scales;
derive its action self terms AND its symplectic cross terms with the
retained mode. It is a prepared state correction, not an observable
redefinition or an asserted autonomous optical eigenmode.

Candidate A: radial localization of the axial Fourier mode, exploiting
the massive Helmholtz pressure and finite-time off-diagonal propagation.
Candidate B: a compact axial envelope with N R≫1 and R≪a, retaining
the zero-axial-frequency leakage and complete pressure response. The
exact local mode and its vorticity/Kelvin generator are derived first.
No measured frequency, spectral discretization or numerical soft-splitting
gate is used. Selection is exact Euler compatibility, finite action norm,
controlled actual angle/spin error, and a regular macro K→0 preparation.

Requires: exact uniform-rotation linearized Euler/Kelvin evolution;
0123's full material spin and collar-kernel identities; parent's 0125
explicit 3Ω mode and physical sideband. Pass licenses: the localized
correction and its quantified physical rows/action cross terms. It does
not license equality at every time to the noncompact polynomial mode,
global torus transfer, or autonomous continuum constitutive closure.

Oracle: analytic residual and kernel bounds plus exact SymPy checks of
the local Euler mode, compact divergence/Kelvin completion, physical
observables, and the correction's quadratic action expansion. No
empirical comparator or production numerics.

## Failure-derived route and frozen result

The exact local single-N mode revealed a cleaner representation:
compact axial primitives with zero global first moment. The original
odd fundamental's inverse-square axial correction is not L2 when its
global first moment is nonzero. A finite weak outside-tag return cancels
that moment without changing the core slope or tag moment; its complete
norm and pressure contribution remain in the retuned matching equation.

`localized-three-omega.md` establishes the resulting arbitrary-profile
3Ω Euler/Lin field, its fully compact radial-potential preparation,
and finite-time controlled core-angle, material-spin and shape rows.
Two actual profiles give independent angle/spin controls. A third
off-core profile independently prescribes the leading KKS cross with
the actual tapered fundamental. All controls are finite geometry and
profile integrals, fixed before the macro wave number is taken small.

The O(K²) action cross is evaluated rather than discarded because the
correction's self-energy is O(K4). After KKS-cross matching, the remaining
Hamiltonian cross is an explicit pressure-residual integral with a
collar-distance bound. The physical finite-time error is O(K² ε_loc),
not O(K4) at fixed localization. Reference-time angle/spin matching is
exact for the complete measured response matrix; constant preparation
does not impose exact equality of all future time traces.

`verify.py` passed all 22 exact checks on its first execution, exit 0,
3.580131109 seconds, captured in `first-run.txt`. It independently
derives the full local Euler and Lin equations, compact potential,
pressure-collar source, actual spin and shape, finite-tag independent
controls, KKS cross integrals and macro-order action expansion. Exact
taper rank minor: `-555751/1048576`. No numerical spectral design.
Ruff passes and the new-proof whitespace check has no diagnostics
(no-index diff exit 1 denotes the new file).

Route verdict: established for the compact-profile construction and
its explicitly controlled finite-time scope. Evidence scope:
`EXACT_LOCAL_MODE_PLUS_COMPACT_PREPARATION_CONTROLLED_FINITE_TIME_
OBSERVABLES_AND_ACTION_CROSSES`.
The single-N massive-pressure and generic axial-envelope methods were
considered comparison candidates; their separate localization estimates
were superseded, not claimed executed. Global fixed-curl torus transfer
and an autonomous constitutive closure remain distinct parent work.
Only this attempt directory was changed; frozen 0123 and canonical
modules were not edited. This child boundary is frozen for integration.
