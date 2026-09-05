# 0058 — nonlinear orientation-ergodic contrast, without erasing fluctuations

Main owns this new N6 continuation. Base and positive parent objective remain
0035/P251. The existing first-moment verifier is correct but does not alone
establish vanishing conservative stiffness. Execute a representation change:
average the complete finite-angle microscopic orientation energy over an
independent local Haar ensemble BEFORE taking coarse-angle derivatives.

Frozen theorem: a product Haar law for the complete unresolved orientation
states is invariant under independent local coherent frame shifts. For any
integrable finite-range orientation energy, its average is independent of
those shifts. Thus the coherent conservative torque and couple stiffness
vanish, while the fluctuation energy generally stays positive. Independence
is essential: perfectly correlated uniform marginal phases give a nonzero
relative-angle stiffness. This explicit counterexample is the mutation.

The premise is a measurable statistical closure on states and correlations,
not an assumed zero constitutive modulus. It does not follow from isotropic
one-point means, time reversal, or topology alone. The no-retained-coherent-
angular-current premise remains necessary for the dynamical no-spin closure;
the static Haar theorem does not erase a collective rotational current,
transient memory, or viscosity.

Exact Haar change of variables is the general proof; exact finite Fourier
integrals are the symbolic oracle for sign, measure, independent-cell and
fluctuation bookkeeping. No new simulation or numerical remainder is needed.
Canonical implementation extends the new micropolar.py with the finite-phase
integration helper and tests; existing Monte Carlo signed-response evidence
is reused at its original narrower scope.

Result: general Haar change-of-variables proof established in
haar-contrast.md; the new canonical integration helper and correlated-law
counterexample pass (five total micropolar API tests). The first consumer
edit misplaced the pre-existing signed-response function tail in the new
function, causing a NameError AFTER the three new exact checks passed.
The saved failure is an implementation artifact; restoring the function
boundary gives all four affected checks passing. The unchanged Monte Carlo
design and earlier receipt remain valid. Ruff and diff checks pass.
