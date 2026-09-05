# Independent review: nonlinear product-Haar contrast

Reviewer: `/root/smooth_core_review`, distinct from the author of 0058
and its API/consumer changes. Date: 2026-09-05. One bounded scientific/code
review under AGENTS.md and the physics skill. The reviewed object is the
static product-Haar theorem and its exact circle oracle, not a new
simulation or a dynamical angular-current closure theorem.

## Decision and supported result

Established as stated: a conditional product-Haar law for the complete
unresolved orientation state makes its averaged integrable energy
independent of independent local coherent frame shifts. Thus its static
coherent torque and stiffness vanish. This does not imply that its average
fluctuation energy vanishes. The exact independent-phase example retains
strictly positive energy, while the correlated-uniform-marginal example
retains nonzero relative-angle stiffness.

No load-bearing correction is requested. The full-state and conditional
correlation premises are substantive statistical assumptions, not
consequences of one-point isotropy, circulation reversal or topology.
The no-retained-coherent-angular-current condition remains separate;
the static theorem does not establish a zero dynamical susceptibility
or absence of rotational probability current.

## Evidence and general Haar proof

I read the README, complete `haar-contrast.md`, the new
`uniform_phase_average` helper and test, the affected CST006 function
boundaries and call site, and the saved API and consumer receipts. The
API receipt passes five tests; the repaired targeted consumer receipt
passes the unchanged signed-response check and three new exact checks.
The historical function-boundary NameError is preserved separately.
The source repair restores the old function's local variables and keeps
the nonlinear checks in their own function. It changes no scientific
equation or tolerance.

For fixed z and coherent frames Q_a, right multiplication
`R_a -> S_a=R_a Q_a` preserves normalized Haar measure on each SO(3)
factor. The assumed conditional joint measure is their product, so the
change of variables can be made independently at every site. The average
is therefore exactly Ebar(I,z), even for finite Q and coupled interactions
in E. Integrability supplies the finite expectation and legitimate product
integration. Constancy itself gives zero derivatives of Ebar; the stated
dominated derivative conditions additionally justify identifying them
with averaged microscopic derivatives when that representation is used.

No expansion around a small angle, first-moment factorization or
cancellation of positive energy terms is needed. Local discrete torque
and neighbor-angle curvature vanish before a continuum approximation.
Dependence on z, such as displacement or strain, remains and is not
assigned a modulus by this argument. A full thermodynamic or infinite-cell
limit would retain its own convergence requirements; the finite/local
static identity does not assume that limit to establish its conclusion.

The words complete orientation state are essential: if an unaveraged
hidden frame remains locked to a visible core while only that core axis
is randomized, the law and frame shift in the theorem may no longer
describe the actual energy. Holding the stated conditional law fixed
under Q, or transforming the remaining conditional variables covariantly,
is the explicit closure that licenses the change of variables. It is
not merely a declaration that the desired coefficient is zero.

## Exact fluctuations and correlation counterexample

For independent uniform theta_1,theta_2, the nonconstant Fourier mode in
`K[1-cos(theta_1-theta_2+q_1-q_2)]` integrates to zero. The constant
mode remains K>0, with zero coherent gradient and Hessian. This example
proves compatibility of positive fluctuation energy with zero conservative
response; positivity is not claimed for every arbitrary integrable energy
without an energy-sign premise.

When theta_2=theta_1 is imposed before integration, both marginal laws
remain uniform but the relative phase is fixed. The integral is then
`K[1-cos(q_1-q_2)]`. Its relative-angle stiffness at zero is K; its
two-coordinate Hessian is K times `[[1,-1],[-1,1]]`, so common rotation
still costs no energy while relative rotation does. This exposes exactly
the false inference from identical uniform one-point data to a product
law. Independence is needed for the general all-interaction guarantee,
not inferred from those marginals.

The SO(2) example is an exposing subgroup oracle. The SO(3) theorem rests
on its own Haar change of variables, not on treating this circle test as
a numerical or algebraic proof of every rotation-group interaction.

## API and consumer scope

The helper performs one normalized integral over [0,2 pi] for each
distinct symbol and leaves unresolved symbolic integrals explicit. This
is exactly the product-circle measure. Its docstring correctly places
periodicity as the additional premise for frame-shift invariance;
the integration routine itself need not infer or enforce periodicity.
Correlations must first be pulled back to their actual shared variables,
as the new test and consumer do. Duplicate phase symbols are rejected,
avoiding an accidental extra integration over one variable.

The new test checks the retained constant energy, zero coherent Hessian,
correlated law and invalid duplicated phases. The consumer calls this
shared implementation and names the general SO(3) proof separately.
Its new checks do not alter the existing Monte Carlo observable or claim
that a finite sample establishes nonlinear Haar invariance. Existing
Monte Carlo evidence remains evidence for its original signed first
moment only.

Finally, a uniform orientation density can carry a nonzero probability
current. Independent orientation averaging also differs from the coherent
circulation-paired ensemble of 0048, whose physical angles are tied and
whose reaction momenta are independently varied. The positive rotor
action and this zero static response therefore concern different declared
ensembles and are mutually consistent. The proof correctly retains the
additional angular-current premise for a dynamical no-spin closure,
without erasing transient memory or viscosity by a static argument.

## Frozen hashes and disposition

- `haar-contrast.md`: `1a6e74ceba87ed265141f0907002c5527059bc44c474ef4003cc6e06a153dff2`
- `src/substrate_framework/micropolar.py`: `7187523783e3dfe30ce062c5b547e5d0524631584dcde13326646d45547a2693`
- `tests/test_micropolar.py`: `d78808625a9af88d6c8151e6dbba4f2f4543273e7e73dfbd1547f16b25664b09`
- `verify_cst006.py`: `d64b4e53eadefb3a519056c87cbe59612ad52542da33129c678ea91890f22a5d`
- `new-consumer-repaired.stdout.txt`: `6ed7845a0bdf0ffdd23edc62ed3166d1596fd8f2fffaab1e052a35e8d515c26c`

Acceptance is recommended for this conditional exact static-response
theorem and its integration helper. It introduces no small-ratio numerical
remainder and no new Monte Carlo claim. The separate dynamical closure
and parent's full microscopic continuum completion are not inferred.
