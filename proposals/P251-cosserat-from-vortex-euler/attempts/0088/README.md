# 0088 — exact Fourier material-stiffness witnesses

Owner `/root`; active P251 continuation. Freeze this representation-changing
candidate before execution: a smooth multiwave Beltrami field on T^3 with
curl eigenvalue N=1,2,3 and a divergence-free real trigonometric material
displacement Xi. Evaluate the ACTUAL material Jacobi stiffness

    K(Xi)=average[Xi.Hess(p0).Xi-|u0.grad Xi|^2],
    p0=-|u0|^2/2, rho=1.

The periodic cell is [0,2pi]^3, with normalized volume measure. An eventual
strictly positive rational Rayleigh witness establishes this one selected
material direction, not Euler spectral stability, a periodic orbit, a
fixed-Kelvin finite-dimensional closure, or an EPS-localized mode. The
single Beltrami plane wave has constant p0 and hence K<=0; multiwave
pressure correlations are the new candidate mechanism. A finite-basis
negative result is not a continuum no-go.

Analytic specification: stationary Euler and pressure are exact finite
Fourier polynomials; derivatives and products are Fourier convolutions,
and integration extracts the zero coefficient without aliasing. The
constraint is k.Xi_k=0 with Xi_-k=conjugate Xi_k. The mass norm is
average |Xi|^2. There is no discretized Euler background, outer truncation,
quadrature or empirical comparator. A finite Fourier test space specifies
admissible witnesses rather than approximating a claimed lowest spectrum.

Search, if used, is explicitly hypothesis generation: float64 symmetric
matrix eigendirections, fixed BLAS threads, followed by independent exact
rational evaluation of any candidate. No search eigenvalue earns a
scientific sign verdict. The exact positive Rayleigh identity is the
stronger final oracle; perturbing its rational coefficients supplies
sensitivity when appropriate. This removes roundoff, mesh extrapolation,
zero-mode gauge, eigenpair residual and lambda_min/lambda_2 requirements
from the FINAL proposition because no eigenvalue or full stability claim
is made. Those skill prescriptions remain applicable if the claim changes.
The small-ratio-numerics skill was fully read before this design freeze.

If a witness exists, the next construction is its compact localization,
physical core-angle jet, and actual material/Kelvin reduction; 0084 owns
that analytic bridge. No old positive orbit-H cage coefficient transfers.

The first search (saved `first.stdout`) produced no positive candidate in
the nine ABC/carrier/basis combinations. Near-zero printed eigenvalues are
search diagnostics, not exact null or exhaustion claims. Failure-generated
expansion: use curl eigenvalue5 and non-axis-aligned Pythagorean wavevectors
(4,3,0), (4,-3,0), (0,4,3), (0,4,-3), (3,0,4), (-3,0,4), in two-, three-
and six-wave sets with fixed quarter-turn phases. `expanded_search.py`
freezes these geometries before their search output; the witness criteria
and exact-verification requirement are unchanged.

The expanded search also produced no positive candidate in its six selected
spaces. Both commands exited zero; these outputs remain hypothesis-generation
diagnostics with no lowest-spectrum or no-go claim. The continuation changed
the physical background rather than tightening numerical tolerances: the
primary force-free energy-principle literature supplied a two-wave ABC trial
in0089 (one ABC amplitude exactly zero, unlike the three equal amplitudes in
the initial search). Exact pressure and curl evaluations independently give
positive K for its rational divergence-free displacement. Main0090 transfers
that genuine material stiffness to a compact same-EPS-field cage.

Route verdict: the searched finite-basis candidate sets are blocked at the
missing positive witness, with evidence scope hypothesis-generation only.
The materially different source-derived candidate is established in0089;
neither verdict closes the parent material-action obligation. No rounded
search eigenvalue is promoted or used to exclude the continuum class.

Post-run implementation correction: the expanded driver originally imported
NumPy before the shared search module set its thread environment variables.
It now imports the environment-setting module first. The historical expanded
output therefore does not assert a measured one-thread run. No replay is
needed: it supplies neither a sign verdict nor the positive witness, and
0090 uses the independently exact0089 result instead.
