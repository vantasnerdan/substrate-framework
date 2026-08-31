---
name: small-ratio-numerics
description: Read BEFORE freezing the design of any numerical verifier whose gates include a soft Hessian eigenvalue, lambda_min or Morse-index claim, a stability-window edge, a force, or an energy difference or splitting that could sit within about three orders of magnitude of the discretization, quadrature, or roundoff floor — including when reproducing prior work or reusing committed machinery, and whenever results move with box size, mesh, basis order, quadrature nodes, or thread count. Field-tested methods for computing with small ratios — soft modes, weak forces, tiny splittings — drawn from Skyrme, Einstein–Skyrme, boson-star, and multi-scale continuum practice — error budgets, zero-mode gauges, eigenpair residuals, observed-order extrapolation, asymptotic matching, and precision ladders.
---

# Small-Ratio Numerics

A synthesis of how communities that live in the small-ratio regime — nuclear Skyrme, baby-Skyrme, magnetic skyrmions, Einstein–Skyrme, boson stars, liquid-crystal defect energetics — handle a problem that defeats naive single-grid minimization. The recurring pattern: a strong short-distance scale sets the core and the bulk of the energy; a much weaker long-distance scale carries the physics you care about; the question (a soft Hessian eigenvalue, a tiny energy difference, a weak force) lives in the gap between them. Once $\varepsilon$ or $\lambda_{\min}/\lambda_2$ drops below ~1e-3, the weak signal is usually smaller than truncation error on any single grid, and the calculation reports spurious instability, wrong force signs, runaway frequencies, or "minima" that move with the box.

The field's response was not more precision — it was splitting the question into pieces that each have a well-conditioned answer. That split, and the techniques around it, is what this skill shares. It is advice, not a work plan; pick the pieces that fit the situation. (Origin: the methods note in issue #155 and the literature it cites.)

## Eligibility before precision

Do not form a Hessian, generalized eigenproblem, fixed-$J$ functional, force,
splitting, or radiation threshold until the mathematical object, ensemble,
admissible representation, and observable have earned their upstream licenses.
The background must satisfy the correct Euler–Lagrange equations in the claimed
function space; record both its residual and a conditioning-aware forward-error
estimate. Build the full second variation by differentiating the complete
functional before restricting zero-valued fields or blocks. Declare the
constrained tangent space, gauges and moduli, kinetic-metric rank and sign,
asymptotic operator, and essential-spectrum threshold.

For a collective clock, first establish a globally defined compact action,
fixed period and generator normalization, action invariance, vacuum fixation,
conserved Noether charge, and nondegenerate Legendre map. Only then may “fix
$J$” mean minimization on a physical fixed-charge manifold. A localized
positive-norm eigenmode establishes `MODE_EXISTS`; it does not establish a
cyclic phase, Noether charge, fixed-$J$ clock, or nonlinear periodic orbit.

## Constrain the obvious moduli, then minimize

An unconstrained energy should not be asked to be strictly convex. Topological charge, translations, rotations, isorotations, and a conserved angular momentum $J$ or Casimir are exactly flat or extremely soft; leaving them free makes the Hessian look marginal even when the physical profile is stable. After the eligibility gate licenses the topology, compact action, conserved and normalized $J$, and Legendre map, fix those quantities and minimize the **full fixed-$J$ Routhian** on that manifold—the frequency $\omega$ is then an output rather than a coordinate that can run away. Critical frequencies (e.g. $\omega \le \min(\mu, 1)$ for isospinning hedgehogs) are diagnosed on the constrained ansatz, then checked against symmetry-breaking perturbations. Virial/Derrick identities make a cheap independent monitor of the minimizer: if the virial residual exceeds the energy difference you care about, that difference is not yet a result.

This is the rigid-rotor / collective-coordinate split only after those licenses
hold. Restoring spin or isospin on a separately optimized static profile is a
kinematic approximation; it does not establish fixed-$J$ dynamics or a
fixed-$J$ stationary solution unless the profile is solved from the complete
Routhian.

## Soft eigenvalues: spectrum, extrapolation, eigenvector

$\lambda_{\min} > 0$ on one mesh in one box is necessary, never sufficient, when it sits three or more orders below $\lambda_2$.

- **Build the second-variation operator and extract a few eigenvalues, not just a sign.** Radial Schrödinger/Sturm–Liouville problems for hedgehogs; sparse Hessians in 3-D. Symmetrize the discrete operator explicitly and record $\|H - H^\top\|$ — discretization breaks symmetry at exactly the 1e-12 scales in play, and a symmetric solver fed a slightly non-symmetric matrix quietly answers a different question. Equilibrate so the condition number reflects the physics, not the units. Shift-invert Arnoldi (ARPACK and descendants) is the default when conditioning is extreme; place the shift slightly *negative* of zero so the near-zero cluster is separated from it, and deflate licensed zero modes rather than asking the solver to resolve a cluster at the origin. For a generalized pencil $Hv=\lambda Kv$, first restrict to a $K$-positive constrained subspace and normalize with $v^TKv=1$. Report $r=Hv-\lambda Kv$ and $\|r\|/(\|H\|\|v\|+|\lambda|\|K\|\|v\|)$. A forward eigenvalue or sign enclosure additionally requires $K$ conditioning and spectral separation, or a controlled transform such as $K^{-1/2}HK^{-1/2}$. If $K$ is singular or indefinite, no physical-frequency verdict is available until the admissible positive subspace is derived. “The solver converged” is not an error bound.
- **The licensed zero modes are an error gauge.** Only analytically established symmetry directions preserved by the boundary conditions and representation qualify. A taper, mask, wall, gauge choice, or restricted chart may lift a nominal zero mode physically or numerically. For licensed modes, their discrete nonzero scale measures error on that mesh, box, and operator; a $\lambda_{\min}$ below it is unresolved.
- **Extrapolate in mesh and in domain — after verifying the order.** $\lambda(h) = \lambda_\infty + A h^p + \cdots$, $\lambda(R) = \lambda_\infty + B e^{-mR}$ or $R^{-q} + \cdots$. Use a crossed $h\times R$ design (and basis/order where relevant): vary $h$ at fixed $R$ and $R$ at fixed $h$, tracking branch identity throughout. Co-varying mesh and box is not a continuum test. Richardson extrapolation is only valid in the asymptotic regime. A free-exponent three-parameter fit requires at least four rungs plus a held-out or competing-model check; three rungs are admissible only when $p$ was independently fixed before opening the values. Zero modes from broken continuous symmetries can converge as slowly as $h^{1/8}$ on naive lattices — precisely the case where assuming design order lies; improved stencils restore it. If $\lambda_{\min}(R)$ collapses as the box grows, the mode was a finite-domain gift.
- **Look at the eigenvector.** A bulk shape mode, a boundary-layer splitting mode, and a mesh-scale oscillation are different objects; only the first says something about the continuum object. A soft direction concentrated near the outer boundary calls for a larger domain or a different boundary condition, not a verdict.
- **Treat the two signs asymmetrically.** One admissible negative Rayleigh
  witness can establish instability of the stated candidate. Random positive
  directions cannot establish positivity or Morse index zero; that requires a
  converged lowest spectrum plus control of the omitted tangent-space complement
  and essential spectrum.
- **Separate bound modes from the scattering continuum.** Any finite box discretizes the continuum; a "spectrum" not compared against the asymptotic linear operator mixes vibrations with box-quantized radiation.
- **Where a proof is possible, prove it.** The $B=1$ Skyrmion's linear stability was reduced to a radial Schrödinger operator with no bound states — the gold standard the numerics approximate.
- **Reporting habit**: always quote $\lambda_{\min}/\lambda_2$ and its continuum limit, alongside the eigenpair residual and the numerical zero-mode scale. A minimum three orders softer than everything else is a different claim from a well-conditioned well — and a legitimate one (see the last section).

## Weak forces: pair moments, don't subtract energies

Self-energy subtraction at large separation ($E(R) - 2E_1$) is the most common route to a wrong sign — and it is one instance of the general floating-point rule: **never compute a small quantity as the difference of two large ones**. Each energy is known to a few digits, the interaction is many orders smaller, and the difference is noise, box artefact, or orientation error. The cure is always a reformulation in which the integrand *is* the small quantity — the mutual interaction density, the perturbation functional, the difference equation for a splitting between near-degenerate configurations — evaluated directly.

For well-separated objects the field's reformulation is **asymptotic matching**: verify the isolated object has a multipole expansion at infinity; identify the leading moments (Skyrmions have no monopole — dipole or higher); the linear interaction of two well-separated objects is a *pairing* of those moments, scaling as $R^{-(M+N+1)}$ for moment orders $M$, $N$. Relative rotation in space/internal space can flip the sign; the pairing is often attractive for $|N-M| \le 2$. Compute the isolated moments to comfortable precision, evaluate the pairing analytically, and — if at all — use a full two-body run only as a sign-and-coefficient check at moderate separation. This transfers unchanged to any linear far field: massive/massless pions, Maxwell, GEM, linearized gravity. The kernel changes; the strategy does not.

Before interpreting that pairing as a force or bound-state mechanism, license
the **two-object ensemble** rather than borrowing the one-object clock. For two
phases or collective coordinates, derive the joint symmetry rank and the full
$2\times2$ inertia/kinetic matrix, show it is positive and nondegenerate on the
physical quotient, distinguish fixed charges from fixed frequencies, and
relax the jointly coupled two-body functional. A cross term computed on two
frozen isolated profiles is a field-pairing coefficient; it is not yet a
fixed-charge force, a nonlinear equilibrium, or a synchronized clock.

Two arithmetic details decide whether the far field is usable. $e^{-mR}$ at the box sizes this regime needs sits in denormal territory or underflows outright — and denormals are flushed to zero under common fast-math settings — so carry tail amplitudes in log space and match there, reserving `expm1`/`log1p` for the near-cancelling combinations. And moment and pairing integrals are integrals of exponentially decaying functions: sampled on the bulk grid they waste the whole strategy. Use a compactified coordinate or double-exponential (tanh-sinh) quadrature for the matching integrals, and enter the quoted truncation error into the error budget below.

An imposed taper, hard support, pinned wall, compact extension, mask, or core
cutoff is a construction input. Independence from a box outside that imposed
support is tautological and cannot establish localization by the unconstrained
field equations. Such a calculation is a kinematic or sensitivity diagnostic
until the equations select the decay and width.

## Separate scales instead of refining one grid

When the weak sector is a correction — gravitational coupling, small pion mass, quartic boost response, an $\varepsilon$ in front of a curvature-squared term — put $\varepsilon$ in the equations *on purpose*: solve the dominant theory to high accuracy and freeze that core; treat the weak sector as perturbation or slow collective-coordinate dynamics of the core's moduli; match an inner nonlinear core to an outer linearized field, where the force, frequency shift, and interaction sign live. The gravitating-soliton version is adiabatic continuation: start from the flat-space minimizer, step the coupling upward with Newton–Raphson at each step, constraint violations and virial identities monitoring accuracy — and parameterize the branch by pseudo-arclength, not by the coupling itself. Natural-parameter stepping stalls exactly at the fold where a stable and an unstable branch coalesce; arclength continuation rounds it, turning an eigenvalue drifting through zero from a failed run into a traversable **bifurcation** point on the branch. Matched asymptotics and EFT reductions are the analytic form of the same split — cheaper than high-precision 3-D, and they make explicit the order at which the weak signal first appears (a response starting at fourth order in amplitude is invisible to a linearized two-defect calculation; that is a statement about the expansion, not the sign).

## Dynamics that tolerate soft directions

Gradient descent on a stiff-plus-soft energy is slow along the soft manifold and noisy across it. Two field-tested replacements:

- **Arrested Newton flow**: evolve a second-order-in-time equation with the static energy as potential, zeroing velocity whenever energy rises. Efficient in very high dimension (GPU Skyrme searches with $10^5$ random starts). Near-degenerate minima joined by shallow barriers are clustered with a metric on observables (energy, size, inertia eigenvalues) rather than collapsed by an aggressive energy tolerance — the "gray zone" is inspected, not auto-deleted.
- **Path methods** (nudged elastic band and relatives): when the question is the barrier between near-minima, or whether a soft direction leads to fission, compute the path. Softness of a local Hessian can mean a long valley rather than an unphysical object — a different statement needing a different diagnostic.

For time-dependent checks (quasi-normal modes, radiation tails, genuine dynamical instability), trust the linear spectrum and a resolved nonlinear evolution together; linear QNMs that disagree with a careful nonlinear run are not believed. High-order differences, small CFL, and — for late-time tails or tiny growth rates — extra floating-point precision are routine there.

## Precision as a targeted tool

Arbitrary precision goes where the *quantity of interest* is tiny — not everywhere. The ladder, cheapest first: (1) non-dimensionalize so the core balance is $O(1)$, small parameter in a coefficient rather than field amplitude; (2) cancellation-free formulation plus compensated summation — long reductions (energy integrals, inner products) at 1e-13 relative move with summation order alone, so route the reductions feeding the quantity of interest through pairwise/Kahan/`math.fsum` accumulation, which removes that noise at source; (3) Richardson extrapolation in $h$ and $R$ with the observed order verified; (4) independent discretizations (finite difference vs spectral vs finite element; 3-D vs symmetry-reduced ODE) — agreement of extrapolated $\lambda_{\min}$ and far-field moments is stronger evidence than extra digits on one code; (5) tight residual control on iterative linear solvers (the Hessian solve, not the energy sum, is usually the precision bottleneck); (6) quad or arbitrary precision on the soft eigenproblem's residual, matching coefficients, or pairing integral; (7) for the final sign call, interval/ball arithmetic (Arb, mpmath's `iv`) on the last-mile computation — the pairing integral, the matched coefficient, the radial eigenproblem — turning "the sign is probably right" into an enclosure. Rung 7 is the numerical endpoint of "where a proof is possible, prove it."

Record storage, accumulation, operator assembly, and solver precision for every
load-bearing replay. Do not infer an end-to-end dtype from prior prose or from
committed machinery: inherited calculations are exactly where an unmeasured
precision floor hides. If reduced precision enters (for example, float32 fields
for a GPU search stage), treat it as a search tool only: the ~6e-8
representation floor caps derived quantities near 1e-7 relative, and stencil
operators can amplify it. Accumulate reductions in float64 or better, promote
the converged configuration, and Newton-polish it before measuring anything.
The low-precision solution is an initial guess, not a measurement surface.

Tolerances are derived, never defaulted. A library default like `np.allclose`'s `atol=1e-8` declares everything in this skill's regime equal to zero, and a purely relative tolerance is meaningless when the target is zero — every absolute tolerance ties to a measured scale. Each claimed quantity carries an **error budget**: itemize background forward error and its propagation into the observable, branch-identification and continuation-history uncertainty, representation and constraint-projection error, observable-definition or normalization error, Hessian/operator-construction error, roundoff ($\kappa \cdot \varepsilon_{\text{mach}}$ through the worst solve), truncation ($A h^p$ with measured $p$), basis/model error, domain and boundary-condition error ($B e^{-mR}$ where applicable), quadrature, solver/eigenpair residual, and evaluator noise (see reproducibility below). Require the signal to exceed their sum by a stated margin (an order of magnitude is customary). A small background residual is backward error, not a forward bound on the result. Quote the budget with the result; a small number without its floor is a mood, not a measurement.

The cheapest independent monitor in the literature: a topological density that does not integrate to its integer at the same accuracy as the claimed energy is a global error bar on everything else.

## Cross-checks that do not share the soft direction

A claim resting on one soft eigenvalue should be confirmed by at least one observable that is not that eigenvalue. The common set: topological charge to near machine precision; far-field moments against the analytic multipole/Yukawa/GEM tail; the virial/Derrick residual; force from the mutual interaction term or the analytic pairing — never from $E(R)-2E_1$; linearized spectrum versus a short nonlinear perturbation (does the soft mode oscillate, radiate, or grow?); a symmetry-reduced 1-D code at extreme resolution against the 3-D code; and a roundoff-jitter test — perturb inputs at machine-$\varepsilon$ scale (or let a stochastic-arithmetic tool such as Verificarlo, Verrou, or CADNA do it systematically) and confirm the sign and leading digits of the quantity of interest hold still. Agreement supports a shallow real minimum only after the background, ensemble, admissible tangent space, full representation, and essential-spectrum licenses are independently established. Correlated checks on the same invalid background or restricted ansatz cannot create those licenses.

## Reading a marginal candidate

A configuration stable only above a critical box size, with $\lambda_{\min}$ three orders below $\lambda_2$, is treated as a **hypothesis about the continuum** — neither a failed idea nor a finished theorem. The questions the field asks, in the order that saves labour:

1. Is the soft eigenvector a bulk deformation, or a boundary/mesh mode?
2. Is $\lambda_{\min}$ above the in-situ error scale — the numerical zero modes and the itemized error budget?
3. Does $\lambda_{\min}$ extrapolate to a strictly positive continuum value as $h \to 0$, $R \to \infty$, with the observed convergence order matching design?
4. Does a virial identity close at the claimed accuracy?
5. Does the far-field expansion exist and match the linearized theory?
6. Can the weak interaction sign be read from that far field without a subtracted energy?
7. With a small parameter $\varepsilon$: does the sign survive an $\varepsilon$-expansion about a controlled core, or only when $\varepsilon$ shares the core's grid?
8. If two near-minima exist, is the barrier computed, or only local curvature?

Effort on 6–8 before 1–4 is clean is wasted; a two-body sign from a single-scale 3-D run is not evidence until 5–7 are clean. That ordering is the main labour-saving device. The labour that looks like "many attempts" in a single-scale minimizer is, here, spent once on a controlled core, a trustworthy tail, and a pairing — after that, small ratios are bookkeeping.

## Reproducibility in this regime

One practice specific to execution environments, learned the hard way here: results at the 1e-13-relative level and below depend on BLAS thread count through reduction order alone (measured on our certified evaluator: bit-stable at fixed thread count, 1.6e-15 relative drift between 1 and 2+ threads). A check that passes as a script and flips under importlib/harness invocation is usually this. The thread count is one member of a family that moves results at the same scale: compiler flags (`-ffast-math` re-associates and contracts to FMA; flush-to-zero kills the denormal tails above), SIMD width differing across CPUs, non-deterministic GPU reductions and atomics, and OpenMP reduction order. Two remedies outrank pinning: compensated summation on the reductions that matter (precision ladder, rung 2) removes the sensitivity at source, and deterministic modes exist where the environment cannot be frozen (MKL's conditional numerical reproducibility, deterministic cuBLAS/framework flags). Pin and record what remains — threads, seeds, library versions, hardware, flags — at module level (so imports see them), measure the evaluator's noise floor once across thread counts and invocation paths, and quote it alongside any tolerance: a claimed accuracy below the floor is a prediction about the runner, not the physics. Iterative chains (continuation ladders, root-finding on stiff residuals, eigenvalue sign calls) amplify the floor; validate their outputs independently of the chain that produced them.

## Sources

The papers behind these practices (methods, not exhaustive bibliography):

- Gudnason & Halcrow, *A smörgåsbord of Skyrmions*, [arXiv:2202.01792](https://arxiv.org/abs/2202.01792) — arrested Newton flow, gray-zone near-degeneracies. Lightly bound Skyrme / NEB: [arXiv:2305.18126](https://arxiv.org/abs/2305.18126).
- Vibrational modes of Skyrmions, *Phys. Rev. D* **98**, 125010 (2018), [doi:10.1103/PhysRevD.98.125010](https://doi.org/10.1103/PhysRevD.98.125010); Creek, Donninger, Schlag & Snelson, *Linear stability of the Skyrmion*, [arXiv:1603.03662](https://arxiv.org/abs/1603.03662); QNM/Roper-like vibrations: [arXiv:1710.00837](https://arxiv.org/abs/1710.00837).
- Manton, Schroers & Singer, *The interaction energy of well-separated Skyrme solitons*, [arXiv:hep-th/0212075](https://arxiv.org/abs/hep-th/0212075); baby-Skyrme dipole picture: [arXiv:2101.07552](https://arxiv.org/abs/2101.07552) and the Piette–Zakrzewski–Manton line.
- Gravitating/multi-scale continuation: *Phys. Rev. D* **109**, 045002 (2024), [doi:10.1103/PhysRevD.109.045002](https://doi.org/10.1103/PhysRevD.109.045002); eigenvalue coalescence at critical coupling: Bratek, [arXiv:math-ph/0505043](https://arxiv.org/abs/math-ph/0505043).
- Isospinning critical frequencies: [arXiv:1309.3907](https://arxiv.org/abs/1309.3907).
- Continuation through folds: Allgower & Georg, *Numerical Continuation Methods*; validated/ball arithmetic: Johansson, *Arb*, [arXiv:1611.02831](https://arxiv.org/abs/1611.02831); stochastic arithmetic: Verificarlo, [arXiv:1509.01347](https://arxiv.org/abs/1509.01347); compensated summation: Higham, *Accuracy and Stability of Numerical Algorithms*, ch. 4.
