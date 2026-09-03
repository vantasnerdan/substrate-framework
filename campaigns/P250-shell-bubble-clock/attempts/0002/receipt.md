# G1 Analytic Specification Receipt — sigma_0 and the wall profile

Status: receipt frozen BEFORE any production numerical output (owner
directive: numerics follow the sorted concept; S1 numerical_representation
requires "analytic specification receipt, crossed resolution refinement,
itemized error budget").

## Typed chain

- mathematical object: the heteroclinic stationary solution
  u(x) = (m, c, b, f)(x) of the rotating-frame wall EL system
  2K u'' = dV_w/du, K = diag(1/4,1/2,1/2,1/2), at the certified deep-branch
  crossing omega_*^2, connecting vacuum A = (1,0,0,0) (x -> -infinity) to
  the deep bulk B = (0, c*, b*, f*) (x -> +infinity), on the aligned
  real-psi wall slice of accepted C-M5W-001.
- symmetry/conservation license: mechanical first integral T - V_w = const,
  exact (module identity `first_integral_conservation_identity` == True);
  both ends are zeros of V_w with T = 0, so const = 0 on the connecting
  orbit: T = V_w pointwise on the kink.
- exact variational functional: sigma_0 = inf over H1 paths A->B of
  int sqrt(2 V_w) ds_g (C-M5W-001); on the minimizer,
  sigma_0 = int (T + V_w) dx = 2 int V_w dx = 2 int T dx.
- admissible function space and representation: u in C^2(R, R^4),
  exponentially attracted to A and B; discretized on [-L, L] with either
  Dirichlet-at-tail or linearized-tail (exponential-fit) boundary
  conditions; collocation (scipy solve_bvp).
- analytic scale and asymptotic structure (computed exactly this attempt,
  before any BVP run): Hessians in the (m,c,b,f) coordinates,
    H_A = diag(5, 10, 22 - 4 w2, 6 - w2) at w2 = omega_*^2
        = diag(5, 10, 15.3442, 4.3361) (positive definite),
    H_B positive definite with eigenvalues (of (2K)^{-1} H_B)
        6.3889, 8.4046, 10.7356, 49.1637.
  Tail decay matrices D = (2K)^{-1} H; masses (smallest sqrt eigenvalue):
    A end 2.08231945194, B end 2.52762966155; global m_min = 2.0823.
  V_w(B) at the 55-digit certified point = -3.2e-50 (numerical zero).
  Nondegeneracy license inherited from attempt 0001 (Krawczyk uniqueness in
  the 1e-12 rational box; Gershgorin margins 8.1699/4.87896/13.1268).
- observable: sigma_0 (slice-sector kink tension) and the four profile
  functions u(x); monotonicity branch identity.
- irreducible numerical remainder: the VALUE of sigma_0. Attempt 0001
  refuted closed-form routes (lex-Groebner timeout; degree-108 Sturm
  prohibitive; PSLQ degree-5 falsified at 90 dps), so no exact value is
  available; the profile is the minimizer of the exact functional.
- numerical approximation: solve_bvp collocation on a crossed (nodes x L)
  ladder at fixed omega_*^2; three independent tension routes on a refined
  quadrature grid with math.fsum; exponential-fit vs Dirichlet boundary
  treatment as independent domain-error probes; multi-start branch checks.
- permitted verdict: SIGMA0_SLICE_VALUE — numeric evidence for the
  slice-sector kink tension at omega_*^2 and the constructed profile, with
  an itemized error budget. Does NOT license: off-slice minimality (full
  space), bag physics (G2), stability (G4), or any claim promotion by
  itself (claim candidate C-M5W-006 formalizes after results exist).

## Frozen numerical design freedoms (frozen before outputs)

1. Frequency and bulk data: omega_*^2 and (c*, b*, f*) at the 55-digit
   certified values from attempts/0001/maxwell_point.json. Their certified
   uncertainty (1e-12 box) is a measured budget item (b4), not a design
   freedom.
2. Production domain route: L-continuation through L = 2, 3, 4.5, 6, 8, 10,
   12 with tol 1e-6 .. 1e-10 (see amendment); warm-chained node ladder
   n in {401, 801, 1601, 3201} at L=10 and L=12 for the h-crossing.
   (amended from the original direct-solve ladder — see below.)
3. Boundary treatments: (i) Dirichlet at both ends (error O(e^{-m L}));
   (ii) exponential-fit (stable-manifold projection; error O(e^{-2 m L})).
   Agreement of the two at the finest rung is an independent check.
4. Tension routes: R1 = int (T + V_w) dx, R2 = 2 int V_w dx,
   R3 = 2 int T dx, each with math.fsum; their spread is the primary
   consistency signal.
5. Initial guesses: the L=2 rung starts from the converged field-space
   string; direct long-box guesses (tanh-arc, linear ramp, smoothed step)
   are basin probes only.
6. Environment: single-thread pinned at process start (OMP/MKL/OPENBLAS=1)
   before numpy import; recorded versions; deterministic evaluation.
7. Observed order: Richardson fit on sigma_0(n) with order verified from
   consecutive differences before extrapolation is quoted.
8. Error budget items (each quoted with the result):
   b1 collocation residual (solve_bvp residual RMS, tol 1e-11);
   b2 truncation A h^p (observed p);
   b3 domain B e^{-2 m_min L} (fitted from the L ladder, exponential-fit
      BC) with the Dirichlet run as cross-probe;
   b4 certified-box drift: re-solve with bulk B endpoints perturbed by the
      1e-12 certified box (and omega_*^2 endpoints) -> sigma_0 drift;
   b5 quadrature: refined-grid trapezoid vs spline quadrature difference;
   b6 evaluator noise: permuted-grid fsum identities (exact) and an
      input-epsilon jitter of the integrand evaluation;
   b7 branch identity: BC agreement at L_FIX plus monotone branch identity
      across continuation rungs;
   b8 first-integral drift: max_x |T - V_w| along the converged profile.
   Acceptance: total budget <= 1e-3 * sigma_0 with every item reported.

## Production-route amendment (implementation repair, recorded pre-output)
Landscape check (executed before any solver run): V_w >= 0 on the sampled
corridor (min 0 exactly at the bulk states; interior minimum 0.015 at
m = 0.1), so no forbidden region blocks the kink. Direct long-box
collocation solves fall into multi-kink trains (Newton basin structure —
documented as basin probes), and both hand-rolled string and NEB normal
flows proved unstable here. The production route is therefore
L-CONTINUATION of the Dirichlet kink, warm-started per rung (short boxes
admit only the 1-kink), with warm-chained h-ladders at L=10 and L=12, the
exponential-fit BC as the independent boundary-treatment cross-check
(b7), certified-box drift at production settings (b4), and compensated fsum
tension routes. The converged field-space string (string_guess.json)
supplies the L=2 initial guess and an independent order-of-magnitude
cross-check only (its own sigma carries ~0.5% path-convergence uncertainty
and is NOT part of the budget).

## Explicit non-goals
No off-slice statement; no claim that the slice kink is the full-space
minimizer beyond the C-M5W-001/002 exclusion structure (shear/orientation/
psi_I positive on clock-active branches, so the slice carries the lowest
tension among sectors with these asymptotics — recorded as the reason the
slice value is THE wall tension at this scope); no bag; no stability.
