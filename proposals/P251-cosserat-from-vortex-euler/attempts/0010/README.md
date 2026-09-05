# Attempt 0010 — N2: factor-2 Coriolis defect cascade, corrected m=2 dispersion

## Route

Fresh re-derivation of the full finite-k linear system from exact identities
(divergence-reground directive), SymPy carrying every factor, cross-checked
by three independent routes: (i) 3D polar-component elimination with exact
basis-rotation terms, (ii) pure-2D contour-dynamics linearization (Biot-Savart
sheet perturbation, no 3D momentum equations, no Bernoulli), (iii) numeric
root of the exact transcendental dispersion (mpmath, exact I2/K2).

## Defect cascade (all prior-record defects; mechanisms named)

1. **Attempts 0006-0008 momentum solution (single-Omega form) DEFECTIVE.**
   The linearized Euler r-hat/theta-hat equations about solid rotation carry
   the Coriolis parameter 2*Om, not Om: (v0.grad)v' = Om d/dth v' contributes
   the basis-rotation terms (-Om v_th, +Om v_r) AND (v'.grad)v0 = Om zhat x v'
   contributes them again. Correct equations:

     -I wt v_r  - 2 Om v_th = -P'/(rho)
     -I wt v_th + 2 Om v_r  = -I m P/(r rho),     D_tile = 4 Om^2 - wt^2

   0006's D = wt^2 - Om^2 dropped one basis term. (Regression: the corrected
   SymPy solution reproduces 0006's recorded forms under Om -> Om/2 exactly,
   pinpointing the dropped factor.)

2. **Interior Poincare coefficient.** Incompressibility elimination on the
   corrected system gives  r^2 P'' + r P' + (lam^2 r^2 - m^2) P = 0  with

     lam^2 = k^2 (4 Om^2 - wt^2)/(w wt).

   On the physical branch (wt -> -Om): mu^2 = -lam^2 = 3 k^2 (log-free power
   series; modified Bessel I2 inside). Supersedes 0008's lambda = Om k/wt
   (factor defect, same missing 2) and the reground-session candidate
   4 Om^2 k^2/wt^2 (that route treated div(D_t v') = D_t(div v') in
   curvilinear coordinates — invalid).

3. **Dispersion.** Corrected k->0 result (three routes agree):

     ** omega in { Omega (m-1) }  physical;  omega = m Omega branch REFUTED **

   - m=1 translation neutrality preserved (omega = 0). PASS.
   - The 0006 "gauge branch" omega = m Om (kinematically advected sheet,
     zero perturbation fields) VIOLATES the frozen-vorticity tangential
     condition v_th,out - v_th,in = +2 Om eta (residual exactly 2 Om eta):
     an advected sheet deformation with no velocity response drags the
     vorticity step inconsistently. The branch is an artifact of 0006's
     incomplete BC set (tangential condition omitted).
   - Kirchhoff probe, corrected reading: pattern speed omega/m =
     Omega(m-1)/m; m=2 gives Omega/2 = zeta/4 — the exact epsilon->0
     Kirchhoff elliptical-patch rotation rate. 0006's Omega/4 FAILED this
     limit; 0006's probe claim ("pattern speed Omega is the epsilon->0
     Kirchhoff rate") was incorrect.

4. **Finite-k dispersion (m=2).** Exact transcendental (I2 inside, K2
   outside, corrected Bernoulli p_out = I rho wt(r) Phi with the advective
   Doppler wt(r) = w - m Om a^2/r^2, kinematic pair, tangential condition).
   Expansion (symbolic series, log-parts identically zero on-branch, both
   residual collapses proportional to (Om + 12 delta); numeric root confirms):

     ** omega = Omega (m-1) - Om (k a)^2 / 12 + O((k a)^4) **

   Numeric: shift/k^2 = -0.083228 -> -1/12 = -0.083333 at k = 0.02..0.2;
   F(root) ~ 4e-29 (true zero of both residuals). No log anomaly at O(k^2)
   (K2 logs cancel exactly on-branch — verifier check).

## Consequences up-ladder

- **0008's soft-sign (Krein) finding RECLASSIFIED**: the negative fixed-eta
  static energy and the block on C_tw were computed on the defective
  single-Omega fields. The corrected system is a standard stable surface
  mode (real frequency, positive stiffness expected). The canonical-energy
  receipt obligation stays open until the corrected static/wave-energy
  recomputation lands (attempt 0011); if that lands positive-definite, the
  Krein receipt requirement is discharged as defective-input.
- **0007 outside piece and 0008 inside piece**: machinery survives, numbers
  recompute on the corrected branch (0011).
- **0009 bend channel**: K-identity and exterior energy unaffected (pure
  exterior math); c_inner matching must use corrected interior dynamics.
- Verifier debt: verify_cst002 must include the single-Omega mutation — the
  exact defect that survived 0006-0008's probe batteries — as a must-FAIL
  mutation, plus the tangential-condition check (gauge-branch killer).

## Route verdicts

- route_verdict: established (corrected dispersion, three-route);
  evidence_scope: SYMBOLIC_EXACT with numeric corroboration.
- Prior-route verdicts updated: 0006 dispersion route refuted-with-mechanism
  (missing basis-rotation term); 0008 C_tw-block route refuted-with-mechanism
  (defective input fields).
- Next route: corrected static-configuration energy and canonical mode
  energy on the corrected branch -> C_tw, M_eff extraction (0011), then
  frozen verifier verify_cst002.py with the defect-derived mutation battery.
