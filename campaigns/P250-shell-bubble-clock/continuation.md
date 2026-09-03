# P250 Continuation Brief — read this first after compaction

Campaign: issue #195, proposal P250, branch `research/P250-shell-bubble-clock`,
PR #196 (OPEN, `Fixes #195`). Head at pause: `091c636` + continuation commits.
Claims C-M5W-001..005 are ACCEPTED (release v0.172.0, 257 claims valid).
Review record: `proposals/P250-shell-bubble-clock/reviews/C-M5W-001-005-review.md`
(closed; MC-1..MC-4 all applied and verified by reviewer ClaimReview).

## Owner directives in force (2026-09-02)

1. Numerics are NOT forbidden; they must follow a concrete well-sorted
   concept. The concept is now sorted (accepted claims), so the gap-closing
   numerics below are authorized.
2. Gap-closing work lands on the SAME branch and SAME PR #196 (owner
   decision; the terminal boundary is extended by owner authority).
3. Do not comment on issue #195 or Discussion #186 (owner's standing
   comment) until a PR is approved.
4. PR #196 is NOT to be merged by the authoring agent: distinct merger
   required. The user (repository owner) merges or directs otherwise.

## Exact physics state (everything needed to resume; verified in SymPy)

Canonical action source: `src/substrate_framework/m5_exterior_clock.py`
(C-M5C-001..004). Wall module: `src/substrate_framework/m5_wall_clock.py`.
Tests: `tests/test_m5_wall_clock.py` (15), consumer replay
`tests/test_m5_exterior_clock.py` (9). All green at head.

Aligned real-psi wall slice: S = diag(m, c+b, c-b), psi = f >= 0.
- Kinetic metric K = diag(1/4, 1/2, 1/2, 1/2) in (m, c, b, f).
- Slice potential V0 = V_M5 + 2c^2 + 2b^2 + 6(b-f^2)^2 + W(f),
  V_M5 = -r^2/2 - (m^3 + 2c^3 + 6cb^2) + r^4 + 1/2, r^2 = m^2+2c^2+2b^2,
  W(f) = 3f^2 - 4f^3 + 2f^4.
- Inertia iota = f^2 + 4b^2. Rotating frame: V_w = V0 - (w^2/2) iota.
- First integral on any stationary profile: T - V_w = const,
  T = (m'^2 + 2c'^2 + 2b'^2)/4 + f'^2/2.
- Static slice potential has the unique zero at (1,0,0,0); V_M5(-NN^T) = 2.
- Vacuum: (m,c,b,f) = (1, 0, 0, 0). Deep bulk state: m = 0 exactly with
  (c*, b*, f*) = (0.30280764518677380369, 0.65773437772324925193,
  0.81436149699856776719) — certified unique in a 1e-12 rational box
  (two-step Krawczyk, strict inclusion) with positive-definite fixed-omega
  Hessian (Gershgorin lower bounds 8.1699, 4.87896, 13.1268).
- Deep-branch crossing (Maxwell frequency):
  omega_*^2 = w_of(f*, b*) = 32 f*^2 - 12 f* + 6 - 24 b*
  = 1.6639457000591502988561930002961614442197723...
  Rigorous enclosure (interval, width ~1.2e-43):
  [1.66394570005915029885619300029616144415736458,
   1.66394570005915029885619300029616144428226086].
- Maxwell system (m = 0 branch, differentiate at FIXED omega, then
  substitute w^2 = w_of; substituting first corrupts db — past defect MC-3
  class):
  pA = 8c^3 - 3c^2 + (8b^2+1)c - 3b^2 = 0
  pB = 2(8b^3 + 8bc^2 - 6bc - 2b w^2 + 7b - 6f^2) = 0
  pC = 2 V_w = 0,  w^2 = 32f^2 - 12f + 6 - 24b.
- Exact rational witnesses (GLOBAL): V_w at (0, 31/100, 13/20, 41/50) is
  -13739/18750000 < 0 at w^2 = 5/3 and -33854777/25000000 < 0 at 45/16,
  hence omega_c^2 < 5/3 < 45/16 with omega_c^2 := 2 inf(V0/iota) and
  omega_c^2 <= omega_*^2 definitionally. Identification omega_c^2 =
  omega_*^2 requires the global SOS certificate (G3) — NOT asserted.
- Excluded directions (exact, positive on clock-active branches): shear
  2(8b^2-3b+8c^2-3c+4m^2-3m+1); orientation 48 b f^2 (true rotation angle;
  12 b f^2 in the half-angle parameterization — same physics); psi_I enters
  only quadratically-positive. Orbit-fixed locus: psi = 0, tangent
  isotropic; the whole S1 orbit fixes it; phase jumps across it are exactly
  free (uniform ramp bound L^2 (delta theta)^2 eps / 4 -> 0).

## Gaps to clear, in priority order (owner-selected)

STATUS: G1 DONE (attempt 0002, 2026-09-02): sigma_0 = 0.72929841787(21),
budget total 2.0e-7 PASS, monotone profile in profile_L12_n3201.csv,
production route = warm-started Dirichlet L-continuation. Next: G2.
Solve the planar wall BVP at omega^2 = omega_*^2:
  unknowns u(x) = (m, c, b, f)(x); EL: 2K u'' = dV_w/du;
  BC x -> -infinity: u -> (1,0,0,0); x -> +infinity: u -> (0, c*, b*, f*);
  both ends are zeros of V_w so the first integral reads T = V_w on the
  profile; tension sigma_0 = int(T + V_w) dx = 2 int V_w dx = 2 int T dx.
Design freedoms (already frozen in the proposal's analytic receipt):
  SciPy solve_bvp on a truncated domain [-L, L], tanh-scaled initial guess,
  crossed h-by-L refinement with observed order, second boundary treatment
  (pole/exponential fit at the ends), itemized error budget (background
  residual, truncation Ah^p, domain Be^{-mass L}, quadrature, evaluator
  noise), deterministic single-thread reductions, compensated summation for
  sigma_0. small-ratio-numerics skill BINDS (tension is the small quantity;
  never infer sigma_0 by subtracting bulk energies — use the wall-localized
  integrand). Deliverables: sigma_0 enclosure, profile figure data,
  monotonicity/branch-identity checks, attempt record under
  proposals/P250-shell-bubble-clock/attempts/0002/.
Exact cross-checks to carry: sigma_0 = 2 int V_w dx = 2 int T dx must agree;
endstate residuals vs the certified bulk states < 1e-30-ish; first integral
constant along the solution.

### G2 — the bag (after G1)
Radial reduction: S = S(r), psi = f(r) on the slice; gradient density gains
the spherical angular term (derive the EXACT radial operator from
Tr(sum_i dS_i dS_i) before any code — analytic step first, per contract);
bulk states as in G1; the bag = 1D radial BVP at fixed omega in
(omega_*, omega_edge) with Q = omega I and R read from the profile;
verify R = 2 sigma/p against the solve and dE/dQ = omega along the family.
Deliverables: R(omega) table, radius/energy per charge, profile figure data.
New claim candidates (C-M5W-006+): bag construction + numerical tension —
formalize after results exist.

### G3 — global SOS certificate
Prove V0 - (omega_*^2/2) iota >= 0 on the slice (then full space) — closes
omega_c^2 = omega_*^2 and upgrades C-M5W-004's exclusion. Try exact SOS via
the known decomposition structure (lock square + tensor bound + W) before
reaching for SDP tools.

### G4 — bag stability (issue-sanctioned successor of S3)
Full fixed-Q Hessian spectrum about the constructed bag (G2 output);
small-ratio machinery: symmetrized operator, K-positive restriction,
lambda_min/lambda_2, zero-mode gauges, crossed refinement.

## Known gotchas (paid for already — do not re-pay)

- `maxwell_system()` differentiates at FIXED omega then substitutes; do not
  reorder (MC-3/MC-1 class).
- mpmath iv endpoints are intervals themselves; keep interval code
  all-interval; no mp.matrix/m midpoint casts inside iv mode.
- `np.trapz` is REMOVED: canonical code uses
  `substrate_framework.numerics.trapezoid_integral`; mutable scripts use
  `np.trapezoid`. Preflight for `getattr(np, "trapz")` too.
- Edits: re-anchor by content search, not line numbers (three stale-anchor
  incidents this campaign); after 2 failed patch rounds, full rewrite.
- YAML lists: quote any list item containing ": ".
- After ANY claims.yaml change: run scripts/render_docs.py AND
  scripts/render_memory.py, then scripts/validate_repository.py.
- Validation and commit are separate process invocations. Scoped validation
  covers: fixed checks + tests/test_m5_wall_clock.py +
  tests/test_m5_exterior_clock.py.
- Attempts are append-only; capture verifier stdout on FIRST execution.
- No issue/discussion comments (owner directive 3 above).

## Key artifact paths

- Proposal manifest: proposals/P250-shell-bubble-clock/proposal.yaml
- Memory contract: memory/vantasner/proposals/P250-shell-bubble-clock.md
- Attempt 0001: proposals/P250-shell-bubble-clock/attempts/0001/
  (derivation.md = full route history incl. superseded elimination routes
  and the falsified PSLQ shortcut; certify.py; maxwell_certificate.json;
  maxwell_point.json with 55-digit coordinates)
- Review record: proposals/P250-shell-bubble-clock/reviews/
  C-M5W-001-005-review.md
- Validation receipt: proposals/P250-shell-bubble-clock/validation-receipt.yaml
- Release: governance/releases/v0.172.0.yaml (+ current.yaml)
- PR: https://github.com/vantasnerdan/substrate-framework/pull/196
