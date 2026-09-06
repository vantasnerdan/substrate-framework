**The full nonlinear rotational-Euler objective remains open. The campaign produced useful prepared linear-response results, but this fresh audit also identifies a missing proof in the compact-ring existence argument used by C-CST-018. My earlier unqualified statement that its geometry obligation was established was too strong.**

This discussion restores Federico's broader question, records the audit requested by Dan, and identifies concrete continuation routes. It distinguishes a supported conditional implication from a proved supplier, and a prepared linear response from an autonomous nonlinear constitutive law. The current registry still lists C-CST-018 as accepted in v0.183.0; this audit challenges its unconditional geometry inference. No historical campaign, registry entry, or release has been silently rewritten.

## 1. Which objective are we trying to complete?

Federico's original request was to derive Cosserat equations by coarse-graining Euler, with rotations carried by vortex structures. [Issue #198](https://github.com/vantasnerdan/substrate-framework/issues/198) records that intake; his [September 3 specification](https://github.com/vantasnerdan/substrate-framework/issues/198#issuecomment-5526924549) gives a linear micropolar target. Dan now supplies Federico's later request for the **full nonlinear rotational theory associated with Vikulin**, derived from Euler for the energetic fluid.

These are different targets. [Issue #200](https://github.com/vantasnerdan/substrate-framework/issues/200) and [PR #201](https://github.com/vantasnerdan/substrate-framework/pull/201) ended at a prepared, nonzero-wave-number, finite-window **linear-response and second-variation** statement. The terminal claim expressly excludes unrestricted nonlinear invariant-manifold and finite-amplitude stability results. Searching the P251 corpus and its effort memory found no Vikulin reference. Consequently the campaign did not establish Euler-to-Vikulin equivalence, regardless of whether its narrower theorem can be repaired.

Two distinctions matter for specifying the next objective:

- **Cosserat theory is not inherently linear.** Geometrically nonlinear formulations already allow finite rotations; a quadratic material energy is a separate constitutive choice. [Böhmer, Downes and Vassiliev, *Rotational elasticity*](https://arxiv.org/abs/1008.3833), supplies an explicit example.
- **An asymmetric effective stress is not by itself a violation of angular momentum.** With intrinsic spin and a couple current, the full angular balance controls the skew stress. The microscopic Euler momentum flux remains symmetric; deriving the coarse spin/current split is part of the work.

The exact Vikulin reference Federico intends is still unconfirmed. A primary candidate, [Erofeev, Pavlov and Vikulin, *Do rotational waves really exist?* (2018)](https://mpm.spbstu.ru/en/article/2018.59.7/), describes rigid rotating blocks in an elastic host and a one-dimensional sine-Gordon reduction. Those host elastic properties are inputs there. They are not already derived from Euler, and a block-chain reduction does not by itself specify an arbitrary three-dimensional finite-rotation continuum. Federico's intended paper/equations would resolve this target ambiguity.

Likewise, [Enciso–Peralta-Salas's knotted-tube theorem](https://arxiv.org/abs/1210.6271) is a stationary geometric existence result. It provides vortex topology in a decaying Beltrami field, not a rotational constitutive law or an automatic compact-velocity, finite-density assembly.

## 2. What was reviewed, and what survives?

The audit is pinned to main `b6fc902a0942d07996f12a81028fbd3f7c909a43`, release v0.183.0. It inventories **267 attempt directories and 268 Python scripts**, maps the original C-CST-001–007 routes, and audits the eleven accepted C-CST-008–018 claims at their proof and dependency boundaries. Three fresh **Sol High** reviewers separately reconstructed predecessor claims, geometry, and response/action/current arguments before comparing their findings with the earlier favorable reviews. Root checked the objective history, primary literature, exact separating examples, and the central geometry transfer.

This is a campaign-wide dependency audit with detailed examination of the proofs that determine the terminal result; it is not a claim that every line of every exploratory script has been formally verified. Historical session metadata confirms mixed model/effort settings, including Luna Low, but does not establish which model caused an error. The evidence below is mathematical. Substantive work in this audit used Sol High reviewers; Luna is reserved for optional literature discovery.

| Part of the campaign | Fresh assessment |
|---|---|
| Original straight-line/free-director mechanism and `alpha = L_v T / 6` | The unsupported route is preserved historically, not accepted as the terminal coefficient derivation. Free-frame rotation does not supply the claimed line-stretch locking energy. |
| C-CST-008–016 | No new equation-level contradiction found at their individually stated scopes. Exact ensemble identities, restricted quadratic actions, and prepared response constructions remain useful. Several coefficients or controls depend explicitly on the preparation. |
| Straight compacton, translation balance, radial nondegeneracy, radial edge estimates | Substantive positive geometry results survive. They do not yet supply the full perturbed free-boundary inverse. |
| C-CST-017/018 response, physical gain, KKS/Jacobi forms, observation connection and bulk current | The audited implications remain supported at their prepared linear scope, conditional on the actual background and at fixed nonzero wave number. |
| C-CST-018 finite-radius compact carrier | **Proof gap:** the coefficient-uniform, smoothly parameter-dependent bordered inverse needed for nonlinear existence is asserted without being constructed. |
| Full nonlinear Euler-to-rotational/Vikulin law | Not established. Finite rotations, nonlinear closure, Euler-derived interaction, and a controlled nonlinear evolution limit remain to be constructed. |

The earlier **2,677 passing tests** remain useful implementation regression evidence. They do not certify an analytic PDE existence theorem. This audit used source reconstruction and two small exact symbolic checks instead of repeating the full suite.

## 3. The new defect: a local estimate was promoted to an existence theorem

The critical transition occurs in [0263's trace estimate](https://github.com/vantasnerdan/substrate-framework/blob/b6fc902a0942d07996f12a81028fbd3f7c909a43/proposals/P251-cosserat-from-vortex-euler/attempts/0263/trace-estimate.md#L254-L320) and [perturbation/gluing argument](https://github.com/vantasnerdan/substrate-framework/blob/b6fc902a0942d07996f12a81028fbd3f7c909a43/proposals/P251-cosserat-from-vortex-euler/attempts/0263/perturbation-gluing.md#L60-L124).

The exact radial analysis retains the full normal operator and gives a useful high-mode source-to-edge bound proportional to `|m|^(-2/3)`. The later argument extends the problem through the edge, invokes local hypoelliptic estimates, and then claims an all-order tame inverse uniform under nonradial Hanzawa/profile perturbations. That last implication needs additional mathematics.

In particular, the cited [Bramanti–Zhu theorem](https://msp.org/apde/2013/6-8/apde-v6-n8-p01-p.pdf) supplies local interior estimates. Its authors explicitly leave higher weighted-order estimates untreated. It supplies neither this moving-boundary Poisson map nor the smooth tame inverse family needed by a nonlinear inverse theorem. Qualitative regularity and uniformly bounded low-order constants do not establish that family.

The missing construction is precise: define the actual domain, trace and parameter spaces; include the free graph, center gauge and translation/profile borders in one operator; prove its inverse with quantified loss and parameter bounds; then justify the nonlinear inverse theorem. Small spectral subspaces can persist under perturbation, while an exact translation kernel need not persist at a nonsolution iterate. The displayed Neumann formula does not resolve that distinction by itself.

[Review 0268](https://github.com/vantasnerdan/substrate-framework/blob/b6fc902a0942d07996f12a81028fbd3f7c909a43/proposals/P251-cosserat-from-vortex-euler/attempts/0268/review.md#L31-L52) repeated the inference without supplying those missing objects. The conditional implementation explicitly says it does not solve the finite-radius existence problem. Thus review repetition and green API tests did not close the analytic gap.

**This is a missing proof, not a counterexample to compact Euler rings or a no-go for the objective.** The minimum repair is to construct the stated inverse, preferably in explicit finite-regularity spaces if a no-loss estimate is available, or in a fully specified tame scale if necessary. Until then, the joint theorem's compact-carrier conclusion must be treated as conditional in scientific reuse. A separate registry correction or completed proof is needed to reconcile accepted status with this finding; this discussion itself does not perform a promotion or demotion.

## 4. Were there mistakes or shallow reasoning?

Yes, there are concrete examples; the campaign also contains real repairs and substantial calculations.

- **Incorrect perturbation scaling:** two `O(K)` cross blocks were initially treated as harmless around an `O(K²)` acoustic block. Their Schur contribution is also `O(K²)`. [0250's correction](https://github.com/vantasnerdan/substrate-framework/blob/b6fc902a0942d07996f12a81028fbd3f7c909a43/proposals/P251-cosserat-from-vortex-euler/attempts/0250/schur-correction.md) exposed this; the later cellwise construction resolves the relevant averaged cross row.
- **A test copied the wrong normalization:** the first 0270 API made spin inertia depend on kinetic energy rather than the tag second moment. The [repair record](https://github.com/vantasnerdan/substrate-framework/blob/b6fc902a0942d07996f12a81028fbd3f7c909a43/proposals/P251-cosserat-from-vortex-euler/attempts/0270/repair-review.md) adds an exposing sensitivity test. A green test that repeats a mistaken formula is not independent validation.
- **Missing centering terms:** the original tensor oracle omitted central-moment corrections. [0267 identified and checked the repair](https://github.com/vantasnerdan/substrate-framework/blob/b6fc902a0942d07996f12a81028fbd3f7c909a43/proposals/P251-cosserat-from-vortex-euler/attempts/0267/review.md). The repaired terms preserve the even Schur gain.
- **An unproved theorem transfer survived review:** the newly identified 0263→0268 inverse gap is the consequential outstanding example. The weakness is asserting that the difficult uniform inverse follows from local regularity, then reviewing that assertion as if it were a proof.

The fresh response audit also rejected one tempting false alarm: 0250's bounded normalized quadratic-gain map and 0265's singular full-history inverse are different maps. The latter genuinely costs `O(|K|^-2)`; the final construction fixes nonzero `K` before choosing carrier accuracy. That ordering is supported. An audit should correct its own mistaken objection as readily as an earlier mistaken proof.

## 5. Why even the repaired linear theorem would leave the nonlinear goal open

The prepared construction chooses a linear target history, synthesizes it from actual Euler/Lin directions, and matches quadratic forms on that selected section. This demonstrates a capacity of Euler's linearized dynamics. It does not show that one fixed coarse-graining of an open class of Euler states autonomously follows a particular nonlinear law.

A small exact example makes the distinction unavoidable:

`V1(theta) = k(1-cos(theta))`

`V2(theta) = k(1-cos(theta)) + beta(1-cos(theta))²`.

Both have linear stiffness `k`; their finite-angle restoring torques differ first at cubic order. Thus a matched linear optical gap cannot select sine-Gordon or a unique nonlinear interaction.

Nor is an orientation variable alone automatically closed. The exact incompressible Euler strain `u=(a x,-a y,0)`, with `p=-a²(x²+y²)/2` at unit density, changes a material blob's covariance eigenvalues. It deforms the blob rather than merely rotating its frame. This unbounded, nonperiodic example exposes what the local Euler equations alone allow; it is not a counterexample within a compact or finite-energy subclass or to the specially prepared construction. A nonlinear theory needs shape/inertia variables or a demonstrated slaving mechanism.

## 6. Concrete paths toward the entire objective

**A. Derive exact nonlinear coarse balances, retaining shape and memory.** Start with the actual Euler momentum, intrinsic angular momentum, covariance, pressure torque and boundary transport. For a commuting filter, the exact subfilter stress is `tau = overline(u⊗u) - U⊗U`; its evolution belongs to a hierarchy. Fix the observation map and derive the unresolved terms explicitly. The first decisive test is whether two Euler states with the same proposed coarse state have different coarse accelerations. If they do, add the missing internal variable or retain memory. A local law requires a further closure or scale-separation theorem.

**B. Construct a genuine finite-rotation Euler reduction.** Use a physical `R(X,t) in SO(3)`, deformation, and necessary shape variables. Pull back the full Euler action on a justified volume-preserving vortex family, including pressure and ambient corrections. Compute cubic and quartic terms and the residual normal to that family; then prove control under unforced nonlinear Euler over a declared amplitude and time regime. Repair the custom compact ring first, or change to a published compact-shell supplier and redo the response on that same field. Restricting an action alone does not prove dynamical closure.

**C. Reconstruct the intended Vikulin equations and derive their inputs.** Freeze the specific primary source, variables, finite-angle potential, host assumptions and stress/current convention. Derive its inertia and interaction from Euler before comparing coefficients. A one-axis sine-Gordon reduction can be a useful bounded target, but the full three-dimensional objective still requires finite-rotation compatibility and closure of shape and transverse modes.

My recommended order is **A and target identification in C first, with B as the constructive continuum route**. This earns exact nonlinear information before committing to another long preparation campaign. The compact-inverse repair is a separate obligation; solving it would restore the narrower theorem, not automatically finish the nonlinear objective.

A successful full derivation would provide one observation map and one encoder, independent of the desired future history, for a stated open class of macroscopic initial data. Unforced Euler evolution would converge to an objective nonlinear rotational continuum, with compatible action and angular-current limits and coefficients fixed by the microscopic construction. The time regime can be finite; global regularity of arbitrary three-dimensional Euler need not be assumed.

These are concrete plausible routes, not a promise that the full equivalence exists and not an exhaustion verdict. Federico's exact Vikulin reference is the remaining target clarification.

## Audit record

Detailed reports, source line anchors, campaign/claim inventories and exact-check receipts are preserved in **P252, attempt 0001**:

- [Geometry proof audit](https://github.com/vantasnerdan/substrate-framework/blob/c6e950df1e5b97c40ea8e8c631a7725ee60c9980/proposals/P252-nonlinear-euler-objective-audit/attempts/0001/geometry-audit.md).
- [Response, action and current audit](https://github.com/vantasnerdan/substrate-framework/blob/c6e950df1e5b97c40ea8e8c631a7725ee60c9980/proposals/P252-nonlinear-euler-objective-audit/attempts/0001/response-audit.md).
- [Predecessor claim-by-claim audit](https://github.com/vantasnerdan/substrate-framework/blob/c6e950df1e5b97c40ea8e8c631a7725ee60c9980/proposals/P252-nonlinear-euler-objective-audit/attempts/0001/predecessor-audit.md).
- [Target reconstruction and nonlinear routes](https://github.com/vantasnerdan/substrate-framework/blob/c6e950df1e5b97c40ea8e8c631a7725ee60c9980/proposals/P252-nonlinear-euler-objective-audit/attempts/0001/target-and-routes.md).

The audit adds no accepted claim and leaves immutable P251 records intact. The publication receipt links the exact audit checkpoint and this discussion.

Audit source checkpoint: [`c6e950df1e5b`](https://github.com/vantasnerdan/substrate-framework/commit/c6e950df1e5b97c40ea8e8c631a7725ee60c9980).
