# 0035 — smooth-continuum continuation contract

Parent objective: the original P251 N1–N7 construction, including finite-core
Euler dynamics, reaction-bearing translation and stationary EPS compatibility.
User reaffirmed on 2026-09-04 that commits are checkpoints, not stopping points.
The original exact objective remains frozen; no asymptotic scope reduction was
authorized. Base: campaign v0.171.0, checkpoint 3626fbf; accepted main v0.174.0.
No accepted claim changes are proposed by this attempt.

## Candidates and typed obligations

Continue the existing candidate universe append-only:

1. Smooth desingularization of the six-vortex relative equilibrium and its
   isovortical slow optical sector (0036). Primary theorem applicability is a
   construction input only when its actual hypotheses match the object.
2. Exact planar stationarization by a uniform-vorticity shift, compared with
   stationarity in rotating coordinates and with a genuinely 3D EPS field.
3. Material-parcel center and locked-inertia decomposition of the Euler action,
   keeping pressure transfer, shape motion and the full advected constraints.
4. Stationary 3D Beltrami/EPS construction with a directly derived collective
   action if planar continuation cannot supply the required compatibility.

Selection retains the proposal's criteria: Euler equations, frame covariance,
isotropy, parameter economy, correctly typed limits, explanatory reach. No
empirical comparator or fitted coefficient is used. Existing EPS sources and
classical Euler are permitted inputs. New desingularization theorems are
candidate imports to be checked and explicitly inventoried, not accepted canon.

For (1), requires: finite-core equilibrium in the specified domain and full
isovortical variations; pass licenses: finite-core existence and only the
spectral persistence actually proved; does not license: material translation,
3D EPS embedding or continuum homogenization; unlocks: coupled action.
For (2), requires: planar incompressibility and a relative equilibrium; pass
licenses: actual stationary Euler field with explicitly changed background;
does not license: unchanged far field, finite total energy or arbitrary knotted
3D EPS tubes; unlocks: compatibility audit of this candidate.
For (3), requires: a smooth material flow and transported finite parcel; pass
licenses: exact kinetic/action and pressure-reaction identities; does not
license: an invariant finite-dimensional ansatz or closed Cosserat law;
unlocks: elimination/closure of the residual shape sector.
For (4), requires: the exact EPS theorem hypotheses and collective-variable
map; pass licenses and downstream closure follow only from the derived map.
Maximum verdicts are route-scoped exact identities/applicability results until
the full dependencies are met; failure scope never propagates to the parent.

## Analytic specification and verification

Use Euler's variational equations and exact calculus first. Define all fields,
constraints and complete variations before a Hessian reduction. Distinguish
material mass center from circulation centroid. Distinguish changing the
physical background from coordinate rotation. Existing sphere moments and
collective-coordinate APIs are context, not a license to insert missing mass.
Small-ratio skill prescriptions apply to any eventual spectral remainder;
currently there is no licensed production numerical remainder. Exact symbolic
identities require no mesh, zero-mode floor or eigenvalue tolerance. Any
exploratory computation remains explicitly hypothesis-generation evidence.

Ownership: main owns 0035 and 0037 plus proposal/memory; analytic worker owns
0036 only. Existing attempts stay immutable. Reusable successful definitions
will enter importable modules and targeted tests after their meaning is fixed.
No terminal PR refresh or promotion follows from a checkpoint.

Status: active. Route verdicts and first-run outputs are recorded with each
construction. Next action: validate the extended proposal, then execute the
parallel finite-core and Euler-action constructions.
