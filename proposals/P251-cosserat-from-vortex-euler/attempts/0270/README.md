# 0270 — compact-ring materialization contract

Owner: `herdr optical-review pane w3:p3`, independent of root-authored 0265
and C-CST-018.  This is a bounded additive API transaction; no claim
promotion, central/governance edit, or commit is authorized.

## Frozen objective

Materialize importable exact axisymmetric Grad–Shafranov definitions for a
supplied symbolic streamfunction and labels, including:

1. cylindrical velocity and pressure expressions;
2. literal Euler residual components derived from the cylindrical metric and
   centrifugal factors;
3. residual/domain checks with an independent analytic solution; and
4. positive isotropic tag-inertia and density normalization.

The API is conditional on an actual supplied solution.  It must not claim to
compute the finite-R inverse, establish PDE existence from a symbolic
placeholder, or replace the 0268 finite-R theorem.

## Inventory and disjoint writes

After central activation, inspect the existing inventories
`src/substrate_framework/euler_compact.py` and `src/substrate_framework/euler_joint.py`
and reuse any applicable definitions rather than duplicating them.  Create
`src/substrate_framework/euler_compact_ring.py` only if the inventory lacks
the requested compact-ring abstraction.  Add focused tests only in the
corresponding test module, differentiating the actual cylindrical Euler
residual and exposing wrong metric/centrifugal factors.  All new review,
receipt, and validation evidence stays under `attempts/0270`; source writes
are limited to the additive module/test surface above.

## Contract and checks

The implementation must state its streamfunction convention, cylindrical
coordinate domain (`r>0` where required), label regularity, and pressure
normalization.  Tests must include an independent analytic Grad–Shafranov
solution, exact/known residual checks, domain rejection, and positive finite
isotropic inertia/density normalization.  A symbolic placeholder with
nonzero residual is not accepted as a solution certificate.

Central/schema activation of 0270 is required before substantive inventory
inspection or source changes.  No full suite, central/registry edits, or
commits.  The transaction verdict concerns only additive algebra/API
materialization and its focused tests; finite-R inverse and downstream
response/action/current joins remain outside scope.

`maximum_verdict`: conditional exact compact-ring definitions and tested
normalization, with all PDE/existence and parent-join boundaries explicit.
