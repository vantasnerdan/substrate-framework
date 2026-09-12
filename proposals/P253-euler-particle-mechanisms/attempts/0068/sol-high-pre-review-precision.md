# Sol-High pre-review precision receipt

Date: 2026-09-06

The final author pass repaired two presentation/domain issues before independent
review:

1. The gauge action now separates the free field Lagrangian from the material
   coupling.  The joint material-map action includes that coupling exactly
   once; the Euler--Lorentz and Maxwell equations are unchanged.
2. The characteristic speed applies to the Maxwell subsystem.  The joined
   incompressible system retains its elliptic pressure projection and therefore
   has no strict finite propagation cone for the complete state.
3. The retained-state projection is a bounded finite-rank idempotent with
   smooth range. Both `Pi` and `Q=I-Pi` are bounded on the dimensionally typed
   product `H^s` and `H^(s-1)` spaces and reduce the full joint
   divergence/Gauss constraint kernel; the finite-dimensional range supplies
   the stated `H^(s-1)` to `H^s` smoothing. Orthogonality across the
   dimensionally mixed `(u,chi,E,B)` state is neither assumed nor used. The
   resolved path and unresolved datum lie in their exact constraint
   subspaces.
4. Gauge functions have vanishing spacetime boundary terms; the Coulomb
   asymptotic has compact-support/finite-first-moment hypotheses. The static
   calculation is separated from a persistent charged-carrier theorem, and
   the next small-coupling carrier continuation is explicit.
5. The small-coupling route freezes the material tag leaf while
   `Q_g=g integral chi` varies with `g`; only variations at a fixed branch
   member use a fixed `Q_g`. The coefficient ledger records the gauge-field
   and tag normalization freedom and the invariant combinations
   `epsilon*mu` and `g^2/epsilon`.

The arbitrary-tag defect generated that continuation after preregistration:
an invariant core label, sub-gauge-speed carrier, comoving elliptic Maxwell
solve, and `O(g^2)` KKS/Hessian correction on each `Q_g` level while the tag
leaf is held fixed. The Noether translation row remains a cokernel identity to
derive. The route is recorded as blocked on the 0066 carrier complement,
stabilizer-preserving tag profile, and joint localization coercivity, not as
an established branch.

The exact algebra/API predicates did not change, so their existing captured
receipts were retained rather than rerun.  Current post-edit hashes are recorded
in `artifact-hashes.sha256`.
