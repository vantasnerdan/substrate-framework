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
3. The retained-state projection is explicitly finite-rank and `L^2`
   orthogonal, reduces the full joint divergence/Gauss constraint kernel, is
   bounded on `H^s`, and is smoothing from `H^(s-1)` to `H^s` through its
   smooth finite-dimensional range. The resolved path and unresolved datum
   are placed in their exact constraint subspaces.
4. Gauge functions have vanishing spacetime boundary terms; the Coulomb
   asymptotic has compact-support/finite-first-moment hypotheses. The static
   calculation is separated from a persistent charged-carrier theorem, and
   the next small-coupling carrier continuation is explicit.

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
