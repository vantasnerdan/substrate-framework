# Persistence artifact (beacon, 0157/persist, FIRED-ON-#87)

`persist/persist_defect.py`, exit 0, tracked path verified
(check-ignore negative before first run). Frozen pre-compute
thresholds all met — no tuning post-freeze.
Perturbed defect: E0=29.04788 -> Efin=28.95146 = E_clean to 5dp
(<2% bar), charge 0.99308 (<1% bar). Trivial-background same
bump: 0.10830 -> 0.00000 (<5% bar). Claim scope (honest): linear
texture stability of the far field + charge-readout robustness;
topological protection is analytic input (winding invariant under
continuous deformation), not numeric output. Core disk r<1 pinned
as declared sub-continuum regularization.
Mid-build honest failures banked: (1) raw closed-loop diffs
telescope to 0 identically — charge lives in the unwinder
accumulator (end-minus-start after 2pi-scaled unwrap); (2) energy()
cut-wall artifact noted (28.95 vs continuum ~0.3 — comparisons
only, never absolute); (3) free-u gradient flow ejects the defect
through the boundary (surface-annihilation channel, consistent
with D1 moral) — hence the w-flow formulation.
Fence honored: demonstrated, ALIVE never declared.
