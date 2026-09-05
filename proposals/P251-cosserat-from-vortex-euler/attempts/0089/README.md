# 0089 — primary-source-guided positive material stiffness witness

Owner `/root/smooth_core_review`; constructive child of P251 / issue #198,
this directory only. Base accepted release `v0.171.0`; exact functional
import from 0084, not a new MHD dynamics import. Parent objective unchanged.

Frozen positive deliverable: extract a concrete smooth divergence-free
periodic displacement from primary force-free field evidence, repair only
printed coefficient rounding if necessary, and evaluate its Euler Jacobi
stiffness exactly. Fixed candidate source is East–Zrake–Yuan–Blandford,
arXiv:1503.04793v3, Appendix B.2, Table 2. Their rounded finite Fourier
coefficients are trial inputs, not empirical fitted Euler coefficients.
The sign is recomputed from the actual field and exact functional. Other
tables labeled compressible do not qualify as divergence-free witnesses.

Candidate route A: exact rational divergence-free projection of Table 2.
Candidate route B: its registered odd Fourier lattice on the two-wave ABC
field, using exact algebra if printed rounding obscures the sign. The
periodic-boundary route is preferred to the separately suggested cylindrical
Lundquist route because it avoids importing confining-wall surface energy.
No empirical comparison, solver instability theorem, or cylinder boundary
assumption enters the selected result.

The strongest oracle is rational Fourier convolution and Parseval applied
to `K/rho=lambda<F.curl F>-<|curl F|²>`, `F=Xi cross u`. No numerical
eigenvalue, truncation-limit claim, MHD growth rate or Euler stability
conclusion is inferred. A fixed positive rational witness requires no
soft-eigenvalue numerics. Parent 0088 owns broader witness search; this
attempt returns a single reproducible candidate with explicit conventions.

## Completed receipt

`positive-witness.md` specifies the exact projected rational Table II trial
and its source applicability. It establishes
`K/(rho V)=5008301/1250000000>0` by both the curl functional and independent
pressure-Hessian/advective evaluation. The periodic material metric is
`M/(rho V)=14496029/125000000`. Eight exact checks pass on the first run;
`stdout.txt` preserves the execution. Ruff and scoped diff checks pass.

`route_verdict: established`; `evidence_scope: symbolic_verified` for the
fixed positive Fourier witness, not MHD dynamics or Euler spectral stability.
The sibling 0090 owns compact localization and transfer to the same EPS
background. Actual material/Kelvin joining remains in the parent work.
