# 0095 — exact Kelvin complement and physical-observation repair

Owner `/root/smooth_core_review`; P251 / issue #198 child, this directory
only. Original issue #198 was reread through `gh issue view`: its objective
is an exact conditional Euler coarse-graining with named independently
falsifiable ensemble premises. Proposal N3 expressly includes affine
transport and excludes nonaffine relaxation. That original distinction is
retained; this attempt does not add a nonlinear invariant-manifold demand.

Inputs: 0091's actual material action and spin map, 0084's Kelvin/Lin
identity, 0085's compact-velocity isovortical pair, and 0087's defined
volume-preserving mean. This author owns no input or canonical changes.
No old geometric mass is treated as reduced Euler inertia.

Frozen positive deliverable: derive the exact phase-space complement and
its physical field/mean map; identify when the declared nonaffine exclusion
licenses a conditional constitutive action and what it does not say about
unrestricted Euler reconstruction. Execute a complete algebraic repair of
restriction-versus-reduction using an exposing coupled Hamiltonian example,
including its actual retained momentum. No numerical sign or fitted closure.

Candidate A: first perform full Kelvin reduction, then keep and eliminate
the nonaffine complement in the same symplectic/Hamiltonian action.
Candidate B: explicitly freeze that complement as the original Cauchy–Born
variational closure, compute the discarded residual and bound its effect.
The two are compared, not conflated. Selection is by physical observation,
unchanged Euler/Kelvin conventions and exact action identities. Main owns
the compact-route continuum joining; 0097 owns tagged cotangent reduction.

## Completed receipt

`complement-and-observation.md` derives the exact residual, full retained
complement evolution, frequency-dependent Schur action, and matching
physical observation map. It constructs the actual material Lin/tag/mean
composition and gives a finite-time residual error bound. The original
exclusion of nonaffine relaxation is identified as the pedigree of the
conditional constitutive action, not treated as a newly discovered defect.

An explicitly positive four-dimensional Hamiltonian exposes the difference
between frozen and dynamically eliminated nonaffine motion: inertia `1/2`
becomes `(35-omega²)/(63-2omega²)`, with actual momentum transformed by the
same expression. This is an algebraic oracle, not claimed as an Euler
profile or counterexample to the particular EPS compact sector. Its full
physical field correction and initial Taylor discrepancy are computed.

Fourteen exact checks pass on the first execution, retained in `stdout.txt`;
the verifier uses the importable noncommuting Schur-jet API. Ruff and scoped
diff checks pass. `route_verdict: established`; `evidence_scope:
REPRESENTATION_SCOPED` for the Euler application. The parent can use the
fixed-projection conditional constitutive action without an invented
nonlinear-invariance requirement, while keeping any stronger free-Euler
or actual-mean identification tied to the computed complement/field map.
