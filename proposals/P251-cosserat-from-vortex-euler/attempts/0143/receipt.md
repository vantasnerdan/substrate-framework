# 0143 source and exact construction receipt

Baseline bf508e1 (preregistration began at a93e78e); accepted v0.175.0
unchanged. Global stationary isotropic compact-cell assembly is
established; rigid-pose elasticity on that family is refuted by its
exact stationary Kelvin degeneracy. Parent objectives remain active.
The strongest evidence is the analytic construction and source-licensed
smooth localization, not a sampled or numerical null spectrum.

Primary source: Constantin--La--Vicol,2019, Theorems1/2, equations18--23
and the annular cutoff following108; source PDF SHA256
1c945acecb8242fd686d546f4d0876522d39122f1666e2d265473304d3af275b.
The nonsmooth central pre-cutoff circle is explicitly excluded. The
2026 symmetry preprint was discovered during the registered localization
search; its exact Theorem1.1 hypotheses are recorded as a route-class
limitation, not imported as a parent obstruction or global smooth no-go.

First execution first.stdout preserves5 passing predicates and one
translation-identity failure, exit1. Diagnosis: SymPy left derivative
substitution nodes unevaluated; even simplify printed a distributive
zero with rho factored differently. Evaluating those derivatives with
doit before simplification gives the exact unchanged full Euler
linearization. repaired.stdout preserves8/8, exit0; Ruff passes.
No tolerance, physical term, domain or scientific predicate was weakened.

    PYTHONPATH=src python proposals/P251-cosserat-from-vortex-euler/attempts/0143/verify.py

SHA256 global-localizable.md:
1d4e81fdb4b0303c79b02e1c5c9be50cd33e7170fe05ed3a869ace59c28b4312

SHA256 verify.py:
d583a3639bdfcd33ef272936fb40647d9832635e531f11bcce8985d3cdff9e5b

Checkpoint support: bf508e1 was pushed. selector.stdout confirms the
already passing0136 scoped workflow (tests/test_euler_forcefree.py).
GitNexus refreshed atbf508e1:46026 nodes,72166 edges,758 clusters,2 flows.
Its new-symbol impact still omits the test/attempt imports found by rg;
the explicit two-consumer boundary in0136 remains the validation basis.

Continuation:0141 same interacting array optical response;0142 exact
constant-factor EPS-compatible optics;0144 generic3D coupled mean
response. No terminal PR or campaign-exhaustion decision follows.
