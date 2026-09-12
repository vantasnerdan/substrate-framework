# P0 source/foundation map — tool-cited (beacon 0108)

Every row cites the tool output that produced it (receipt IDs in `tool-receipts.md`).
Direct re-verification in THIS pass vs inherited reviewed scope is stated per row.
Primary PDFs not re-fetched here are named as missing-source boundaries, not asserted.

## Retained campaign state (pause findings, verified at source this pass)

| State | Content | Tool citation |
|---|---|---|
| 0042/0045 retained | Exact finite-rank resolved/complement Euler state, same-field pressure, causal unresolved map `W[v,w0]`; joined result after one bounded correction: transported bounded region / finite-moment tag, Sobolev finite-rank projection, `P_obs,Q_obs` distinct from velocity projection, ambient momentum conditional on `L1`/decay | T4 (`memory grep 0042` hit: `codex/proposals/P253-euler-particle-mechanisms.md:145`); T11 (`attempts/0042/` exists with `result.yaml`, correction receipts) |
| 0032/0038/0039 Gavrilov | Fixed compact Gavrilov carrier carries an exact whole-space DA linear semigroup lower bound `||S(jT_*)||_ess >= lambda_+^j`; adverse evidence for that carrier as P2 sector; linear theorem, not nonlinear instability | T5 (`memory search` hit `codex/proposals/P253-euler-particle-mechanisms.md`, hyperbolic-carrier section) |
| 0095/0104 thin ring | A–F + H established: fixed-member compact-interior Cao transfer, returned DA phase, fixed-`J_loc` observed essential-norm growth `||J_loc S_g(jT_*)|| >= c_obs lambda_+^j` for each fixed sufficiently thin member + small signed charge under `C_mon(delta)g^2 < c_HF delta`. G (exact fixed-row leaf) and I (complementary charge) BLOCKED. Nonlinear fixed-leaf closure missing | T9 (`attempts/0104/review.md`, `verdicts.yaml`, full unit table) |
| 0107 interrupted | Route A (relative-momentum leaf `J_ren = Delta I_h + Delta P_EM`, compact witnesses, corrected initial curve), Route B2 (compact-edge interface transfer open), Route C (unwritten Hill Hessian proof) are author drafts, not conclusions | Resume file `attempts/0107/pause-state.md` (task-mandated resume input) |

## Supplier reconciliation (frozen-issue table order)

| # | Source | Exact hypotheses (as usable here) | Perturbation class / norm | Transfer to establish | Standing in this pass |
|---|---|---|---|---|---|
| S1 | Choi 2020, Hill stability (arxiv 2011.06808) | Translating Hill vortex; orbital stability from energy+impulse+vorticity constraints | Axisymmetric, swirl-free class in the source; full-3D extension unearned | Perturbation class + stability norm; extension beyond axisymmetric/swirl-free when needed | No coercive full-3D Hill Hessian asserted: `euler_p2_principal.py` header states its helpers assert no such Hessian, evidence = active 0095 (T8). Route-C Hill result unwritten at pause (0107 pause-state). |
| S2 | Cao–Lai–Qin–Zhan–Zou 2022/23, steady rings | Thin-core vortex-ring family; variational stability mechanism; leading jet `mu=3/8·k·r·L/pi`, `c=k·L/(4·pi·r)`, `I_z=pi·rho·k·r²`, ordered `(mu,c)/(kappa,I_z)` Jacobians | Fixed-member compact-regular core tube; weighted global input `X_in^s` (`<x>^{2σ}` A2, `1/2<σ<3/2`) → local same-`H^s` observed output | Actual conserved-state neighborhood, full physical modes, internal identity, interactions on the SAME family | SELECTED carrier. Fixed-member compact-interior transfer established (0104 Unit C, T9); reusable jet encoded in `euler_cao_schur.py` (T7/T8). Endpoint/collar/uniform-global-C1 NOT claimed. |
| S3 | Dávila–del Pino–Musso–Wei 2022, leapfrogging | Smooth Euler interacting-ring solutions + derived asymptotic interaction law | Finite-window reduced law; exact-solution vs reduced-law distinction load-bearing | Join selected stable carrier; preserve exact-vs-window distinction | Finite-window law only; join with S2 carrier is prospective (no new construction this pass). Primary PDF not re-fetched here — missing-source boundary. |
| S4 | García–Hassainia–Hmidi 2026, periodic leapfrogging (arxiv 2603.21644) | All-time periodic interaction in translating frame for constructed axisymmetric family | Parameter/regularity hypotheses of that family; periodic existence ≠ general stability | Review preprint + hypotheses; no stability inference | Not re-verified this pass — missing-source boundary; no conclusion drawn from it. |
| S5 | Slobodeanu 2015/2019, steady Euler ↔ Faddeev–Skyrme+mass | Exact STEADY correspondence with strongly coupled quartic model + potential | Steady admissible variations only | Global domain, potential/input accounting, time-dependent action + stability/interaction transfer | Time-dependent transfer UNEARNED (frozen-issue warning preserved). No construction this pass. Primary not re-fetched — boundary. |
| S6 | Slobodeanu 2019, Euler flow on S³ + FS solution | Explicit topological example | Closed 3-sphere geometry | Sphere ≠ isolated Euclidean particle; transfer needed | No transfer claimed. Primary not re-fetched — boundary. |
| S7 | Faddeev–Niemi 1997, knot solitons | Topological-soliton mechanism + energetic structure | Chosen field-model stabilizers | Which terms/constraints/dynamics actually follow from Euler; topology ≠ particle/quantization | No Euler realization claimed. Primary not re-fetched — boundary. |
| S8 | Gavrilov 2018 + Constantin–La–Vicol, compact existence | Compact-velocity steady existence (if that geometry is useful) | Existence class only | Stability + interaction remain separate; strict compact velocity is a choice, not the particle definition | Adverse linear evidence on the fixed compact Gavrilov carrier (T5); no P2 persistence claimed. Primaries not re-fetched — boundary. |
| S9 | Choi–Jeong 2021, filamentation near Hill (arxiv 2107.06035) | Geometric deformation coexisting with orbital stability | Exposing example | Chosen identity/stability observable must detect this failure mode | Design constraint on P1 observables (recorded in P1 file); primary not re-fetched — boundary. |

## Where structure would originate (no invention; open slots named)

- Localization: S2 compact-core ring + impulse-defined far field (`euler_impulse.py` dipole, T8). No wall/trap/pin imported.
- Relativistic structure: NOT derived. `euler_scale_causality.py` header: similarity weights do not turn a carrier band into a causal cone (T8). P4 bridge unearned.
- Quantum structure (spin-1/2, statistics, charge coupling): NOT derived. No module asserts them; P4/P5/P6 unearned. Any S5/S7 suggestive structure stays a hypothesis until time-dependent Euler transfer exists.
- Interaction: S3 finite-window law is the candidate supplier; same-family derivation from the full action/current is the open P3 construction. `euler_impulse.py` cross-energy is leading-order algebra, not a mechanical force (module docstring, T8).
- Scale selection: open. Similarity weights (`euler_scale_causality.py`, T8) expose the scaling boundary; no substrate parameter fixes R/kappa.

## Missing-source list (blocked direct re-verification; independent S2/P1 work continues)

S3, S4, S5, S6, S7, S8, S9 primary PDFs; HF primary already independently checked at pause
(SHA-256 `c6a35c44…` per 0104 review, T9 — inherited, not re-fetched here). No conclusion in
this file depends on the un-fetched primaries.
