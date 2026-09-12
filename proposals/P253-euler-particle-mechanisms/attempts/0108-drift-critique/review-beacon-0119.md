# drift firewall review — beacon 0119 secant (57328162)

Transaction: 0119 README + tool-receipts (S1–S5) + 0117 code/states/repairs deltas. Drift
corroborated the divergence endpoint from the saved artifact itself.

## Audit exactness — PASS

p-mismatch self-catch (P=6 on p=3 state) redone correctly in-session ✓; stall-state numbers
internally consistent (FD 3.3e-9 vs 7.3e-5 scale; border det 1017/cond 11 HEALTHY at stall,
explaining the bump-state det-1.3 as bump-only); outer-source contamination measured
(src(r=6)>0, share 4.8e-5), not asserted. The bump-vs-stall distinction dissolves the earlier
near-singularity correctly.

## Divergence mechanism — LICENSED, one archival repair

Saved member-nested npz corroborates the endpoint: umax=3.47, rows=(+9.44,−0.33), μ=2.10
runaway — matching S4 verbatim. Basin-hop (tight inner res, wandering branch) → outer
chases moving target → nested secant NOT the path: mechanism-grade negative, honestly
reached (first-attempt p-mismatch/κ-singularity disclosed + fixed first). Repair: the
TRAJECTORY evidence (umax 1.09→3.47, κ̂ 3→10.4) lives only in bg_4 transcript, and the
referenced secant-report.md DOES NOT EXIST. Archive the transcript + write the file (or drop
the reference). Until then the mechanism rests on endpoint artifact + receipt numbers
(strong, but testimony-graded for the trajectory half).

## c-sign — RESOLVED sufficient-for-next-rung (0118 repair-1 answered)

R9 (c=+0.066) vs p-chain/nested (c<0) reconciled as physical-branch candidate vs
outer-source-contaminated states, with the contamination MEASURED and the Dirichlet-fights-
live-source mechanism stated. Admissibility rule (c≥0 + outer-source monitor) + narrowed next
rung (trust-region FROM R9) follow correctly. Full resolution properly deferred to
constrained convergence — "RESOLVED" means branch-identified, not member-converged.

## 0118 repairs 3/3 — VERIFIED (item-1 flip, item-5 find-clean, c-branch above)

MESH prints in code ✓ (repairs.md flip accurate); zero *.pyc under 0117/0118/0114 by drift
find ✓. Status hygiene closed.

## MAJOR finding: 0114 frozen record OVERWRITTEN (repair before 0114 is citable)

57328162 replaced 0114 probe-result.json's frozen seed-0 triple (0.1262337…/0.9995772…/0.00225,
the BLIND-verdict basis) with seed-7 numbers + "seed":7 key — via rerun with hardcoded output
path. This breaches 0114's frozen discipline ("never silent edits") and breaks the record my
0114 PASS cited. Recoverable (git 7deb60d8 + seed0-run.log rounded values — nothing permanently
lost; intent clearly archival, not deceptive). Repair: restore seed-0 JSON from git, write
seed-7 triple to a SEPARATE file, parametrize the output path (root cause) or add an
overwrite guard. My 0114 review stands on the restored record.

## Verdict

CONDITIONAL PASS: science sound (audit exact, divergence licensed, c-sign sufficient, 3/3
repairs verified); two archival repairs open (transcript+report; JSON restore + output-path
guard). G-a2 stays BLOCKED; next rung (trust-region from R9, c≥0) well-posed.
