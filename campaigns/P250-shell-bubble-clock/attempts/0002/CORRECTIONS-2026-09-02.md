# CORRECTIONS — 2026-09-02 (appended post-review; original files untouched)

Individual review of claim C-M5W-006 found two load-bearing defects and two
misstated diagnostics in the attempt-0002 record. Per append-only attempt
discipline this note records the corrections; the original artifacts above
are unchanged.

1. Tail asymptotics at the wrong bulk point: `sigma0_wall.py` (lines 79-81)
   substitutes `m = B_STATE[0] = c*` into the B-end Hessian `H_B`, but the
   certified wall endpoint is B = (0, c*, b*, f*). The archived
   `sqrt_lam_B = [0.830988, 3.356512, 3.594512, 7.060428]` is therefore not
   the linear asymptotics at the endpoint; the receipt's recorded values
   [2.5276, 2.899, 3.277, 7.012] are correct and are reproduced by the
   corrected evaluation (see `../0005/run_repair_stdout.txt`). The
   exponential-fit boundary treatment built on the wrong basis is
   superseded; the corrected cross-check agrees with Dirichlet to 3.7e-13.

2. Headline anchored to a failed solve: the archived headline row
   (h_crossing_L12, n_req = 3201) has solve_bvp status = 1 (mesh ceiling)
   and the acceptance predicate ignored solver status. The repaired
   headline is anchored to the status-0 rung (n_req = 1601):
   sigma_0 = 0.729298417862520, route spread 4.9e-13, budget total
   5.836e-10 (b1-dominated), bracket notation 0.72929841786(58). The
   archived headline sigma_0 = 0.7292984178707203 is superseded; it differs
   from the repaired value by 8.2e-12, which shifts downstream consumers
   (chi, E_crit) by <= 3.5e-11 relative — invisible at their recorded
   precision.

3. Diagnostics misstated in downstream prose (not in this attempt's files):
   `fi_drift` is an ABSOLUTE statistic max |T - V| on the collocation grid
   (4.08e-14 on the status-0 anchor), and monotonicity is a sampled
   criterion (zero reversals above 1e-8), not a continuum-exact statement.

Repaired evidence: `../0005/evidence_repair.py` with captured stdout
`../0005/run_repair_stdout.txt`, outputs `sigma0_results_repaired.json` and
`profile_L12_anchor.csv` (2001 samples of the n=1601 anchor; the archived
`profile_L12_n3201.csv` name is misleading — it holds the n=1601 solution).
