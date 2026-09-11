# Receipt: D3 multi-period (FROZEN F-bar, coherence-corrected pre-run)

Command: `python3 run_d3.py` (94.5 s). Env: CPython, numpy; A3 backbone N=64.
F-bar: 11-d3/00-fbar.md (freeze + anti-rigging correction, both pre-compute).
Numbers: closure 9.9e-11→1.2e-10 all 5 periods (H=5 ✓); transverse phase
0.00016 ×5, std 0.00e+00; overlaps 1.0000. Verdict per frozen bar:
ALIVE-leaning (coherent substrate).
TAUTOLOGY CAVEAT (strength-killer, disclosed not buried): under exact closure
the per-period monodromy is the SAME matrix, so identical phases are
AUTOMATIC — coherence was guaranteed by H=5, not discovered. The verdict is
technically per-bar but substantively near-vacuous: what D3 actually measured
is (a) horizon H=5 (orbit integrates cleanly 5 periods — genuine, supports
any multi-period program), (b) transverse phase 0.00016 at N=64/dt=0.005 vs
0.384 at NQ=200/dt=0.004 — THIRD marginality signature (pair position moves
O(1) across discretizations; Krein story holds).
In-run repair (disclosed): flow→flow_frac (shoot solves window-exact; plain
flow gave 1.9e-3 first-period mismatch — the A3 window lesson re-learned).
Proposed D3b (NOT run — needs charter/f-bar): adiabatic-Γ-drift Berry run
(ramp Γ over periods; non-periodic parameter drift makes Berry accumulation
nontrivial). The geometric question stands open; D3-periodic cannot see it.
