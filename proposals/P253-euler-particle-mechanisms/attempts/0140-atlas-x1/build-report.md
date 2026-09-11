# X1 BUILD report — sign computation attempted, BLOCKED (atlas)

Charter: shepherd X1 BUILD under frozen design (falsifier/stop/cost) +
drift pin (sign = reproduction) + N4-order line. Method: paper + banked
numbers only, per scope. No new solves executed.

## Computation performed (banked A3 monodromy, all m-keys, both eps, N64+N128)
Eigendecomposed every banked per-m block (12×12): ALL multipliers satisfy
|ρ| = 1.000000 exactly, phases in ± conjugate pairs (e.g. m2: ±2.2605;
m1 e5e-07: ±0.4282, ±2.9811; m3: ±0.6820; m4: ±0.9269 — full table in
build log below). m1 e1e-06 N64 block is identity-like (all args 0).

## Informative result: single-point data CANNOT supply the sign (structural)
Hamiltonian monodromy spectra pair-cancel: every +arg has a −arg mirror,
so no preferred fast-rotation sense exists at any single slow point. A
Berry/Hannay sign needs M_m(R) around a slow circuit (pair braiding),
which one-point data cannot supply at any analysis depth. This is
symplectic structure, not a data-quality gap — reanalysis cannot fix it.

## Verdict: X1 BLOCKED (not dead) — two named missing pieces
- B-0071 (frozen D-0071, anti-fudge honored): joint slow-collective orbit
  chart at working order. The κ-fit parametrization is a fit, not a chart;
  not substituted. Shared with R-B — one construction serves both.
- B-CIRC: monodromy at ~8 slow-circuit stations (assemblies via run_a3.py
  machinery, cost unverified — may exceed assembly scope; needs shepherd
  costing before anyone runs).
Falsifier did NOT fire (no sign computed either way); flat-connection
death not evidenced either. X1 stands scoped+pinned, awaiting B-0071.

## Coordination
N4-zero-test-first respected: no N4 duplication attempted (N4 state per
cipher N1/N2/N4/N5 ranking — N4 first, status open; X1 extends, not
duplicates, when both move).

## Options for shepherd
(a) Charter 0071-chart construction (serves R-B revival + X1 circuit);
(b) stand X1 down to conditional-sketch until the chart exists elsewhere.
Recommended: (a) — the chart is the lane's highest-leverage shared object
(0071 named it missing; R-B, X1, and 0127-reopen all price it).

## Build log (reproduction)
`python3 -c` eigendecomposition over Mono_{e1e-06,e5e-07}_{m01,m02} +
Mono_N128_e1e-06_m06 (keys m1–m6 as available); |ρ|−1 = 0 to 6dp all;
phases paired ± to 4dp all. Script: rerun the one-liner in session
receipt; inputs are cipher banked npz (unmodified).
