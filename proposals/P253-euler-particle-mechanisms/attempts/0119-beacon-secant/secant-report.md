# Secant report (beacon 0119 — closes the obligation's third item)

## Nested secant v1 (unthreaded, unguarded): DIVERGED, mechanisms recorded

Inner Newton (hardcoded P=6) from p3-warm state fell to the trivial root
(pde_res 2.6e-17, umax=0); rows went singular (κ̂≈0.015, rbar garbage
+4.94); outer secant drove μ→−15.8, c→1.7, umax→11.6, pde_res→1e7.
Two independent defects: (a) p-mismatch (inner P=6 vs warm p=3);
(b) no κ-floor guard (secant on singular rows). Both fixed before v2.

## Nested secant v2 (p-threaded, κ-guarded, p=3, p3-warm): DIVERGED with
mechanism (bg_4 delivery, this session)

Inner converges TIGHT (pde_res 1e-8..1e-10) but basin-hops branches:
umax 1.09→3.47, κ̂ 3→10.4 across 8 outers; outer secant chases a moving
target (μ 0.35→2.10); final rows (+9.44, −0.33). Branch non-uniqueness
at fixed (μ,c) DEMONSTRATED: tight PDE residual does not select the
member branch. Nested secant is NOT the path (retired with evidence).

## Standing direction (unchanged)

Trust-region bordered Newton from the R9 state (c=+0.066, kap=1.018,
rbar=0.996) with c≥0 bound + outer-source admissibility (0119 README).
G-a2 still BLOCKED. No verdict altered by this report.
