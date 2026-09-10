# 0124 NO-GO certificate: state-error-transport lemma route (beacon)

Charter: shepherd ruling on 0123. Instrument: sharp_dQ.py (Weyl-exact
3x3 dQ + exact nonlinear transfer + P2 quadrature probe).
## Measured (all exit 0; CORRECTED 2026-09-10 — see note)

CORRECTION: G1 first printed 0.87/0.69 via a spurious /6 in
dzeta=(jf/6)du (zeta = EPS^-2 s^6 = f, so dzeta/du = jf exactly).
Caught by the asymptotic-slope check (exact/sc -> 5.21, not 0.87);
fixed in sharp_dQ.py, rerun below. Exact-transfer (G2') and sweep
numbers never used /6 — unaffected.
- Fitted state: G1 = 4.16 (charter gate <= 0.9: FAILS; margin x2.6
  vs lam 10.77); G2'-exact central = 9.10 (FAIL).
- Trust-r3 (best, res 6.6e-3): G1 = 5.21 (charter gate: FAILS;
  margin x2.1 vs lam 11.16); G2'-exact central = 9.18 (FAIL).
- S2 quad probe: P2-R rms 8.9e-3 vs P1-R 10.3e-3 (x1.16 — clean,
  quadrature is NOT the floor).
- Floor robustness: res ~1e-2 across 3 meshes (coarse/fitted x2),
  P1/P2, 140+ iters (Newton + gradient), two warm starts.

## Mechanism (named)

Free-boundary motion: ||du||_oo ~ 0.2 moves the src=0 kink across
nodes; s^6 responds O(1)-nonsmoothly, so exact dQ ~ 9 ~ lam_min while
linearized dQ ~ 4-5 (margins only x2-3 even linearized). Linear
transport needs ||e||_oo ≲ 0.03 (measured C-threshold, 0125 sweep).
Estimate-side repair EXHAUSTED (G1 already sharp — no further
looseness to harvest).

## Scope (what is and isn't closed)

- CLOSED: 0121 perturbation lemma via state-error transport from any
  p=6 bordered-Newton member at res ≳ 6e-3.
- OPEN (specified, not executed): (i) solver breakthrough to res ~1e-3
  (no candidate — Newton stalls, gradient floors, quad clean);
  (ii) witness not requiring pointwise state accuracy (discovery);
  (iii) lemma as CONDITIONAL with explicit residual hypothesis
  (shepherd scope call).
- BANKED: sharp G1 instrument replaces the 0121 bound as the tight
  linearized statement (5.21 trust / 4.16 fitted; margins x2.1/x2.6 —
  holds but thin, and transfer fails regardless); Q healthy at both
  states (11.16 / 10.77, doublet + axial preserved).

## Stop-rule execution
Charter gate fails (G1 4.2/5.2 vs 0.9) + G2' exact fails + S2 clean -> no-go certified.
No mesh rebuilds, no solves spent in 0124 (assemblies only).
G-a2 H-side inherits this wall (0121 deferral hardened, not lifted).
