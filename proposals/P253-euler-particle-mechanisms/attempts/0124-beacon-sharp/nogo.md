# 0124 NO-GO certificate: state-error-transport lemma route (beacon)

Charter: shepherd ruling on 0123. Instrument: sharp_dQ.py (Weyl-exact
3x3 dQ + exact nonlinear transfer + P2 quadrature probe).

## Measured (all exit 0)

- Fitted state: G1 = 0.69 (PASS <= 0.9, margin x15.5 vs lam 10.77);
  G2'-exact central = 9.10 (FAIL).
- Trust-r3 (best, res 6.6e-3): G1 = 0.87 (PASS, margin x12.9 vs lam
  11.16 — instrument reproduces feed-trust-r3.log digit-exact);
  G2'-exact central = 9.18, half = 3.49 (FAIL).
- S2 quad probe: P2-R rms 8.9e-3 vs P1-R 10.3e-3 (x1.16 — clean,
  quadrature is NOT the floor).
- Floor robustness: res ~1e-2 across 3 meshes (coarse/fitted x2),
  P1/P2, 140+ iters (Newton + gradient), two warm starts.

## Mechanism (named)

Free-boundary motion: ||du||_oo ~ 0.2 moves the src=0 kink across
nodes; s^6 responds O(1)-nonsmoothly, so exact dQ ~ 9 ~ lam_min while
linearized dQ ~ 0.7. Linear transport needs ||du||_oo ≲ 0.02 (x10
below the robust floor). Estimate-side repair EXHAUSTED (G1 already
sharp — no further looseness to harvest).

## Scope (what is and isn't closed)

- CLOSED: 0121 perturbation lemma via state-error transport from any
  p=6 bordered-Newton member at res ≳ 6e-3.
- OPEN (specified, not executed): (i) solver breakthrough to res ~1e-3
  (no candidate — Newton stalls, gradient floors, quad clean);
  (ii) witness not requiring pointwise state accuracy (discovery);
  (iii) lemma as CONDITIONAL with explicit residual hypothesis
  (shepherd scope call).
- BANKED: sharp G1 instrument replaces the 0121 bound as the tight
  linearized statement (0.87 trust / 0.69 fitted); Q healthy at both
  states (11.16 / 10.77, doublet + axial preserved).

## Stop-rule execution

G1 pass + G2' fail -> S2 probe (done, clean) -> no-go certified.
No mesh rebuilds, no solves spent in 0124 (assemblies only).
G-a2 H-side inherits this wall (0121 deferral hardened, not lifted).
