# Attempt 0017 — N2 pose arbitration: clean sympy rederivation confirms c1 = 1/2 - gamma; grid instruments failed with named mechanisms

## Route

Arbitrate the m=1 matching-pose split discovered in 0016 (Kelvin-(96)
two-condition pose c1 = 1/4 - gamma vs campaign four-condition pose
c1 = 1/2 - gamma) by two independent instruments: (a) a pose-independent
grid eigen-solve of the linearized-Euler vorticity system on the Rankine
profile; (b) a from-source sympy rederivation of the whole interface
system with explicit sheet displacement eta.

## Instrument (a): grid eigen-solve — FAILED twice, mechanisms named

- v1 (sparse shift-invert): ARPACK returned the shift sigma back as
  "eigenvalues" (every scanned w equaled its own seed; mode structures
  were garbage). All v1 numbers, including the apparent m=2 agreement at
  9e-4, are artifacts of seeding sigma at the prediction. VOID.
- v2 (dense standard form, exact K_m Robin exterior): the dense spectrum
  contains the discretized exterior shear continuum (a near-neutral
  ladder with O(1e-3) spacing) plus ka-independent artifact modes
  (e.g. w = -0.24200026 at every ka), inside which the physical branch
  cannot be identified without pose-committed mode matching. The
  arbiter property (pose independence) is therefore lost: VOID as
  evidence, preserved as instrument record.

Mechanism: the linearized exterior problem supports a dense grid
spectrum; identification requires either continuity tracking seeded in
an isolation window (none exists for the slow branch on this grid) or a
spectrally cleansed discretization. Route closed blocked-with-mechanism.

## Instrument (b): from-source sympy rederivation — DECISIVE

First sympy setup (single-Omega hand components) was itself defective:
it dropped the cylindrical basis-rotation terms of (v'.grad)v0.
Corrected lab-frame components (verified against the rotating-frame
Coriolis form):

    e_r : -I wt vr - 2 Om vth = -P'
    e_th: -I wt vth + 2 Om vr = -I m P/r
    e_z : -I wt vz = -I k P

give D = (w - m Om)^2 - (2 Om)^2 and lam^2 = k^2 (4 Om^2 - wt^2)/wt^2
(Poincare form, as in the 0015 verifier), with

    v_r  = I (Om m r P' + 2 Om m P - r w P')/(r D)
    v_th = (- Om m^2 P - 2 Om r P' + m w P)/(r D).

Exterior: potential phi = C K_m(kr), pi_out = I wt C K_m (axial
momentum). Four conditions with EXPLICIT eta (kin-in, kin-out,
pressure continuity, frozen sheet strength [v_th] = 2 Om eta), eta = 1
normalization, small-ka branch expansion wt = -Om - Om k^2 (lg + c1)/2
with lg = log(2/ka), K_1 carried to its Euler-gamma term:

    R_tang ~ k^2 * Om * (-2 c1 - 2 gamma + 1)/(denominator)  ->  c1 = 1/2 - gamma

sympy-verified (sp.simplify(sol - (1/2 - gamma)) == 0). This is now
derived from source, independent of every inherited 0008-0015
convention. Cross-check: the 0015 verifier's amplitude bookkeeping
(C = -I wt/dK, B-form, r1, r2) is in consistent eta=1 units and matches
the rederivation exactly — the verifier's 29/29 state is vindicated.

## Pose verdict (evidence scope)

- The campaign four-condition pose is internally consistent, source-
  derived, and exact at k = 0 (displaced-vortex identities re-verified:
  v'_th,in = -Om eta, v'_th,out = +Om eta, jump = 2 Om eta = frozen
  delta coefficient; single-valuedness across the displaced surface
  equivalent to the jump condition at the unperturbed radius).
- The restated Kelvin-(96) two-condition pose omits exactly the
  frozen-strength condition; its own roots carry c1 = 1/4 - gamma
  (robust bisection on the transcribed equation, three decades).
- The T-identity instrument (tangential-jump ratio of the (96)-pose
  field) remains inconclusive at finite ka because two independent
  convention slips poisoned successive hand transcriptions (factor-2 in
  v_th; lambda-doubling in the (86) ratio). It is NOT evidence either
  way and is recorded as an open construction.

C-CST-002 m=1 branch statement STANDS as recorded: bend channel
w* = -(Gamma k^2/4pi)(L + 1/2 - gamma) [NUMERIC_EVIDENCE, ASYMPTOTIC,
source-derived]. The split against (96)-as-restated is a documented
matching-pose difference (REPRESENTATION_SCOPED), with the frozen-
strength condition's k=0 exactness as the tiebreaker evidence.

## Defects caught this attempt (failure-generated candidates)

1. Hand component equations missing basis-rotation terms (single-Omega
   defect) — caught by the lambda^2 mismatch against Poincare; repaired
   by rederivation. Named mutation class for verify_cst002: M6
   "single-Omega Coriolis" (already present as M2) is joined by
   M7 "missing basis-rotation term" for any future operator build.
2. Series-order bug: residuals are O(k^2 lg); series cuts at O(t^0)
   return 0=0. Any future symbolic constant extraction must expand to
   the k^2 coefficient explicitly.
3. Instrument artifacts: shift-invert echo and continuum-polluted dense
   spectra — both recorded above; neither may be cited as evidence.

## Status

- Route verdict: established as stated for the arbitration route's
  algebraic arm; the numerical arm blocked with mechanism named.
- evidence_scope: REPRESENTATION_SCOPED on the m=1 pose split (unchanged
  from 0016, now with the split mechanism pinned to the tangential
  condition at linear order).
- Owed: N2 -> N3 handoff (C_tw recomputation under the J-species,
  ensemble moduli), claim-by-claim review, terminal campaign PR after
  the full ladder per the frozen obligations.
