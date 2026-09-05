# Attempt 0023 — N6 no-Cosserat contrast side established: verify_cst006 7/7

## Route

Orientation-ergodic contrast ensemble (obligation N6): exact isotropy
argument (parts/n6_part1_ergodic.py) + simulated decorrelated ensemble
confirmation (parts/n6_part2_simulation.py), consolidated in
`verify_cst006.py` (7 checks, 2 mutations, exit 0).

## Result

Decorrelation premise: <n_i t_j> = <n_i><t_j> = 0 (independent isotropy of
tangent and frame distributions). Consequences, all exact:

- phase-averaged frame map <L(phi)> = 0 over the uniform phase circle;
- the couple moduli (c_tr, c_s, c_a) and the spin stiffness are identically
  zero -- the decorrelated ensemble has NO kappa dependence in e_seg: no
  couple operator, no spin stiffness;
- the effective medium reduces to the Navier-Cauchy sector (lambda, mu) with
  the SAME stretch coefficients as N3 (decorrelation does not touch the
  tangent distribution).

Contrast recorded: the locked ensemble (N3) carries the nonzero couple
operator (c_s = B L_v/10, c_a = B L_v/6, c_tr = -B L_v/30, alpha = L_v T/6).
Under full tangent isotropy the NET couple stress vanishes in BOTH ensembles;
the N6 contrast is the couple OPERATOR itself (response to applied wryness),
not the net mean.

## Simulation confirmation (declared observable)

N = 200000 decorrelated segments (seed 20260904), quadrupole probe
w = n_1^2 - n_2^2 against the coherent frame projector L(phi):

    m_net[1] = -7.88e-4   tol(5 sigma/sqrt(N)) = 4.09e-3   within: True
    m_net[2] = -1.54e-3   tol(5 sigma/sqrt(N)) = 4.08e-3   within: True
    scale-relative |m_net|/scale = 1.5e-3 << 1.1e-2 declared

## Checks and mutations

- <L(phi)> = 0 (sympy integration); locked counterpart L = I (contrast real).
- Moment factorization premise encoded and checked.
- Simulation within declared tolerance (both components).
- M1 biased sample (7 sigma/sqrt(N) offset) flagged by the tolerance check
  (verifier sensitivity). M2 correlated-frames joint moment (delta/3) rejected
  against the product measure.

## In-run defect

- First M1 bias (5.1 sigma/sqrt(N)) was inside the recomputed tolerance
  margin (the pre-existing mean offset ate it); raised to 7 sigma/sqrt(N).

## Status

- route_verdict: established (contrast side exact + simulated confirmation
  within declared tolerance)
- evidence_scope: EXACT (isotropy argument) + SIMULATION_EVIDENCE (MC within
  declared tolerance)
- Next unlocked: N7 EPS bridge (declared existence layer; Beltrami realized
  example vs declared ensemble premises).
