# Receipt: A3 resonance scan (exploratory PoC scope, model-capped)

Command chain (all from `receipts/a3-scan/`, env CPython 3.12.2 numpy 1.26.4):
`python3 run_a3.py gate` / `orbit` / `mono [aa [N]]` / `m0 [aa [N]]` /
`newton3` / `verdicts`; plus `python3 <repo>/0108-cipher-radical/receipts/poc2-filament/run_poc2.py`
verbatim rerun. Design: `0120-cipher-m2b1/A3-design.md` (C1–C8 + D4).
Lint: pycodestyle long-line warnings only (receipt precedent: bytes frozen).

## Banked numbers
- gate: radial 0.090663083528 vs PoC-2 mutual, azim 3.5e-17 — PASS.
- Shot orbit (smooth 3x3 single shooting, quadratic 1.9e-3→1.9e-5→2.9e-12):
  N=200: (0.777140, 1.188978) T=4.08800 res 2.9e-12; N=64: (0.77714040, 1.18897773)
  T=4.088 res 2.9e-12 (same orbit class; dt=0.004 converges to banked
  (0.773723, 1.185226) exactly — quadrature/dt shift only).
- PoC-2 verbatim rerun reproduces receipt digit-for-digit (Newton 6.4e-11,
  eigs 0.92751474±0.37378072i, res2 6.42e-11).
- m≥1 (shot orbit): all |ρ|−1 ≤ 3e-6 (m=1..6); eps-leg (halving, m=1,2):
  4.3e-9/7.9e-9; T-window leg (±1 step, m=1,2, N=64): max 1.000008 — ROBUST.
  SOFT3-deflated m=1: n_grow(1e-4)=0. Verdict: PASS-in-model (stable).
- W_MEASURED=0.76 (FD velocity-gradient operator norm × R̄²/Γ; definitional
  cousin — not identity — of Ruban's straining-flow W); Λ_Saff=5.05, Λ_ln=2.98.
- m=0 WINDOW FRAGILITY (main finding): PoC-2 flow_shape monodromy at fixed
  step counts — 1021 steps: {0.9275±0.3738i, 1} (banked PASS); 1022 steps:
  {1.0±0.006i, 1}; 1023 steps: {1.46, 1.0, 0.685} (saddle). ONE RK4 step
  (0.1% of period) swings elliptic→marginal→unstable. Banked PoC-2 elliptic
  PASS is therefore window-luck (int(T/dt) truncation of fp-dusty T), NOT a
  stability result. DOWNGRADE FILED: PoC-2 Floquet verdict → UNRESOLVED by
  fixed-T Floquet (receipt itself stands as recorded computation; its
  verdict does not). Needs phase-free (section-map) stability or analytic
  framing. m=0 A3 verdict: UNRESOLVED (method-limited), not stable, not unstable.

## Verdicts vs frozen predicates
- Ruban weak-m=1,2-growth test: NO licensed growth anywhere (|ρ|−1 ≤ 8e-6 ≪
  1e-4 floor) → prediction neither confirmed nor excluded (may live below
  floor). O(1) instability excluded in all resolved channels.
- a-ladder / Buttà-control / Rankine-variant legs: NOT RUN — reframed as
  premature until a window-robust m=0 diagnostic exists (ladder on a fragile
  diagnostic measures windows, not bands). Proposed next: section-based m=0.

## Method lessons (reuse)
1. Fixed-T Floquet on drifting-pair orbits is window-fragile — always run the
   ±1-step leg before consuming any multiplier verdict (would have caught this
   at PoC-2).
2. m=0 FD must use shape coords (axisym leak made an all-zero block); GAUGE2
   (azimuthal per-ring) is exact kernel — name, don't compute.
3. Crossing-detect Newton branch-jumps under orbit change — smooth single
   shooting (unknowns incl. T) converges quadratically where section-Newton
   slides into valleys.
4. Eval kernels cache modules — fresh interpreters for changed code; background
   via harness auto-background (foreground long call), never nested `(…) &`.
5. Stray CWD-relative saves (Mono npz at repo root) — moved in receipt dir.

## OUT (unchanged scope)
Live-field Euler check (filament law ⇒ continuum) uncrossed; all verdicts
-in-model. a-ladder/Buttà-control/Rankine legs queued behind window-robust m=0.

## Landing (N=128, shot orbit, deflated verdicts)
m=1..6 deflated|rho| all 1.000000 (6dp), n_grow(1e-4)=0 everywhere;
m=1 soft-attrib=4 (x/y-translate pairs, as designed); eps-leg 3.15e-9/3.77e-9.
Two resolutions (N=64/128) + T-window (±1 step) + eps-halving all agree:
m≥1 PASS-in-model (stable) is the most over-determined verdict of the campaign.
m=0 stays UNRESOLVED (window-fragile); PoC-2 downgrade stands as filed.

## R1–R5 closure (drift re-verdict)
- R1: Mono_N128_e1e-06_m06.npz (N-tagged copy of landing npz; |ρ|−1 ~1e-10–1e-9).
- R2: window1021-1023.npz + numbers above (downgrade evidence banked as artifact).
- R3: no run.log exists — compute ran under harness auto-background with output
  delivered to transcript (captured verbatim in the numbers above); empty
  dead-launch logs were committed once by accident and removed the next commit.
  Rerun commands reproduce every number (see header).
- R4: D5 reframe in A3-design.md.
- R5: cross-pointer appended at the 0108 PoC-2 receipt (see that file).

## Section-m0 landing (D6r2 window-exact flow_frac + fractional-T legs, N=200)
Shot orbit re-converged by window-exact shooting to BANKED values
(0.773723, 1.185227) T=4.08803 res 6.4e-11 — the 0.4% shift was itself
truncation artifact. Krein arc (section-flow 3x3, Teff exact):
T−2: 0.9547±0.2976i (|.|=1) · T−1: 0.9773±0.2117i (|.|=1) ·
T+0: 1.000183/0.999817/1.0 · T+1: 1.2367/1.0/0.8086 ·
T+2: 1.3498/1.0/0.7409. Elliptic pair collides at 1 then splits reciprocal-real
across ~0.1% of period: the orbit sits AT a Hamiltonian-Hopf point in m=0.
eps-leg at T+0: 1e-6/5e-7 agree; 2.5e-7 roundoff-fogged (flow roundoff/eps ~ 4).
N=64 own-shot-T gives {1,1,1} flat — cross-resolution spread ±2e-4 VETOES any
growth reading of N=200's 1.000183 (D-08: floor necessary, not sufficient).
m=0 verdict: MARGINAL-COLLISION (mechanism mapped both sides; exact-period pair
within 2e-4 of unit circle; growth UNLICENSED). Banked elliptic PASS mechanism
explained: 1021-step window lands on the arc's elliptic side.
