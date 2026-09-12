# Receipt: S1 κ-fit (FROZEN F1 bar) — verdict KILL with named mechanism

Command: `python3 run_kappa.py gate` (PASS 2.6e-4 at eps=2.5e-5) then
`python3 run_kappa.py fit`. Env: CPython 3.12.2, numpy, reduced axisym model
(Saffman self + RM mutual, NQ=200, dt=0.004 — PoC-2 gold backbone).
F1 metric/bar applied exactly as frozen in ../../03-prefire.md.

## Numbers
- Passage A (banked section): κ=2.145, r=0.984.
- Passage B (T/4 later, R≈1.001/1.000): κ=2.041, r=0.753.
- Ratio 1.05 (slopes agree!) but residuals ≫ 0.25 → F1 verdict: KILL.
- Diagnostic (not a re-verdict): residual decomposition — Lorentz-template
  3%, phase-drift 18%, rest 78% (elliptic section-mode mixing). Not a soft
  artifact; kill stands on substance.

## Kill mechanism (named): poloidal-vs-toroidal plane mismatch
In axisymmetry the Magnus-form perturbation (t̂×W, t̂ azimuthal) lives in the
poloidal (R,Z) plane, while a genuine Lorentz-form force v×B (v,B both
poloidal) points AZIMUTHAL (toroidal) — orthogonal planes by symmetry. The
integrated axisymmetric response therefore CANNOT take Lorentz form at any κ;
slope agreement (1.05) is scale coincidence. S1-as-stated is dead.
Design correction logged (not executed): perturbation used mutual-slip per
sketch; template used mutual-field Lorentz form; correct Lorentz template
needs full-v × B (toroidal response) — unmeasurable in axisymmetric
diagnostics by construction.

## Constructive residue (charter separately, real build)
S1's honest successor is 3D non-axisymmetric: Magnus response in m≥1 channels
(A3 3D code + Magnus term; F1 bar generalizes with per-m templates). Not built
here. S1 axisymmetric route: CLOSED (killed, mechanism named).
Gate bug found+fixed in-run: self-field dominated first perturbation
(self comoves — mutual-slip per design); eps lowered 1e-4→2.5e-5 for gate.
Lint: style-only. No imports; verdict per frozen bar, no judgment calls.
