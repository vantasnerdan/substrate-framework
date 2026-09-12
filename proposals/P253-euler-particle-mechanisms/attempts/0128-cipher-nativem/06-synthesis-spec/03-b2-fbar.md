# B2 F-bar freeze (pre-compute; D-08 form) — dynamical charge-observable

Why helices (not rings): S1 proved axisymmetric response cannot take Lorentz
form at any κ (poloidal-vs-toroidal trap). The observable MUST live in 3D:
helical filaments (framing-distinct L) whose response CAN be toroidal.
Construction: single helical ring (R=1, helix amplitude a_h=0.1, L turns) in
A3 3D code (N=64, Saffman-local + RM mutual); Magnus test-force
a_M = ε·(t̂×W) with W = full local fluid velocity (3D vectors); short window
T_win=0.5 (no periodic orbit needed); template T = v×B full-3D Lorentz form
on BARE kinematics (v = filament velocity, B = BS field); slope κ_L =
⟨R,T⟩/⟨T,T⟩ per L ∈ {1,2,3}; ratios vs 2/1, 3/1, 3/2.
- PASS (missing-5 SUPPLIED): all three ratios within 2× of integer ratios AND
  all residuals ≤ 0.25 → dynamical observable exists; P1/F-a upgrade toward
  executable (pending firewall).
- KILL-observable (absent): all residuals > 0.50 → form wrong even in 3D;
  S1-3D direction CLOSED with mechanism (state which channel carries it).
- Gray otherwise → UNRESOLVED + named leg (longer window, N, eps).
- Gates (pre-verdict): eps-linearity (halving, rel <1e-3); N-leg 32 vs 64
  (slopes within 20% or flag N-dependence, no verdict).
Honesty scope: slopes test TRANSFER of single-parameter form across framing
integers (template-geometry coupling caveat recorded in receipt, not in bar).
Run reports (κ_1,κ_2,κ_3, ratios, residuals, legs) + verdict, nothing else.
