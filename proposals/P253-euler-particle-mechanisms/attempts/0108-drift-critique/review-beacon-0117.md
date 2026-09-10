# drift firewall review — beacon 0117 member build (d8508f07, NOT-DONE)

Transaction: design.md (frozen) + build_member.py + feed_member.py + 3 npz + rung-log (R0–R13)
+ tool-receipts (T1–T3). Claimed: STAGE-1 LOCATED, G-a2 NOT DONE (blocked-with-mechanism).
Drift re-ran the feed (not the multi-hundred-second build) and checked every cited mechanism.

## Convergence + iz≈π typing — PASS (honest post-hoc)

0.735→0.020 in 3 Newton iters = bordered formulation working, then systematic floor (ray
autopsy: true stall, ascent direction) — "Formulation VALIDATED" reads at formulation scope,
immediately qualified by floor + NOT-DONE. iz=3.23≈π is design monitor (a) firing (P1 jet
I_z=πρκR² at frozen inputs κ=R=ρ=1), NOT fitted (κ,R frozen inputs; iz output); 2.8% deviation
consistent with coarse mesh + 2% floor. Post-hoc label honest.

## Fixes + FD Jacobian — PASS (mechanisms sound, verified)

ddot-vs-dot (400× stiffness vs premade-laplace convention), sign-flip (FD-proved 200%, fixed to
+M3@jf, 8e-7/3e-8), row-degeneracy (κ/iz rows parallel over thin core at R=1 — physically
correct mechanism; 0080 BR-border fix follows), line-search apply-bug vs true stall
distinguished, J-audit 1e7 false alarm correctly diagnosed as a wrong check (J@d vs quotient).
FD agreement + 1.8e-14 solves + 3.5e-15 load-path agreement are real verification. Rung log
append-only with wrong-attractor Picard collapses preserved. No bug hidden.

## Exploratory feed — numbers REPRODUCED, four repairs (non-fatal)

Drift reproduced λ_ω=4.7172 doublet + 13.35 axial (floor 2.1e-1) via tested leray/make_grid
reuse, H-only with G-await-Maxwell, APPROXIMATE labels intact. Repairs: (1) mesh-args
provenance gap — feed DEFAULTS (80/40) fail `u.size==basis.N` assert; npz needs nr=40/nz=20:
record exact CLI + mesh params in npz, else the receipt is not self-reproducing; (2) ~20%
systematic needs ITEMIZATION per design §Tolerances (member floor 0.21 + interp + smoothing +
box); (3) feed zeta uses δ=1e-3-smoothed positive part (softplus) vs the build's sharp (·)₊^p —
UNLABELED near-free-boundary deviation: label + bound it or use sharp max; (4) design monitors
(b) dipole-match (c) observed-order (d) jitter have NO status lines anywhere — report pending
or polish-rung assignment (kappa_hat monitor (a) is live in-code; the rest silent).

## Polish-rung block — LEGITIMATE

Fine-mesh stalls (0.07–0.21, identical ray signature) + p5 timeout correctly labeled
non-verdict; trust-region from warm start is the appropriate next step for a true stall;
fine→Maxwell→production order correct. Mechanisms banked, no effort-closes-G-a2 claim.

## Verdict

PASS the NOT-DONE build and exploratory feed (repairs 1–4 open, none verdict-changing; pycache
nit recurs — third appearance). G-a2 stays BLOCKED at member-convergence + Maxwell stage; no
verdict narrowed or widened. Next executable: trust-region rung with monitor-status lines.
