# B4-C1 F-bar freeze (pre-compute; D-08 form) + build caps (frozen)

Build caps (from survey, now frozen): paper + sympy ONLY — no fluid compute,
no ensembles, no solves; single computation (S² Hamiltonian cocycle battery
below); wall < 1 h; background priority (yields to 0131 lane needs).
Guillotine stands: R-valued data self-compactifying revives nothing — a
nonzero result must exhibit a GROUP/manifold source (here: SO(3) control vs
non-so(3) pairing), else it dies with the rest.
Construction: ham(S²) with area form; Lichnerowicz-type cocycle
c_h(f,g) = ∫_{S²} h·{f,g}·ω, h = cosθ (non-constant ⇒ nontrivial candidate;
constant-h gives identically zero — first tripwire). Tests on exact trig
harmonics (sympy, exact integrals):
- T0 (cocycle identity): δc_h = 0 on random triples (framework gate; FAIL =
  bug, stop, no verdict).
- T1 (Whitehead control): c_h|_{so(3)} must be a COBOUNDARY (H²(so(3))=0):
  solve c_h = δλ on ℓ=1 basis. FAIL (unsolvable) = framework bug, stop.
- T2 (LIVE): c_h on non-so(3) pair (Y_20, Y_21^c): nonzero ⇒ nontrivial
  class ⇒ candidate integral level (report pairing number; carrier ID open).
- Verdict: PASS-lean (C1 first gate) iff T0+T1 pass AND T2 ≠ 0 (tolerance:
  |T2| > 1e-6 exact-arithmetic nonzero — no D-08 margin needed, exact math);
  KILL (this cocycle class) iff T2 = 0 exactly (extension undetected where
  it could live); gray iff T1 fails (re-examine, no verdict).
  PASS-lean ≠ promote: carrier identification (which orbit? which level sets
  charge?) remains open → C1 advances to carrier gate, nothing is claimed.
Run reports (T0/T1/T2 + pairing number) + verdict, nothing else.
