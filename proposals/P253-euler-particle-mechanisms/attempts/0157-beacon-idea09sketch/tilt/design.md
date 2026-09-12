# Tilt-coupling design — FROZEN PRE-COMPUTE (beacon, 0157/tilt)

Status: FROZEN 2026-09-12. Charter: shepherd DECISION (owner-delegated),
F-C tilt-channel defect coupling ACCEPTED. No compute executes under
this doc until drift/shepheard acknowledge the freeze; one run, no
iteration on physics unless chartered (standing discipline).

## Object

Does the banked screw defect's shear field source the F-C tilt
sector via the DERIVED J(eps), and does tilt back-react on the
defect (core energy, mobility, Peierls-type barrier)?

## Banked inputs (read-only, never edited)

- 0154/04-fc-build.md: J0 = G^2.pi.R^4/(2.l^3) > 0;
  J(eps) = J0.(1 + 4.eta_in - 3.eps_zz); omega^2(k) = (4.J0/I).sin^2(k.l/2);
  U_c(th) = n.G.Gc.Om(th)/4pi; Om2 = -pi.R^2.d.(d^2-2R^2)/(2.(R^2+d^2)^{5/2});
  d(om^2) = n.G.Gc.Om2/(2.pi.I); gates 1-3 measured; FB-C1 pre-registered;
  footnotes travel (I scaling-level, quadrupole untracked, tau priced).
- 0157/dipole-run: screw dipole strain, only eps_xz, eps_yz nonzero,
  tr eps = 0 exactly; units a=1, b0=1, mu=1.
- 0157/persist: core disk r<1 = declared sub-continuum regularization.
- 0157/dynamics: D3 Peach-Koehler glide estimate (mobility baseline).

## Freeze-level observation (shapes T1, not a result)

The DERIVED J(eps) sees in-plane dilatation eta_in and axial
strain eps_zz. Screw-defect strain is pure off-diagonal shear:
eta_in = 0, eps_zz = 0 at O(eps). EXPECTED T1 OUTCOME: null at
derived order. T1 therefore surveys, in fixed order: (a) J(eps_screw)
evaluated exactly (null expected — receipted, not assumed);
(b) O(eps^2)/O(th^2.eps) tilt-strain couplings (existence survey
inside F-C's stated remainder orders — no new formalism);
(c) edge/dilatational-defect fallback SPEC (has eta_in/eps_zz at
O(eps)): if (a)+(b) are silent, the P-kill fires FOR SCREW-IN-F-C
and the charter REDIRECTS to edge-defect T1' — lane survives,
candidate narrows. This branching is frozen here, pre-compute.

## Rounds (bankable intermediates, in order)

- T1 SOURCE: (a) J(eps) on the banked screw field, numeric +
  analytic (expect 0 to discretization); (b) remainder-order survey;
  (c) fallback spec if silent. Banks: source-term verdict.
- T2 BACK-REACTION: static tilt background -> defect core-energy
  shift (vs dipole B1 bars) + mobility shift (vs D3 baseline).
  Banks: back-reaction numbers.
- T3 READOUT STAGING: defect-induced static tilt -> d(om^2) via
  banked gate-3 machinery (STAGED ONLY — FB-C1 NOT fired, no
  measurement claim).

## Frozen falsifiers + stop

- P1: no tilt sourced at O(J-derived) AND none in remainder survey
  -> SCREW-IN-F-C KILLED (stated scope); redirect to edge T1'.
  If edge also silent at O(eps): tilt-channel defect coupling
  KILLED, dipole-dipole runner-up activates.
- P2: tilt back-reaction unbinds core (E_core shift exceeds B1
  refinement bars 0.07-0.08%) or flips mobility sign -> back-reaction
  channel KILLED; source result stands alone.
- STOP: two rounds banked, or P1/P2 fires, or derivation needs
  formalism outside banked F-C remainders (new charter required —
  no silent extension). FB-C1 never fires under this doc.
- FENCES: static tilt sector only (F-C as built); dynamical director
  (0160, D3/D4 owed) explicitly out; no electron/carrier/measurement
  claim (ALIVE fence inherits); additivity assumption (RC7 labeled)
  travels into T3 staging.
