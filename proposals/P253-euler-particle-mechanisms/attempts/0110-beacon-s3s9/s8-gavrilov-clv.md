# S8 Gavrilov + CLV verified (beacon 0110)

Sources: Gavrilov *A steady Euler flow with compact support*,
arXiv:1810.08020v1 — cached PDF re-hashed at 0110 boundary: `fcaca85f…e963`
MATCHES recorded hash (V1). CLV (Constantin–La–Vicol) *Remarks on a paper by
Gavrilov*, cims.nyu.edu/~vicol/CLV1.pdf — full PDF fetched this session
(V8; 529-line extraction, artifact://39).

## Exact consumed scope

- Gavrilov (per 0095 audit, hash re-verified): nonzero smooth compact
  finite-energy STEADY Euler field with meridional flow + swirl (Thm p.1,
  construction (6)–(7) pp.5–6). No stability/Hessian theorem.
- CLV Theorem 1: Gavrilov's `C_0^∞` compact steady result, reproduced via
  localizable Grad–Shafranov + cutoff `ũ = η(p)u` omitting the center point
  (pressure minimum ⇒ support = toroidal shell/annulus, NOT a ball).
- CLV load-bearing qualifications:
  - Swirl REQUIRED: smooth compact steady with `F ≡ 0` is identically zero
    ([10] cited). Every compact carrier in this method carries swirl.
  - Regularity: constructed velocity is HÖLDER continuous, Euler holds
    WEAKLY in the ball (Theorem 2); vorticity `L¹` with infinite derivative
    at the center (107). "Smooth compact" needs the cutoff-shell reading.
  - Theorem 3: Hölder multiscale steady states (Onsager-regularity class),
    knotted/linked Lagrangian trajectories possible — kinematics, not
    persistence.

## Perturbation class / norm

None — pure existence. The 0032/0038/0039 reviewed linear-semigroup lower
bound on one fixed Gavrilov carrier (adverse P2 evidence) is the only
stability-type result and points AWAY from that carrier.

## S2-carrier transfer verdict

- Strict compact velocity is a choice, not the particle definition
  (frozen-issue row confirmed: CLV support is a shell via cutoff; S2 needs
  only compact vorticity + decaying tails).
- No stability or interaction transfers; S8 does not touch the S2 family
  (different construction, swirl-carrying vs S2 no-swirl).
- Carrier-selection consequence: S2 (variational thin ring, no-swirl,
  reviewed linear observed growth) vs S8-comparator (compact, swirl,
  reviewed linear growth bound) — the campaign's two concrete carriers with
  opposite signs for P2. No new construction here.
