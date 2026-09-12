# Edge-T1' mini-freeze — FROZEN PRE-COMPUTE (beacon, 0157/tilt)

Status: FROZEN 2026-09-12. Charter: shepherd EDGE-T1' (T1-PASS #101,
P1-redirect provision of tilt/design.md). No compute under this doc
until acknowledged; freeze commit precedes any compute commit.

## Object

Volterra edge dislocation (line z, Burgers b.x̂) in the banked F-A
medium: compute the in-plane dilatation eta_in field, the derived
J(eps) modulation map dJ/J0 = 4.eta_in (eps_zz = 0 plane strain),
and stage the gate-3 readout (local pendulum-frequency shift).

## Derived observation (banked inputs only, pre-compute)

Banked F-A: W = K[tr(eps^2)/10 - (tr eps)^2/30] -> mu = K/10,
lam = -K/15 -> nu = lam/(2(lam+mu)) = -1; kappa_3D = lam+2mu/3 = 0.
Footnote travels: zero bulk modulus (the kappa-0 gap named in the
dipole freeze); tilt analysis uses eta_in only, bulk stability
untouched. Edge dilatation factor (1-2nu)/(1-nu) = 3/2: the banked
medium gives the edge defect STRONG dilatational content.
Volterra edge (isotropic, plane strain):
  u_x = (b/2pi)[th + x.y/(2(1-nu)r^2)],
  u_y = -(b/2pi)[(1-2nu)ln r/(2(1-nu)) + (x^2-y^2)/(4(1-nu)r^2)],
  eta_in = (eps_xx+eps_yy)/2 = -b(1-2nu).y/(4pi(1-nu)r^2) (dipolar,
  zero angular mean, ~1/r). eps_zz = 0.
MA-4 consistency (banked): NO linear-in-tilt source exists, so the
edge does NOT displace static tilt — the honest observable is the
STIFFNESS modulation dJ/J0 = 4.eta_in (local frequency shift via
gate-3 machinery, STAGED not fired). Any static-tilt displacement
in the numerics = implementation bug, receipted as guard E0.

## Frozen falsifiers + stop (scoped to edge)

- E0 (guard): numeric FD eta_in vs analytic on smooth bulk
  (cut-band + core + wrap-edge masks, T1 discipline); bar 1e-3.
  Breach = grid/cut bug, never a physics verdict.
- E1 (physics): eta_in == 0 within discretization anywhere in bulk
  -> contradicts Volterra structure at nu=-1; fires as
  IMPLEMENTATION kill (re-derive, no lane verdict).
- E2 (physics scope verdict): angular mean of dJ/J0 vanishes
  (dipolar 1/r, predicted) -> NO far-field tilt monopole: far-field
  tilt-charge claim KILLED IN SCOPE; near-field modulation stands
  as the staged observable. A nonzero monopole would contradict
  the analytic dipole structure (guard-grade surprise, disclose).
- E3 (back-reaction staging only): core-shift estimate vs B1 bars
  recorded, NOT fired here (belongs to T2-edge charter).
- STOP: edge-T1' banks the source verdict + E2 scope verdict; then
  EITHER T2-edge charter (back-reaction, needs its own freeze) OR
  O(eps^2)-constant charter OR dipole-dipole runner-up per STOP
  rule. No new formalism: Volterra + banked F-A W + banked J(eps)
  only. FB-C1 never fires; static tilt only; 0160 out (inherited).
