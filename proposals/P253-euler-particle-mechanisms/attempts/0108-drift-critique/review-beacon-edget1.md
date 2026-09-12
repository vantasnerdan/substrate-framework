# drift firewall review — beacon EDGE-T1' (c1e6b14d): PASS + 1 fencing repair

Reran tilt/edge_t1prime.py from checkout (0.2 s): exit 0, numbers EXACT
(E0 1.973e-06; E1 peak 0.11818; E2 mean −8.708e-19; peak |dJ/J₀| 0.47274).
Provenance: mini-freeze 8c984414 → compute c1e6b14d ✓. Tracked. E0–E3
executed as frozen; E3 staged-not-fired; STOP exits preserved; FB-C1
untouched; static-tilt-only + 0160-out fences hold; MA-4 held (no static
tilt displaced, none claimed).

## J-firing genuine (no fit knob exists)

- Coefficient 4 is DERIVED (banked 0154 J(ε), read-only) — fitted to
  nothing here; no tilt response is computed at all. ν=−1 re-derived
  in-freeze from banked μ=K/10, λ=−K/15 (hand-verified: λ/2(λ+μ)=−1 ✓;
  κ_3D=λ+2μ/3=0 ✓; (1−2ν)/(1−ν)=3/2 ✓); η_in is the textbook Volterra
  form (hand-verified factor/structure ✓) with an FD guard at 2e-06.
  "Fires" = derived-order source genuinely nonzero (O(1) local signal),
  zero free parameters. Not fitted — unfittable.
- Monopole scoping honest: E2 pre-registered "mean vanishes → far-field
  monopole killed IN SCOPE, near-field stands"; measured −8.7e-19
  (machine zero); surprise-branch armed, untriggered. The kill is
  scoped exactly as frozen — near-field preserved, far-field dead. ✓

## R1 REQUIRED (fencing, framing only): edge ≠ ALIVE object

The numbers are clean; the FRAME overreaches. ALIVE was banked for the
SCREW sector (F2 + dipole + dynamics + persist). The edge defect has no
F2/dipole/dynamics/persist chain — its persistence, charge, mobility are
all unbanked. "First ALIVE-object medium coupling" transfers screw-ALIVE
to the edge object without a receipt. Repair (wording, receipts + any
downstream citation): the redirect SUCCESS stands as "first DERIVED-ORDER
defect→medium-sector coupling, sourced by the prescribed edge field";
edge-object aliveness explicitly UNCLAIMED (screw-ALIVE does not
transfer). No rerun; re-review = text check. P1-redirect success itself
is unaffected — the frozen branch asked for a derived-order source, and
got one.

## Verdict

EDGE-T1' **PASS** (conditional on the one-line fence). Preserved: a
derived — not fitted — defect→tilt-sector source with its far-field
overclaim pre-killed, and a lane that now owns a real medium coupling
without owning more than it measured.
