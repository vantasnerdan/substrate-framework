# Independent review: C-SIG-001, C-SIG-002, C-SIG-003 (P225)

Reviewer: independent subagent (proposer's modules not trusted; S recomputed from
the bare formula in `reviews/recompute_S.py` with no repo imports).

## Verdicts

### C-SIG-001 — induced inverse Newton coupling from the sG mass tower — ACCEPT
Independent recomputation of the sum over kink (mass 8, mult 2) plus breather
levels E_n = 16 sin(n/16), n = 1..25 (n·h < 8π), at h = 1, ξ = 0, with cutoff
Λ = 8π and per-species term (1/6−ξ)m² log(Λ²/m²)/(2π) gives
S = 120.058077992326 (claimed 120.058077992; agreement to 3.3e-10, i.e. the
rounding of the quoted value) and q = 1/S = 0.00832930209048 exactly as claimed.
The per-species coefficient in `induced_inverse_g_species` matches the standard
one-loop a2/Einstein–Hilbert matching: the heat-kernel R-coefficient
(1/2)(1/6−ξ)m² log(Λ²/m²)/(4π)² equals 1/(16πG), which rearranges to
1/G ⊃ (1/6−ξ)m² log/(2π). Confirmed. The structural premises (ξ = 0, kink
multiplicity 2, breathers-as-loop-species, phonon = n=1 level, tower-edge
self-cutoff Λ = 8π) are declared premises rather than derived facts, but they
are exposed for mutation testing and the verifier mutates them, so this is
honest conditional structure, not a hidden fit.

### C-SIG-002 — static weak-field with derived G — ACCEPT
The verifier confirms the Schwarzschild exterior factor f = 1 − 2m/x, the
Newtonian potential Φ = −GM/r with the derived G (written M/(S·ħ·r) in
normalized units), and the standard deflection 4GM/(c²b). The C-GRV-001
monomial q·ℓ²·c³/ħ has dimensions L²·L³T⁻³/(ML²T⁻¹) = L³M⁻¹T⁻² = [G],
independently checked. Nothing here is fitted; it is pure substitution of the
derived G into accepted C-STG-002/C-GW-001 forms.

### C-SIG-003 — breather quadrupole power — ACCEPT
Dimension check in (M, L, T): [G/c⁵] = (L³M⁻¹T⁻²)/(L⁵T⁻⁵) = M⁻¹L⁻²T³;
[(E_scale/c²)ℓ²ω₀³] = M·L²·T⁻³, squared = M²L⁴T⁻⁶; product = ML²T⁻³ = [power].
Correct. The numerical bracket 379.46463806875 is the midpoint of the accepted
C-SG-010 interval (379.4646380687, 379.4646380688) in `governance/claims.yaml`
— a legitimately imported accepted-claim value with stated double-precision
bounds, not a new fit. Power is positive, scales linearly with G, and the
verifier's mutation tests show sensitivity.

## Test/verifier runs
- `pytest tests/test_sg_spectral_induction.py tests/test_sg_gravity_chain.py -q` → **24 passed** (2.40 s).
- `proposals/P225-sg-spectral-induced-gravity/verify.py` → **ALL 21 CHECKS PASS, exit 0**.
- Test suite includes genuine mutation tests (conformal-coupling annihilation,
  kink-multiplicity shift, action-quantum level removal, onsite-coefficient
  power scaling) — tests can actually fail.

## Defects
1. **(low, performance)** `verify.py` is slow: it exceeded a 180 s foreground
   timeout and completed only after ~6–7 minutes wall clock (sympy exact
   arithmetic plus mutation recomputation). Not a correctness issue, but it
   makes replay expensive.
2. **(low, provenance)** The C-SG-010 bracket bounds are hard-coded as
   `sp.Float` literals in `sg_gravity_chain.py` (lines 50–51) rather than
   imported from an importable C-SG-010 artifact; the literal lives only as
   prose in claims.yaml. Values match the accepted interval, so no
   discrepancy — but a future C-SG-010 revision would not propagate
   automatically.
3. **(informational)** The n_max computation uses a floor with a 1e-12 epsilon
   (`maximum_lattice_level`); harmless at h = 1 (25 < 8π ≈ 25.13 is far from
   the boundary) but fragile for h values that land near-integer multiples of
   the edge.

No defect blocks acceptance. Overall: **all three claims ACCEPT**.
