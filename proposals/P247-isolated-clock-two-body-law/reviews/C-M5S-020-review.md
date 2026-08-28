# Review: C-M5S-020 (de-boxed attractive exponential law)

- Reviewer of record: independent reviewer agent Review020 (bounded to this claim's frozen transaction)
- Transaction: six new proposed registry entries C-M5S-015..020; no existing claim modified.
- Verdict: REQUEST_CHANGES (evidence materialization).
- Blocking finding (repaired in the same transaction): the refinement clause (gamma = 0.2524/0.2571, E-infinity stable to 0.1 percent) existed only as prose; v-law-measurement.json contains the primary R=24 ladder fit only. REPAIR: refinement re-materialized by attempts/0007/v_law_refinement.py into attempts/0007/v-law-refinement.json (2x grid: gamma = 0.271731, E_inf = 8.08826, ssr 4.75e-4; 1.5x box: gamma = 0.265034, E_inf = 8.095008); the claim now quotes these materialized values (exponent shift < 12 percent, E-infinity shift < 0.4 percent) with the artifact and script as evidence. The result.yaml prose values (0.2524/0.2571) are superseded by the materialized record.
- Non-blocking note (applied): v_law.py's docstring said V(P) = -0.5 while the implementation and JSON give zero; fixed.
- Reviewer-verified facts: exact potential form, 8.5e-13 integral match, projector vacuum, primary fit, no-Yukawa interpretation, C-M5S-008 contrast, and dependency uses consistent with the records.
