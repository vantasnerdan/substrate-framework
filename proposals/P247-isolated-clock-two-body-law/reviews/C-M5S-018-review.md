# Review: C-M5S-018 (width quasi-moduli obstruction)

- Reviewer of record: independent reviewer agent Review018 (bounded to this claim's frozen transaction)
- Transaction: six new proposed registry entries C-M5S-015..020; no existing claim modified.
- Verdict: REQUEST_CHANGES (evidence materialization and measurement scope).
- Blocking findings (all repaired in the same transaction):
  1. The cited attempts/0009/width-mode.json did not exist (the repaired mode characterization lived only in result prose). REPAIR: materialized by attempts/0011/materialize_width_mode.py with the repaired-class bottom-mode spectrum, row fractions (m2 fraction 0.99999994), drift overlaps (0.9059657 vs R24->R30), and the family core-width ladder over the actual roots R12..R30.
  2. Core-width clause measured on reprojected seeds, not minimum-branch roots. REPAIR: statement now quotes the measured family ladder (1.775 at R=20 to 2.074 at R=30, +17 percent against +50 percent box growth).
  3. "Accelerating increments" false as unqualified (first increment decreases). REPAIR: restated as monotone growth with increments increasing from R=22 onward (0.739, 0.724, 0.737, 0.777, 0.847).
- Reviewer-verified facts: width-stiffness-vs-R.json exactly contains the three quoted stiffnesses, zero negative directions, and the R=24 overlap 0.8645158 (vs R24->R26 drift); energy ladder traces through C-M5S-015 to attempt 0003, transferred to the repaired class by C-M5S-017/0009 (exact chi=0 S-energy preservation).
- Statement updated to the materialized dual-overlap form (0.8645158 vs R24->R26; 0.9059657 vs R24->R30).
