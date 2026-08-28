# Review: C-M5S-015 (milestone-0 de-boxing adjudication)

- Reviewer of record: independent reviewer agent Review015 (bounded to this claim's frozen transaction)
- Transaction: six new proposed registry entries C-M5S-015..020; no existing claim modified.
- Verdict: REQUEST_CHANGES (wording and dependency repair; no scientific defect).
- Blocking findings (all repaired in the same transaction):
  1. Free-wall scope overextension: free-wall.json contains only the order-16 R=8 endpoint; the no-limit verdict was scoped as if a free-wall ladder existed. REPAIR: scope narrowed to pinned-wall radial minimum ladders; R=8 labeled a pinning-cost endpoint.
  2. Band misattribution: [19.20, 29.96] is the order-24 self-consistent band over R={14,18,24}, not the R=24-root window ([19.1991, 60.0]). REPAIR: attribution fixed in the statement.
  3. Missing dependency: the confined-stability premise uses C-M5S-009; the dependency edge was absent. REPAIR: C-M5S-009 added.
- Non-blocking notes (applied): attempts/0003/result.yaml added to evidence (the +/-20 percent gate lives there); scope-difference note added for C-M5S-012's fixed-J exclusion vs this claim's variational energy; energy-growth increments qualified (monotone; increasing from R=22 onward: 0.739, 0.724, 0.737, 0.777, 0.847).
- All other quoted values verified against artifacts by the reviewer.
