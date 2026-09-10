# 0114 drift-firewall repairs (beacon, append-only)

Firewall: drift 0114 PASS (71667472 verified), 4 non-fatal hygiene repairs.
All landed below; no verdict, criterion, or number altered.

1. Mass-gate relabel (vacuous frozen-scope): design.md Repair R1 —
   amended mass gate = pipeline-integrity, vacuous physically; exposing
   content is slope+R² only. Frozen + amended verdicts unchanged.
2. Seed-7 provenance (argv+archive): `s9_probe.py --seed` argv added
   (default 0); `seed0-run.log` (exit 1, frozen BLIND, numbers reproduce
   run-1 exactly: 0.1262/0.9996/0.0022) + `seed7-run.log` (exit 1 frozen;
   numbers 0.1262/0.9996/0.0008 feed the amended reading) archived.
3. Dead control set removal: rc/zc set advected but never read — deleted
   (2 lines + comment fix). Seed-0 rerun reproduces run-1 numbers exactly,
   proving removal changed nothing.
4. Pycache cleanup: `__pycache__/` removed under 0114 (+ identical 0111
   litter in the same pass); verified zero `*.pyc` remain.
