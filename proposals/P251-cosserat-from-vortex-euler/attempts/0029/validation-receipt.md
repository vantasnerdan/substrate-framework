# Validation receipt: PR199 continuation

Base: 10a0f31bf0122f4c744a81dfe1cfa66506d57ad7.
Validated staged tree before this record: 0f1779eb56e2e5e9cacfe9724b6dd33572a2e464.
Source SHA-256:

- rankine_modes.py: 7129582ada9fd7cbbab2e68316771977f5a8e4643d09c800d33405d63ae254ce
- verify_cst002.py: 75332e7efa954eb2f0debfd1f57addd7d2ffbb5256d6b91be712d679d1b53a2b
- verify_cst003.py: a97675497cc0244d6db80520d4857a78f7e8ddf83a903dc2cedf3c4d5aa8df45

Impact: additive conditional Rankine API, its direct tests and the P251
proposal. No existing canonical API, accepted claim, release, generated
documentation or shared numerical/governance machinery changes. GitNexus's
configured index is the separate main worktree; its report concerns P250 and
cannot describe this diff. Direct source search bounds consumers as recorded
in README.md. The staged-change selector in scripts/validate_changed.py
selects scoped mode and tests/test_rankine_modes.py. Its default CLI compares
commits and cannot see the uncommitted index; the staged selection uses its
own parse_name_status and choose_validation_scope functions.

Commands and receipts:

- 0029/verify_correspondence.py: exit 0, 10 exact checks, stdout.txt.
- repaired verify_cst002.py: exit 0, 30 checks, 0031/cst002.stdout.txt.
- repaired verify_cst003.py: exit 0, 23 checks, cst003.stdout.txt.
- 0030/kida_angle_construction.py: exit 0, 25 checks, 0030/stdout-corrected.txt.
- 0031/verify_rankine.py: exit 0, 28 checks, 0031/stdout.txt. Its equations
  also appear in CST002; these are related checks, not independent evidence.
- 0032/hexagon_action.py: exit 0, 24 checks, 0032/stdout.txt; the Cartesian
  Jacobian is independent of the reduced Hamiltonian normalization.
- 0032/collective_field_map.py: exit 0, 18 checks, 0032/map-stdout.txt.
- 0033/verify_mutual_kernel.py: exit 0, 6 checks, 0033/stdout.txt.
- 0034/verify_longwave.py: exit 0, 11 checks, 0034/stdout-density.txt;
  shared map with 0032 is corroborating algebra, not another microscopic input.
- scripts/validate.sh --pytest-scope tests/test_rankine_modes.py: exit 0,
  7 tests plus all fixed repository checks, validation.stdout.txt. All 1031
  memory records validate; 43 warnings are reported, not failures.
- Ruff on every new script/module/test and revised CST002/CST003: exit 0.
- git diff --check and git diff --cached --check: exit 0.

The workflow used the existing shared virtual environment through PYTHON;
no dependency change was made. Syntax-only lambda-to-def and dummy-variable
lint repairs preserve the captured predicates. The independent reviewer
requested one exterior-field documentation correction, now verified in
independent-review.md; equations and roots are unchanged. Existing N1 and
N4–N7 conditional receipts remain applicable and were not rerun.

No claim is promoted and the original exact smooth-Euler continuum objective
is not declared complete. This receipt preserves executable campaign progress.
