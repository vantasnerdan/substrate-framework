# Tool receipts - 0126 (beacon)

- Registry survey (read-only): schema keys, review/verification/
  epistemic/compatibility value sets, C-EUL-* zero-hit check —
  commands in transcript, no writes.
- Conformance probe (inline, exit 0): drafts.yaml parses; both IDs
  absent from registry; zero extra keys vs accepted entries; all
  evidence paths exist; C-EUL-002 dep resolves in-draft. Key shape
  matches C-VTX-002 (numeric_evidence) exactly modulo accepted_in
  (correctly omitted for drafts).
- No registry, release, proposal, or P251 file touched (git status
  clean outside attempt dirs + herd inbox on landing).
