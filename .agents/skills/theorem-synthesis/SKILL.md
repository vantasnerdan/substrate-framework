---
name: theorem-synthesis
description: Compose accepted Substrate Framework claims into a new higher theorem and carry it through proof, review, registry promotion, and release. Use for synthesis campaigns, composite or synthesized claims, cross-sector theorem discovery, SymPy glue chains, Lean end-to-end formalization, conditional interpretive theorems, or planning a theorem promotion from existing accepted atoms.
---

# Theorem Synthesis

Promote a new positive statement by proving the missing glue between accepted
claims. Keep proof, empirical applicability, and scientific interpretation
distinct so each can advance without becoming ceremony for the others. Protect
the strongest useful theorem supported by the composition; do not make a proof
easier by shrinking the result to a tautology or isolated special value.

## Establish the target

Read `AGENTS.md`, `governance/releases/current.yaml`, the relevant entries in
`governance/claims.yaml`, and `.agents/skills/physics-erdos-loop/SKILL.md`.
Instantiate `memory-templates/theorem-synthesis.md` and a synthesis proposal
manifest before substantive proof work.

State:

- one exact higher theorem;
- at least two distinct accepted dependency claims;
- the structural gap not closed by any accepted claim;
- assumptions, exclusions, and intended consumers;
- `core` or `interpretive` layer, with a named hypothesis for the latter.

Do not re-review accepted dependencies. Audit only that the theorem uses their
exact accepted scope. There is no statement-length cap once the dependency and
assumption boundary is explicit.

If the composition needs a real hypothesis, state it and use the interpretive
layer rather than rejecting an otherwise useful conditional theorem. A missing
bridge is “not yet proved,” not a refutation of either accepted atom.

If the scientific mechanism is genuinely open, compare plausible mechanisms
under preregistered structural criteria. If the theorem statement is fixed,
declare `target_kind: fixed_theorem` and proceed with one sound proof route;
do not invent rival mechanisms to satisfy a form.

Run `scripts/find_synthesis_candidates.py` when cross-sector graph structure can
help discovery. Treat its ranking as a hint, never as evidence or a gate.

## Prove the glue

Choose the strongest practical proof backend for the exact statement:

- Use SymPy for exact identities, substitutions, eliminations, and finite
  symbolic chains.
- Use Lean for finite formal statements and prioritize an end-to-end theorem
  over disconnected attestations. Import shared framework definitions where
  practical, reject proof escapes, inspect the theorem and axiom footprint, and
  audit the map back to the physics statement.

Lean is complementary evidence, not a universal prerequisite. A kernel-checked
theorem does not need a performative input mutation; inspect its statement,
imports, axioms, and physical encoding. Apply mutation sensitivity to custom
translators, generated expressions, symbolic verifiers, and numerical evidence
where it can reveal a false green result.

Measurement, numerical evidence, or simulation may test whether assumptions
describe nature and whether consequences are accurate. They do not prove an
exact implication. Record multiple evidence modalities in
`verification_evidence` without collapsing their scopes.

## Promote the theorem

For the proposed claim, set `category: synthesized` and provide:

```yaml
composition:
  dependencies: [C-SECTOR-001, C-OTHER-001]
  structural_gap: <the missing accepted implication or identity>
  glue:
    method: lean  # or sympy
    artifact: <repository-relative proof path>
    entrypoint: <theorem name or verifier entrypoint>
```

Declare `exclusions`. For `layer: interpretive`, also declare
`hypothesis: {label: H-name, statement: ...}`. Interpretive claims may depend on
core claims; core claims may not depend on the interpretive layer.

Review the new theorem individually while preserving the prior acceptance of
its dependencies. Replay direct consumers and the changed graph boundary. Use
`scripts/validate_changed.py` for proportional validation: an additive leaf
synthesized theorem can remain scoped when no existing contract or consumer is
changed; foundation changes, shared machinery, uncertain impact, or altered
existing claims require the wider replay.

Freeze that boundary before review. An adjacent issue in an accepted atom or
historical artifact becomes a follow-up unless it directly invalidates the new
composition. Recheck only corrected glue and affected consumers; never turn a
synthesis review into recursive re-acceptance of its atoms or validation of the
review process.

Use one substantive review pass and one correction check. Report the strongest
supported theorem first, then any unsupported extension and minimum repair.
Do not commission a second review because counts, review prose, or evidence-role
labels changed.

A synthesis campaign succeeds only when its higher claim is accepted, entered
in the registry, and pinned in a release. Keep supporting implementation,
failed proof routes, and useful subclaims as commits on the synthesis campaign
branch. They do not open a partial PR gate. Submit one campaign PR only after
the theorem is positively complete or the proof-route space satisfies the
scientific-exhaustion certificate in `AGENTS.md`.
