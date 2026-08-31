from __future__ import annotations

import pytest

from substrate_framework.governance import (
    GovernanceError,
    validate_proposal,
    validate_registry,
    validate_release,
)
from substrate_framework.verification import (
    CheckFailure,
    CheckLedger,
    SuccessfulCheckTally,
)


def claim(claim_id: str, dependencies: list[str], accepted: bool = True) -> dict:
    return {
        "id": claim_id,
        "statement": f"statement for {claim_id}",
        "provenance": "proposal-P000",
        "verification": "symbolic_verified" if accepted else "unverified",
        "review": "accepted" if accepted else "unaudited",
        "compatibility": "native" if accepted else "unassessed",
        "epistemic": "active" if accepted else "proposed",
        "dependencies": dependencies,
        "evidence": ["tests/test_governance.py"] if accepted else [],
        "assumptions": [],
        "comparators": [],
        "accepted_in": "v-test" if accepted else None,
    }


def synthesized_claim(
    claim_id: str,
    dependencies: list[str],
    *,
    layer: str = "core",
    accepted: bool = True,
) -> dict:
    result = claim(claim_id, dependencies, accepted=accepted)
    result.update(
        {
            "category": "synthesized",
            "layer": layer,
            "exclusions": ["No substrate mechanism follows without an explicit hypothesis"],
            "composition": {
                "dependencies": dependencies,
                "structural_gap": "the accepted atoms have not been composed end to end",
                "glue": {
                    "method": "lean",
                    "artifact": "formal/SubstrateFramework/Glue.lean",
                    "entrypoint": "SubstrateFramework.compose_implications",
                },
            },
            "verification_evidence": [
                {
                    "method": "lean",
                    "artifact": "formal/SubstrateFramework/Glue.lean",
                    "scope": "the implication glue only",
                },
                {
                    "method": "measurement",
                    "artifact": "tests/test_governance.py",
                    "scope": "physical applicability, not the formal implication",
                },
            ],
        }
    )
    return result


def test_empty_registry_is_valid_bootstrap() -> None:
    assert validate_registry({"schema_version": 1, "claims": []}) == []


def test_accepted_dependency_closure_is_valid() -> None:
    data = {"schema_version": 1, "claims": [claim("C1", []), claim("C2", ["C1"])]}
    assert validate_registry(data) == ["C1", "C2"]


def test_synthesized_claim_composes_two_accepted_dependencies() -> None:
    data = {
        "schema_version": 1,
        "claims": [
            claim("C1", []),
            claim("C2", []),
            synthesized_claim("C3", ["C1", "C2"]),
        ],
    }

    assert validate_registry(data) == ["C1", "C2", "C3"]


def test_synthesized_claim_rejects_duplicate_or_unaccepted_atoms() -> None:
    duplicate = synthesized_claim("C3", ["C1", "C1"])
    data = {"schema_version": 1, "claims": [claim("C1", []), duplicate]}
    with pytest.raises(GovernanceError, match="two distinct dependencies"):
        validate_registry(data)

    proposed = claim("C2", [], accepted=False)
    composite = synthesized_claim("C3", ["C1", "C2"], accepted=False)
    data = {
        "schema_version": 1,
        "claims": [claim("C1", []), proposed, composite],
    }
    with pytest.raises(GovernanceError, match="depends on unaccepted C2"):
        validate_registry(data)


def test_interpretive_layer_is_explicit_and_cannot_feed_core() -> None:
    interpretation = synthesized_claim("C3", ["C1", "C2"], layer="interpretive")
    interpretation["hypothesis"] = {
        "label": "H-substrate",
        "statement": "the accepted effective fields arise from the proposed substrate",
    }
    core_consumer = claim("C4", ["C3"])
    data = {
        "schema_version": 1,
        "claims": [claim("C1", []), claim("C2", []), interpretation, core_consumer],
    }

    with pytest.raises(GovernanceError, match="core claim depends on interpretive"):
        validate_registry(data)

    data["claims"].pop()
    assert validate_registry(data) == ["C1", "C2", "C3"]


def test_accepted_claim_cannot_depend_on_proposal() -> None:
    data = {
        "schema_version": 1,
        "claims": [claim("C1", [], accepted=False), claim("C2", ["C1"])],
    }
    with pytest.raises(GovernanceError, match="depends on noncurrent"):
        validate_registry(data)


def test_dependency_cycle_fails() -> None:
    data = {"schema_version": 1, "claims": [claim("C1", ["C2"]), claim("C2", ["C1"])]}
    with pytest.raises(GovernanceError, match="dependency cycle"):
        validate_registry(data)


def test_mutation_gate_rejects_insensitive_check() -> None:
    ledger = CheckLedger("C-test")
    with pytest.raises(CheckFailure, match="insensitive"):
        ledger.mutation_sensitive("value", lambda _: True, 1, [2])


def test_successful_tally_formats_as_count_but_exits_with_status_zero() -> None:
    ledger = CheckLedger("C-test")
    ledger.check("one load-bearing assertion", True)
    tally = ledger.finish()
    assert isinstance(tally, SuccessfulCheckTally)
    assert tally.passed_count == 1
    assert f"{tally}" == "1"
    assert int(tally) == 0
    assert SystemExit(tally).code == 0


def test_supersession_preserves_historical_claim() -> None:
    old = claim("C1", [])
    old["epistemic"] = "superseded"
    old["compatibility"] = "conflict"
    new = claim("C2", [])
    new["supersedes"] = ["C1"]
    assert validate_registry({"schema_version": 1, "claims": [old, new]}) == ["C1", "C2"]


def test_proposal_requires_immutable_source_baseline() -> None:
    proposal = {
        "id": "P000",
        "base_release": None,
        "source_baseline": "substrate@6d1f4e0",
        "question": "derive a positive root claim",
        "invariants": ["normalized sine-Gordon convention"],
        "allowed_imports": ["real analysis"],
        "candidates": [
            {"id": "A", "description": "closed-form construction"},
            {"id": "B", "description": "independent transform construction"},
        ],
        "selection_criteria": ["exact dependency closure"],
        "claims_proposed": ["C-SG-001"],
        "comparators_blinded_until": "structural review complete",
        "status": "draft",
    }

    validate_proposal(proposal)
    proposal["source_baseline"] = ""
    with pytest.raises(GovernanceError, match="immutable source revision"):
        validate_proposal(proposal)


def campaign_proposal_v2() -> dict:
    return {
        "schema_version": 2,
        "id": "P900",
        "base_release": "v-test",
        "source_baseline": "substrate-framework@abc123",
        "question": "construct the requested positive object",
        "invariants": ["preserve the accepted normalization"],
        "allowed_imports": ["real analysis"],
        "candidates": [
            {"id": "A", "description": "direct construction"},
            {"id": "B", "description": "dual representation"},
        ],
        "selection_criteria": ["dependency closure"],
        "claims_proposed": ["C-P900-001"],
        "comparators_blinded_until": "structural freeze",
        "status": "active",
        "candidate_universe": {
            "scope": "all constructions allowed by the objective and imports",
            "frozen_from": ["user objective", "accepted invariants"],
            "route_families": ["direct", "dual"],
            "append_only_expansions": [],
        },
        "obligation_graph": {
            "nodes": [
                {
                    "id": "O1",
                    "positive_intent": "license the mathematical object",
                    "requires": [],
                    "pass_licenses": ["L-object"],
                    "does_not_license": ["L-numeric"],
                    "maximum_verdict": "OBJECT_EXISTS",
                    "failure_scope": "this candidate representation only",
                    "unlocks": ["O2"],
                    "status": "active",
                    "license_chain": {
                        "object": "typed candidate object",
                        "symmetry_or_conservation": "normalization identity",
                        "ensemble": "statics",
                        "variational_functional": "complete functional",
                        "admissible_space": "declared constrained space",
                        "representation_coverage": "full in-scope representation",
                        "observable": "existence predicate",
                        "numerical_representation": "not yet licensed",
                        "permitted_verdict": "OBJECT_EXISTS",
                    },
                },
                {
                    "id": "O2",
                    "positive_intent": "validate the observable",
                    "requires": ["L-object"],
                    "pass_licenses": ["L-numeric"],
                    "does_not_license": [],
                    "maximum_verdict": "CLAIM_ESTABLISHED",
                    "failure_scope": "the licensed numerical route",
                    "unlocks": [],
                    "status": "exhausted",
                    "license_chain": {
                        "object": "typed candidate object",
                        "symmetry_or_conservation": "normalization identity",
                        "ensemble": "statics",
                        "variational_functional": "complete functional",
                        "admissible_space": "declared constrained space",
                        "representation_coverage": "full in-scope representation",
                        "observable": "scale-relative residual",
                        "numerical_representation": "crossed refinement design",
                        "permitted_verdict": "CLAIM_ESTABLISHED",
                    },
                },
            ]
        },
        "license_registry": [
            {
                "id": "L-object",
                "proposition": "the typed object is admissible",
                "status": "unearned",
            },
            {
                "id": "L-numeric",
                "proposition": "the discretization resolves the observable",
                "status": "unearned",
            },
        ],
        "route_frontier": {
            "active_obligation": "O1",
            "considered": ["direct", "dual"],
            "tried": [],
            "failure_generated": [],
            "remaining": ["direct", "dual"],
        },
        "execution_state": "active",
        "objective_state": "active",
        "exhaustion_certificate": {},
    }


def test_campaign_proposal_v2_accepts_active_continuation_state() -> None:
    validate_proposal(campaign_proposal_v2())


def test_campaign_proposal_v2_requires_machine_terminal_fields() -> None:
    proposal = campaign_proposal_v2()
    proposal.pop("route_frontier")
    with pytest.raises(GovernanceError, match="schema v2 missing fields"):
        validate_proposal(proposal)


def test_campaign_proposal_v2_rejects_unknown_license_edges() -> None:
    proposal = campaign_proposal_v2()
    proposal["obligation_graph"]["nodes"][0]["pass_licenses"] = ["L-unknown"]
    with pytest.raises(GovernanceError, match="unknown license ids"):
        validate_proposal(proposal)


def test_campaign_terminal_success_requires_every_obligation_established() -> None:
    proposal = campaign_proposal_v2()
    proposal["execution_state"] = "terminal_success"
    proposal["objective_state"] = "complete"
    proposal["route_frontier"]["active_obligation"] = None
    proposal["route_frontier"]["remaining"] = []
    proposal["obligation_graph"]["nodes"][0]["status"] = "established"
    proposal["license_registry"][0].update(
        {"status": "earned", "evidence": "proof.md", "earned_by": "O1"}
    )
    with pytest.raises(GovernanceError, match="every obligation established"):
        validate_proposal(proposal)


def test_campaign_terminal_exhaustion_requires_no_routes_and_independent_adversary() -> None:
    proposal = campaign_proposal_v2()
    proposal["execution_state"] = "terminal_exhaustion"
    proposal["route_frontier"]["active_obligation"] = None
    proposal["route_frontier"]["remaining"] = []
    proposal["route_frontier"]["considered"] = [
        "historical direct route",
        "external dual route",
        "failure-derived transform",
    ]
    proposal["route_frontier"]["tried"] = list(
        proposal["route_frontier"]["considered"]
    )
    proposal["route_frontier"]["failure_generated"] = [
        "failure-derived transform"
    ]
    proposal["obligation_graph"]["nodes"][0]["status"] = "exhausted"
    proposal["exhaustion_certificate"] = {
        "historical_routes": ["historical direct route"],
        "external_routes": ["external dual route"],
        "failure_generated_routes": ["failure-derived transform"],
        "equivalence_partition": [{"class": "direct", "members": ["A"]}],
        "adversarial_generation_artifact": "attempts/adversarial.md",
        "adversarial_reviewer": "independent-agent",
        "adversarial_reviewer_role": "author",
        "infinite_class_coverage": ["analytic no-go for continuous variants"],
        "routes_remaining": [],
        "route_verdicts": [
            {
                "route": route,
                "verdict": "route refuted",
                "evidence": "attempts/routes.md",
                "continuation": "method, representation, and concept ladder complete",
            }
            for route in proposal["route_frontier"]["considered"]
        ],
        "review": "coverage checked against the frozen universe",
    }
    with pytest.raises(GovernanceError, match="non-author non-implementer"):
        validate_proposal(proposal)

    proposal["exhaustion_certificate"]["adversarial_reviewer_role"] = (
        "non_author_non_implementer"
    )
    validate_proposal(proposal)

    proposal["route_frontier"]["remaining"] = ["untried route"]
    with pytest.raises(GovernanceError, match="cannot retain routes_remaining"):
        validate_proposal(proposal)


def test_fixed_theorem_synthesis_needs_one_route_not_competing_mechanisms() -> None:
    proposal = {
        "id": "P000",
        "base_release": "v-test",
        "source_baseline": "substrate-framework@abc123",
        "question": "compose two accepted implications into one higher theorem",
        "invariants": ["accepted dependencies are not reopened"],
        "allowed_imports": ["C1", "C2"],
        "candidates": [{"id": "proof", "description": "Lean implication chain"}],
        "selection_criteria": ["exact statement match"],
        "claims_proposed": ["C3"],
        "comparators_blinded_until": "not applicable to an exact theorem",
        "status": "draft",
        "campaign_type": "synthesis",
        "target_kind": "fixed_theorem",
        "structural_gap": "C1 and C2 have no accepted composition",
        "composition_dependencies": ["C1", "C2"],
    }

    validate_proposal(proposal)
    proposal["composition_dependencies"] = ["C1"]
    with pytest.raises(GovernanceError, match="two distinct accepted dependencies"):
        validate_proposal(proposal)


def test_release_requires_dependency_closed_claim_set() -> None:
    registry = {
        "schema_version": 1,
        "claims": [claim("C1", []), claim("C2", ["C1"])],
    }
    release = {
        "schema_version": 1,
        "release": "v-test",
        "source_baseline": "source@abc123",
        "released_at": "2026-08-01T00:00:00Z",
        "accepted_claims": ["C2"],
    }

    with pytest.raises(GovernanceError, match="dependencies outside release"):
        validate_release(release, registry)

    release["accepted_claims"] = ["C1", "C2"]
    assert validate_release(release, registry, require_current_set=True) == ["C1", "C2"]
