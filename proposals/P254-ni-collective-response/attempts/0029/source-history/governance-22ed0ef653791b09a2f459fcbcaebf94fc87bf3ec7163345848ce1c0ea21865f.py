"""Validation and rendering for the scientific claim graph.

The functions here are deliberately importable. Campaign scripts should consume
canonical definitions from the package instead of copying constants or helpers.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

import yaml


class GovernanceError(ValueError):
    """Raised when a registry or proposal violates framework governance."""


REQUIRED_CLAIM_FIELDS = {
    "id",
    "statement",
    "provenance",
    "verification",
    "review",
    "compatibility",
    "epistemic",
    "dependencies",
    "evidence",
    "assumptions",
    "comparators",
    "accepted_in",
}

REQUIRED_PROPOSAL_FIELDS = {
    "id",
    "base_release",
    "source_baseline",
    "question",
    "invariants",
    "allowed_imports",
    "candidates",
    "selection_criteria",
    "claims_proposed",
    "comparators_blinded_until",
    "status",
}

REQUIRED_PROPOSAL_V2_FIELDS = {
    "candidate_universe",
    "obligation_graph",
    "license_registry",
    "route_frontier",
    "execution_state",
    "objective_state",
    "exhaustion_certificate",
}

REQUIRED_RELEASE_FIELDS = {
    "schema_version",
    "release",
    "source_baseline",
    "released_at",
    "accepted_claims",
}


def load_yaml(path: str | Path) -> dict[str, Any]:
    """Load a YAML mapping from *path*."""

    source = Path(path)
    data = yaml.safe_load(source.read_text(encoding="utf-8"))
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise GovernanceError(f"{source}: top level must be a mapping")
    return data


def _as_string_list(value: Any, field: str, owner: str) -> list[str]:
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise GovernanceError(f"{owner}: {field} must be a list of strings")
    return value


def _nonempty_string(value: Any, field: str, owner: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise GovernanceError(f"{owner}: {field} must be a non-empty string")
    return value


def _validate_verification_evidence(claim: dict[str, Any], claim_id: str) -> None:
    records = claim.get("verification_evidence", [])
    if not isinstance(records, list):
        raise GovernanceError(
            f"{claim_id}: verification_evidence must be a list of mappings"
        )
    allowed_methods = {
        "analytic",
        "sympy",
        "lean",
        "numeric",
        "simulation",
        "measurement",
    }
    for index, record in enumerate(records):
        owner = f"{claim_id}.verification_evidence[{index}]"
        if not isinstance(record, dict):
            raise GovernanceError(f"{owner} must be a mapping")
        missing = {"method", "artifact", "scope"} - record.keys()
        if missing:
            raise GovernanceError(f"{owner} missing fields: {sorted(missing)}")
        if record["method"] not in allowed_methods:
            raise GovernanceError(
                f"{owner}: invalid verification method {record['method']!r}"
            )
        _nonempty_string(record["artifact"], "artifact", owner)
        _nonempty_string(record["scope"], "scope", owner)


def _detect_dependency_cycles(graph: dict[str, list[str]]) -> None:
    active: set[str] = set()
    complete: set[str] = set()

    def visit(node: str, trail: tuple[str, ...]) -> None:
        if node in active:
            cycle = " -> ".join((*trail, node))
            raise GovernanceError(f"claim dependency cycle: {cycle}")
        if node in complete:
            return
        active.add(node)
        for dependency in graph[node]:
            visit(dependency, (*trail, node))
        active.remove(node)
        complete.add(node)

    for claim_id in graph:
        visit(claim_id, ())


def _validate_proposal_v2(data: dict[str, Any], source: str) -> None:
    """Validate the machine-readable campaign continuation and terminal gates."""

    missing = REQUIRED_PROPOSAL_V2_FIELDS - data.keys()
    if missing:
        raise GovernanceError(f"{source} schema v2 missing fields: {sorted(missing)}")

    universe = data["candidate_universe"]
    if not isinstance(universe, dict):
        raise GovernanceError(f"{source}: candidate_universe must be a mapping")
    _nonempty_string(universe.get("scope"), "candidate_universe.scope", source)
    frozen_from = _as_string_list(
        universe.get("frozen_from"), "candidate_universe.frozen_from", source
    )
    route_families = _as_string_list(
        universe.get("route_families"), "candidate_universe.route_families", source
    )
    _as_string_list(
        universe.get("append_only_expansions"),
        "candidate_universe.append_only_expansions",
        source,
    )
    if not frozen_from or not route_families:
        raise GovernanceError(
            f"{source}: candidate_universe needs frozen sources and route families"
        )

    licenses = data["license_registry"]
    if not isinstance(licenses, list) or not licenses:
        raise GovernanceError(f"{source}: license_registry must be a non-empty list")
    license_ids: set[str] = set()
    for index, license_record in enumerate(licenses):
        owner = f"{source}.license_registry[{index}]"
        if not isinstance(license_record, dict):
            raise GovernanceError(f"{owner} must be a mapping")
        license_id = _nonempty_string(license_record.get("id"), "id", owner)
        if license_id in license_ids:
            raise GovernanceError(f"{source}: duplicate license id {license_id!r}")
        license_ids.add(license_id)
        _nonempty_string(license_record.get("proposition"), "proposition", owner)
        status = license_record.get("status")
        if status not in {"unearned", "earned", "refuted", "blocked"}:
            raise GovernanceError(f"{owner}: invalid status {status!r}")
        if status != "unearned":
            _nonempty_string(license_record.get("evidence"), "evidence", owner)
            _nonempty_string(license_record.get("earned_by"), "earned_by", owner)

    graph = data["obligation_graph"]
    if not isinstance(graph, dict) or not isinstance(graph.get("nodes"), list):
        raise GovernanceError(f"{source}: obligation_graph.nodes must be a list")
    nodes = graph["nodes"]
    if not nodes:
        raise GovernanceError(f"{source}: obligation_graph must contain a node")
    node_ids: set[str] = set()
    node_records: list[dict[str, Any]] = []
    chain_fields = (
        "object",
        "symmetry_or_conservation",
        "ensemble",
        "variational_functional",
        "admissible_space",
        "representation_coverage",
        "observable",
        "numerical_representation",
        "permitted_verdict",
    )
    for index, node in enumerate(nodes):
        owner = f"{source}.obligation_graph.nodes[{index}]"
        if not isinstance(node, dict):
            raise GovernanceError(f"{owner} must be a mapping")
        node_id = _nonempty_string(node.get("id"), "id", owner)
        if node_id in node_ids:
            raise GovernanceError(f"{source}: duplicate obligation id {node_id!r}")
        node_ids.add(node_id)
        node_records.append(node)
        for field in ("positive_intent", "maximum_verdict", "failure_scope"):
            _nonempty_string(node.get(field), field, owner)
        for field in ("requires", "pass_licenses", "does_not_license", "unlocks"):
            _as_string_list(node.get(field), field, owner)
        if node.get("status") not in {
            "pending",
            "active",
            "established",
            "exhausted",
            "refuted",
        }:
            raise GovernanceError(f"{owner}: invalid status {node.get('status')!r}")
        chain = node.get("license_chain")
        if not isinstance(chain, dict):
            raise GovernanceError(f"{owner}: license_chain must be a mapping")
        for field in chain_fields:
            _nonempty_string(chain.get(field), f"license_chain.{field}", owner)

    for node in node_records:
        owner = f"{source}.obligation_graph.{node['id']}"
        unknown_licenses = (
            set(node["requires"]) | set(node["pass_licenses"])
        ) - license_ids
        if unknown_licenses:
            raise GovernanceError(
                f"{owner}: unknown license ids {sorted(unknown_licenses)}"
            )
        unknown_nodes = set(node["unlocks"]) - node_ids
        if unknown_nodes:
            raise GovernanceError(
                f"{owner}: unknown unlocked obligations {sorted(unknown_nodes)}"
            )
        if node["status"] == "established":
            unearned_outputs = [
                license_id
                for license_id in node["pass_licenses"]
                if next(
                    item["status"]
                    for item in licenses
                    if item["id"] == license_id
                )
                != "earned"
            ]
            if unearned_outputs:
                raise GovernanceError(
                    f"{owner}: established obligation has unearned pass licenses "
                    f"{sorted(unearned_outputs)}"
                )

    frontier = data["route_frontier"]
    if not isinstance(frontier, dict):
        raise GovernanceError(f"{source}: route_frontier must be a mapping")
    for field in ("considered", "tried", "failure_generated", "remaining"):
        _as_string_list(frontier.get(field), f"route_frontier.{field}", source)
    active_obligation = frontier.get("active_obligation")
    if active_obligation is not None and active_obligation not in node_ids:
        raise GovernanceError(
            f"{source}: route_frontier.active_obligation is not an obligation id"
        )

    execution_state = data["execution_state"]
    objective_state = data["objective_state"]
    if execution_state not in {"active", "terminal_success", "terminal_exhaustion"}:
        raise GovernanceError(f"{source}: invalid execution_state {execution_state!r}")
    if objective_state not in {"active", "complete"}:
        raise GovernanceError(f"{source}: invalid objective_state {objective_state!r}")

    active_nodes = [node["id"] for node in node_records if node["status"] == "active"]
    if execution_state == "active":
        if objective_state != "active" or not active_nodes:
            raise GovernanceError(
                f"{source}: active execution needs an active objective and obligation"
            )
        if active_obligation not in active_nodes:
            raise GovernanceError(
                f"{source}: active_obligation must identify an active obligation"
            )
        return

    pending_nodes = [node["id"] for node in node_records if node["status"] == "pending"]
    if active_nodes or pending_nodes or active_obligation is not None:
        raise GovernanceError(
            f"{source}: terminal campaign cannot retain an active or pending obligation"
        )
    if frontier["remaining"]:
        raise GovernanceError(
            f"{source}: terminal campaign cannot retain routes_remaining"
        )
    if execution_state == "terminal_success":
        if objective_state != "complete" or any(
            node["status"] != "established" for node in node_records
        ):
            raise GovernanceError(
                f"{source}: terminal success requires the objective and every obligation established"
            )
        return

    if objective_state != "active":
        raise GovernanceError(
            f"{source}: terminal exhaustion leaves the positive objective active"
        )
    certificate = data["exhaustion_certificate"]
    if not isinstance(certificate, dict):
        raise GovernanceError(f"{source}: exhaustion_certificate must be a mapping")
    for field in (
        "historical_routes",
        "external_routes",
        "failure_generated_routes",
        "infinite_class_coverage",
        "routes_remaining",
    ):
        _as_string_list(certificate.get(field), f"exhaustion_certificate.{field}", source)
    route_verdicts = certificate.get("route_verdicts")
    if not isinstance(route_verdicts, list) or not route_verdicts:
        raise GovernanceError(
            f"{source}: exhaustion_certificate.route_verdicts must be non-empty"
        )
    verdict_routes: set[str] = set()
    for index, route_verdict in enumerate(route_verdicts):
        owner = f"{source}.exhaustion_certificate.route_verdicts[{index}]"
        if not isinstance(route_verdict, dict):
            raise GovernanceError(f"{owner} must be a mapping")
        route = _nonempty_string(route_verdict.get("route"), "route", owner)
        if route in verdict_routes:
            raise GovernanceError(f"{source}: duplicate route verdict {route!r}")
        verdict_routes.add(route)
        _nonempty_string(route_verdict.get("verdict"), "verdict", owner)
        _nonempty_string(route_verdict.get("evidence"), "evidence", owner)
        _nonempty_string(route_verdict.get("continuation"), "continuation", owner)
    considered_routes = set(frontier["considered"])
    if verdict_routes != considered_routes:
        raise GovernanceError(
            f"{source}: exhaustion route verdicts must cover every considered route exactly"
        )
    inventoried_routes = set(certificate["historical_routes"]) | set(
        certificate["external_routes"]
    ) | set(certificate["failure_generated_routes"])
    if not inventoried_routes <= considered_routes:
        raise GovernanceError(
            f"{source}: exhaustion inventories contain unconsidered routes"
        )
    partition = certificate.get("equivalence_partition")
    if not isinstance(partition, list) or not partition:
        raise GovernanceError(
            f"{source}: exhaustion_certificate.equivalence_partition must be non-empty"
        )
    if certificate["routes_remaining"]:
        raise GovernanceError(
            f"{source}: exhaustion certificate cannot retain routes_remaining"
        )
    for field in (
        "adversarial_generation_artifact",
        "adversarial_reviewer",
        "review",
    ):
        _nonempty_string(certificate.get(field), f"exhaustion_certificate.{field}", source)
    if certificate.get("adversarial_reviewer_role") != "non_author_non_implementer":
        raise GovernanceError(
            f"{source}: exhaustion adversary must be a non-author non-implementer"
        )


def validate_registry(data: dict[str, Any]) -> list[str]:
    """Validate claim identity, status axes, dependency closure, and promotion gates."""

    if data.get("schema_version") != 1:
        raise GovernanceError("claims registry must declare schema_version: 1")
    claims = data.get("claims")
    if not isinstance(claims, list):
        raise GovernanceError("claims must be a list")

    by_id: dict[str, dict[str, Any]] = {}
    for index, claim in enumerate(claims):
        owner = f"claims[{index}]"
        if not isinstance(claim, dict):
            raise GovernanceError(f"{owner} must be a mapping")
        missing = REQUIRED_CLAIM_FIELDS - claim.keys()
        if missing:
            raise GovernanceError(f"{owner} missing fields: {sorted(missing)}")
        claim_id = claim["id"]
        if not isinstance(claim_id, str) or not claim_id:
            raise GovernanceError(f"{owner}.id must be a non-empty string")
        if claim_id in by_id:
            raise GovernanceError(f"duplicate claim id: {claim_id}")
        by_id[claim_id] = claim

    graph: dict[str, list[str]] = {}
    allowed_review = {"unaudited", "audited", "accepted", "rejected"}
    allowed_compatibility = {"unassessed", "native", "compatible_extension", "conflict"}
    allowed_epistemic = {"proposed", "active", "qualified", "superseded", "refuted"}
    allowed_verification = {
        "unverified",
        "symbolic_verified",
        "formal_verified",
        "numeric_evidence",
        "simulation_evidence",
    }

    superseding_targets: set[str] = set()
    for claim_id, claim in by_id.items():
        dependencies = _as_string_list(claim["dependencies"], "dependencies", claim_id)
        _as_string_list(claim["evidence"], "evidence", claim_id)
        _as_string_list(claim["assumptions"], "assumptions", claim_id)
        _as_string_list(claim["comparators"], "comparators", claim_id)
        _validate_verification_evidence(claim, claim_id)
        graph[claim_id] = dependencies
        unknown = set(dependencies) - by_id.keys()
        if unknown:
            raise GovernanceError(f"{claim_id}: unknown dependencies {sorted(unknown)}")
        if claim["review"] not in allowed_review:
            raise GovernanceError(f"{claim_id}: invalid review status {claim['review']!r}")
        if claim["compatibility"] not in allowed_compatibility:
            raise GovernanceError(
                f"{claim_id}: invalid compatibility status {claim['compatibility']!r}"
            )
        if claim["epistemic"] not in allowed_epistemic:
            raise GovernanceError(f"{claim_id}: invalid epistemic status {claim['epistemic']!r}")
        if claim["verification"] not in allowed_verification:
            raise GovernanceError(
                f"{claim_id}: invalid verification status {claim['verification']!r}"
            )

        category = claim.get("category", "standard")
        layer = claim.get("layer", "core")
        if category not in {"standard", "synthesized"}:
            raise GovernanceError(f"{claim_id}: invalid claim category {category!r}")
        if layer not in {"core", "interpretive"}:
            raise GovernanceError(f"{claim_id}: invalid claim layer {layer!r}")

        if category == "synthesized" or layer == "interpretive":
            if "exclusions" not in claim:
                raise GovernanceError(
                    f"{claim_id}: synthesized and interpretive claims must declare exclusions"
                )
            _as_string_list(claim["exclusions"], "exclusions", claim_id)

        if category == "synthesized":
            composition = claim.get("composition")
            if not isinstance(composition, dict):
                raise GovernanceError(
                    f"{claim_id}: synthesized claim needs a composition mapping"
                )
            component_ids = _as_string_list(
                composition.get("dependencies"),
                "composition.dependencies",
                claim_id,
            )
            if len(component_ids) < 2 or len(component_ids) != len(set(component_ids)):
                raise GovernanceError(
                    f"{claim_id}: composition needs at least two distinct dependencies"
                )
            if set(component_ids) != set(dependencies) or len(component_ids) != len(
                dependencies
            ):
                raise GovernanceError(
                    f"{claim_id}: composition dependencies must match claim dependencies"
                )
            _nonempty_string(
                composition.get("structural_gap"),
                "composition.structural_gap",
                claim_id,
            )
            glue = composition.get("glue")
            if not isinstance(glue, dict):
                raise GovernanceError(
                    f"{claim_id}: composition.glue must be a mapping"
                )
            if glue.get("method") not in {"sympy", "lean"}:
                raise GovernanceError(
                    f"{claim_id}: composition.glue.method must be sympy or lean"
                )
            _nonempty_string(glue.get("artifact"), "composition.glue.artifact", claim_id)
            _nonempty_string(
                glue.get("entrypoint"), "composition.glue.entrypoint", claim_id
            )
            for dependency in component_ids:
                dependency_claim = by_id[dependency]
                if (
                    dependency_claim["accepted_in"] is None
                    or dependency_claim["review"] != "accepted"
                    or dependency_claim["epistemic"] not in {"active", "qualified"}
                ):
                    raise GovernanceError(
                        f"{claim_id}: synthesized claim depends on unaccepted {dependency}"
                    )

        if layer == "interpretive":
            hypothesis = claim.get("hypothesis")
            if not isinstance(hypothesis, dict):
                raise GovernanceError(
                    f"{claim_id}: interpretive claim needs a hypothesis mapping"
                )
            _nonempty_string(hypothesis.get("label"), "hypothesis.label", claim_id)
            _nonempty_string(
                hypothesis.get("statement"), "hypothesis.statement", claim_id
            )
            for dependency in dependencies:
                if by_id[dependency].get("layer", "core") != "core":
                    raise GovernanceError(
                        f"{claim_id}: interpretive claim must depend only on core claims"
                    )
        else:
            interpretive_dependencies = [
                dependency
                for dependency in dependencies
                if by_id[dependency].get("layer", "core") == "interpretive"
            ]
            if interpretive_dependencies:
                raise GovernanceError(
                    f"{claim_id}: core claim depends on interpretive claims "
                    f"{sorted(interpretive_dependencies)}"
                )

        accepted_in = claim["accepted_in"]
        was_accepted = accepted_in is not None
        is_current = was_accepted and claim["epistemic"] in {"active", "qualified"}
        if was_accepted:
            if not isinstance(accepted_in, str) or not accepted_in:
                raise GovernanceError(f"{claim_id}: accepted_in must be null or a release id")
            if claim["review"] != "accepted":
                raise GovernanceError(f"{claim_id}: accepted claim must have review: accepted")
            if claim["verification"] == "unverified":
                raise GovernanceError(f"{claim_id}: accepted claim must have verifier evidence")
            if not claim["evidence"]:
                raise GovernanceError(f"{claim_id}: accepted claim must cite evidence")
        if is_current:
            if claim["compatibility"] not in {"native", "compatible_extension"}:
                raise GovernanceError(f"{claim_id}: current accepted claim must fit the framework")
            for dependency in dependencies:
                dependency_claim = by_id[dependency]
                if dependency_claim["accepted_in"] is None or dependency_claim["epistemic"] not in {
                    "active",
                    "qualified",
                }:
                    raise GovernanceError(
                        f"{claim_id}: current claim depends on noncurrent {dependency}"
                    )

        challenges = _as_string_list(claim.get("challenges", []), "challenges", claim_id)
        supersedes = _as_string_list(claim.get("supersedes", []), "supersedes", claim_id)
        unknown_relationships = (set(challenges) | set(supersedes)) - by_id.keys()
        if unknown_relationships:
            raise GovernanceError(
                f"{claim_id}: unknown relationship targets {sorted(unknown_relationships)}"
            )
        if supersedes:
            if not is_current or claim["review"] != "accepted":
                raise GovernanceError(
                    f"{claim_id}: only current accepted claims may supersede; proposals use challenges"
                )
            superseding_targets.update(supersedes)

    _detect_dependency_cycles(graph)
    for target in superseding_targets:
        if by_id[target]["epistemic"] != "superseded":
            raise GovernanceError(
                f"{target}: target of supersedes must have epistemic: superseded"
            )
    return sorted(by_id)


def validate_proposal(data: dict[str, Any], source: str = "proposal") -> None:
    """Validate the pre-registered inputs required before a campaign begins."""

    missing = REQUIRED_PROPOSAL_FIELDS - data.keys()
    if missing:
        raise GovernanceError(f"{source} missing fields: {sorted(missing)}")
    if data["status"] not in {"draft", "active", "in_review", "accepted", "rejected", "rework"}:
        raise GovernanceError(f"{source}: invalid status {data['status']!r}")
    if not isinstance(data["source_baseline"], str) or not data["source_baseline"].strip():
        raise GovernanceError(f"{source}: source_baseline must name an immutable source revision")
    for field in ("invariants", "allowed_imports", "selection_criteria", "claims_proposed"):
        _as_string_list(data[field], field, source)
    campaign_type = data.get("campaign_type", "discovery")
    target_kind = data.get("target_kind", "mechanism_selection")
    if campaign_type not in {"discovery", "synthesis"}:
        raise GovernanceError(f"{source}: invalid campaign_type {campaign_type!r}")
    if target_kind not in {"mechanism_selection", "fixed_theorem"}:
        raise GovernanceError(f"{source}: invalid target_kind {target_kind!r}")
    if campaign_type == "synthesis":
        if target_kind != "fixed_theorem":
            raise GovernanceError(
                f"{source}: synthesis campaign must declare target_kind: fixed_theorem"
            )
        _nonempty_string(data.get("structural_gap"), "structural_gap", source)
        components = _as_string_list(
            data.get("composition_dependencies"),
            "composition_dependencies",
            source,
        )
        if len(components) < 2 or len(components) != len(set(components)):
            raise GovernanceError(
                f"{source}: synthesis campaign needs at least two distinct accepted dependencies"
            )
        if len(data["claims_proposed"]) != 1:
            raise GovernanceError(
                f"{source}: synthesis campaign must target exactly one higher claim"
            )
    candidates = data["candidates"]
    uniqueness = data.get("uniqueness_evidence")
    if not isinstance(candidates, list) or not candidates:
        raise GovernanceError(f"{source}: register at least one candidate approach")
    if (
        target_kind == "mechanism_selection"
        and len(candidates) < 2
        and not (isinstance(uniqueness, str) and uniqueness.strip())
    ):
        raise GovernanceError(
            f"{source}: register at least two candidates or cite uniqueness_evidence"
        )
    for index, candidate in enumerate(candidates):
        if not isinstance(candidate, dict) or not {"id", "description"} <= candidate.keys():
            raise GovernanceError(
                f"{source}: candidates[{index}] needs id and description"
            )
    schema_version = data.get("schema_version", 1)
    if schema_version not in {1, 2}:
        raise GovernanceError(f"{source}: unsupported schema_version {schema_version!r}")
    if schema_version == 2:
        _validate_proposal_v2(data, source)


def validate_release(
    data: dict[str, Any],
    registry: dict[str, Any],
    *,
    source: str = "release",
    require_current_set: bool = False,
) -> list[str]:
    """Validate a pinned release claim set and its dependency closure."""

    missing = REQUIRED_RELEASE_FIELDS - data.keys()
    if missing:
        raise GovernanceError(f"{source} missing fields: {sorted(missing)}")
    if data["schema_version"] != 1:
        raise GovernanceError(f"{source}: schema_version must be 1")

    release_id = data["release"]
    source_baseline = data["source_baseline"]
    released_at = data["released_at"]
    release_ids = _as_string_list(data["accepted_claims"], "accepted_claims", source)
    if len(release_ids) != len(set(release_ids)):
        raise GovernanceError(f"{source}: accepted_claims must not contain duplicates")

    validate_registry(registry)
    by_id = {claim["id"]: claim for claim in registry["claims"]}
    if release_id is None:
        if source_baseline is not None or released_at is not None or release_ids:
            raise GovernanceError(f"{source}: null release must have null source and no claims")
        return []
    if not isinstance(release_id, str) or not release_id:
        raise GovernanceError(f"{source}: release must be null or a non-empty id")
    if not isinstance(source_baseline, str) or not source_baseline.strip():
        raise GovernanceError(f"{source}: accepted release must pin source_baseline")
    if not isinstance(released_at, str) or not released_at.strip():
        raise GovernanceError(f"{source}: accepted release must record released_at")

    unknown = set(release_ids) - by_id.keys()
    if unknown:
        raise GovernanceError(f"{source}: unknown claims {sorted(unknown)}")
    release_set = set(release_ids)
    for claim_id in release_ids:
        claim = by_id[claim_id]
        if claim["accepted_in"] is None or claim["review"] != "accepted":
            raise GovernanceError(f"{source}: {claim_id} is not an accepted claim")
        missing_dependencies = set(claim["dependencies"]) - release_set
        if missing_dependencies:
            raise GovernanceError(
                f"{source}: {claim_id} has dependencies outside release: "
                f"{sorted(missing_dependencies)}"
            )

    if require_current_set:
        current_ids = {
            claim_id
            for claim_id, claim in by_id.items()
            if claim["accepted_in"] is not None
            and claim["epistemic"] in {"active", "qualified"}
        }
        if release_set != current_ids:
            missing_current = current_ids - release_set
            historical_only = release_set - current_ids
            raise GovernanceError(
                f"{source}: current set mismatch; missing {sorted(missing_current)}, "
                f"noncurrent {sorted(historical_only)}"
            )
    return release_ids


def render_claim_memory(claim: dict[str, Any], released_at: str) -> str:
    """Render a deterministic accepted-claim memory entry from registry state."""

    status = "active" if claim["epistemic"] in {"active", "qualified"} else "archived"
    frontmatter = {
        "description": f"Accepted framework claim {claim['id']}",
        "author": "framework-registry",
        "created": released_at,
        "updated": released_at,
        "tags": ["substrate-framework", "accepted-claim", claim["id"]],
        "category": "claims",
        "confidence": "established",
        "status": status,
    }
    header = yaml.safe_dump(frontmatter, sort_keys=False).strip()
    dependencies = ", ".join(claim["dependencies"]) or "none"
    assumptions = ", ".join(claim["assumptions"]) or "none"
    comparators = ", ".join(claim["comparators"]) or "none"
    evidence = "\n".join(f"- `{item}`" for item in claim["evidence"])
    theorem_metadata = ""
    if "category" in claim or "layer" in claim:
        theorem_metadata = (
            "\n## Theorem Classification\n"
            f"Category is `{claim.get('category', 'standard')}`; layer is "
            f"`{claim.get('layer', 'core')}`.\n"
        )
        composition = claim.get("composition")
        if isinstance(composition, dict):
            glue = composition["glue"]
            theorem_metadata += (
                f"Structural gap: {composition['structural_gap']}\n\n"
                f"Glue proof: `{glue['method']}` at `{glue['artifact']}` "
                f"entrypoint `{glue['entrypoint']}`.\n"
            )
        hypothesis = claim.get("hypothesis")
        if isinstance(hypothesis, dict):
            theorem_metadata += (
                f"\nConditional hypothesis `{hypothesis['label']}`: "
                f"{hypothesis['statement']}\n"
            )
    return (
        f"---\n{header}\n---\n"
        f"# {claim['id']}\n\n"
        "## Statement\n"
        "The accepted statement is reproduced exactly from the claim registry.\n\n"
        f"{claim['statement']}\n\n"
        "## Status Axes\n"
        "The four governance axes remain independent.\n\n"
        f"Verification is `{claim['verification']}`; review is `{claim['review']}`; "
        f"compatibility is `{claim['compatibility']}`; epistemic status is "
        f"`{claim['epistemic']}`.\n"
        f"{theorem_metadata}\n"
        "## Dependency and Import Closure\n"
        "The registry records the accepted closure and declared non-claim inputs.\n\n"
        f"Dependencies: {dependencies}. Assumptions: {assumptions}. "
        f"Comparators: {comparators}.\n\n"
        "## Provenance and Evidence\n"
        "The accepted release and immutable campaign evidence are the authoritative pointers.\n\n"
        f"Accepted in `{claim['accepted_in']}` with provenance `{claim['provenance']}`.\n\n"
        f"{evidence}\n"
    )


def render_release_memory(
    release: dict[str, Any], registry: dict[str, Any]
) -> str:
    """Render a deterministic accepted-release memory entry."""

    by_id = {claim["id"]: claim for claim in registry["claims"]}
    release_id = release["release"]
    released_at = release["released_at"]
    frontmatter = {
        "description": f"Accepted framework release {release_id}",
        "author": "framework-registry",
        "created": released_at,
        "updated": released_at,
        "tags": ["substrate-framework", "accepted-release", release_id],
        "category": "releases",
        "confidence": "established",
        "status": "active",
    }
    header = yaml.safe_dump(frontmatter, sort_keys=False).strip()
    claim_lines = "\n".join(
        f"- `{claim_id}` — {by_id[claim_id]['statement']}"
        for claim_id in release["accepted_claims"]
    )
    return (
        f"---\n{header}\n---\n"
        f"# Release {release_id}\n\n"
        "## Source Boundary\n"
        "This release pins its predecessor evidence boundary and acceptance time.\n\n"
        f"Source baseline: `{release['source_baseline']}`. Released at: `{released_at}`.\n\n"
        "## Accepted Claim Set\n"
        "The release materializes this dependency-closed claim set.\n\n"
        f"{claim_lines}\n"
    )


def accepted_claims(data: dict[str, Any]) -> Iterable[dict[str, Any]]:
    """Yield accepted claims in stable identifier order."""

    validate_registry(data)
    claims = (
        claim
        for claim in data["claims"]
        if claim["accepted_in"] is not None and claim["epistemic"] in {"active", "qualified"}
    )
    return sorted(claims, key=lambda claim: claim["id"])


def render_claim_index(data: dict[str, Any]) -> str:
    """Render canonical documentation solely from accepted registry state."""

    lines = [
        "<!-- GENERATED: scripts/render_docs.py; DO NOT EDIT -->",
        "# Accepted claim index",
        "",
        "This document is generated from `governance/claims.yaml`.",
        "",
    ]
    claims = list(accepted_claims(data))
    if not claims:
        lines.append("No scientific claims have been accepted into this repository yet.")
    for claim in claims:
        metadata = [
            f"- Accepted in: `{claim['accepted_in']}`",
            f"- Verification: `{claim['verification']}`",
            f"- Compatibility: `{claim['compatibility']}`",
            f"- Dependencies: {', '.join(claim['dependencies']) or 'none'}",
        ]
        if "category" in claim or "layer" in claim:
            metadata.extend(
                [
                    f"- Category: `{claim.get('category', 'standard')}`",
                    f"- Layer: `{claim.get('layer', 'core')}`",
                ]
            )
        composition = claim.get("composition")
        if isinstance(composition, dict):
            glue = composition["glue"]
            metadata.extend(
                [
                    f"- Structural gap: {composition['structural_gap']}",
                    f"- Glue proof: `{glue['method']}` at `{glue['artifact']}` "
                    f"(`{glue['entrypoint']}`)",
                ]
            )
        hypothesis = claim.get("hypothesis")
        if isinstance(hypothesis, dict):
            metadata.append(
                f"- Hypothesis `{hypothesis['label']}`: {hypothesis['statement']}"
            )
        if "exclusions" in claim:
            metadata.append(
                f"- Exclusions: {', '.join(claim['exclusions']) or 'none'}"
            )
        lines.extend(
            [
                f"## {claim['id']}",
                "",
                str(claim["statement"]),
                "",
                *metadata,
                "",
            ]
        )
    return "\n".join(lines).rstrip() + "\n"
