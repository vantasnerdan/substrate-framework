"""P2540025 grants/source coverage; reuse immutable0024 IO without changing it.

Producer recipes bind source/data/settings, not a renewable execution clock.
The pre-execution dependency universe covers future lazy consumers; actual
module inventories are checked against it before/after the new computation.
This is protocol-immutable source provenance, not whole-OS attestation.
"""
from pathlib import Path
from decimal import Decimal

from substrate_framework.qe_saved_native import read_bound_bytes
from material_pilot_io import (
    binding, canonical, commit_block, recover_block,
    runtime_context as previous_runtime_context,
    implementation as previous_implementation,
)
import hashlib
import json
from material_publication import BUSCTL, manager_runtime_limit

HERE = Path(__file__).resolve().parent
CAPS = {
    "material-reference": 300, "physical-factors": 180, "material-dyad-maps": 90,
    "material-ionic-linear": 90, "material-bloch-grid": 240,
    "material-neutral-preparation": 120,
}


def runtime_context():
    context = previous_runtime_context()
    context['manager_runtime_max_usec'] = manager_runtime_limit(context)
    return context


def bound_json(record):
    raw, _ = read_bound_bytes(record, max_bytes=16 << 20)
    return json.loads(raw)


def implementation(extra_paths=()):
    records = {item["path"]: item for item in previous_implementation()}
    for path in (Path(__file__), BUSCTL, *extra_paths):
        item = binding(path)
        records[item["path"]] = item
    return [records[key] for key in sorted(records)]


def verify_consumers(grant, sources):
    """Compare actual consumers to bindings that existed before execution."""
    universe = bound_json(grant["dependency_inventory"])
    expected = {item["path"]: item for item in universe["files"]}
    declared = grant["implementation_sources"]
    if len({item["path"] for item in declared}) != len(declared):
        raise ValueError("duplicate declared implementation source")
    expected.update({item["path"]: item for item in declared})
    for item in sources:
        if expected.get(item["path"]) != item:
            raise ValueError("unbound actual consumer: " + item["path"])


def _context(grant, context):
    for key in ("uid", "gid", "groups", "cgroup", "cpu_affinity"):
        if grant["context"][key] != context[key]:
            raise ValueError("actual manager context differs: " + key)
    if not context["INVOCATION_ID"] or "/user.slice/" not in context["cgroup"]:
        raise ValueError("actual user-manager invocation required")
    if len(context["cpu_affinity"]) != 1:
        raise ValueError("this grant requires one pinned CPU")
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        if context["thread_variables"][key] != "1":
            raise ValueError("single-thread numerical runtime required: " + key)
    for key, cap in (("memory.max", 4294967296), ("memory.swap.max", 0)):
        field = context["resources"][key]
        if field["availability"] != "available" or not field["value"].isdigit() or int(field["value"]) > cap:
            raise ValueError("unadmitted kernel resource: " + key)
    # cpuset.cpus.effective is optional and stays explicitly typed unavailable.
    runtime = context['manager_runtime_max_usec']
    exact_us = int(Decimal(str(grant['remaining_seconds'])) * 1000000)
    if (runtime['availability'] != 'available'
            or not 0 < int(runtime['value']) <= exact_us):
        raise ValueError('own manager deadline exceeds exact invocation allowance')


def admit(path, kind, sources, context, context_only=False):
    """Return grant, immutable recipe, identity. Context-only reads no material.

    Required authority and input/settings fields are documented in
    producer-source-contract.json. Production uses a fresh explicit time grant;
    recovered scientific boundaries never confer a new time allowance.
    """
    if kind not in CAPS:
        raise ValueError("unknown producer kind")
    if context_only:
        return {}, {}, None
    grant = json.loads(Path(path).read_bytes())
    registration = bound_json(grant["registration"])
    if Path(grant["registration"]["path"]).resolve() != HERE / "implementation-registration.json":
        raise ValueError("wrong implementation registration")
    if grant["source_access"] != registration["source_access"]:
        raise ValueError("source access differs from registered authority")
    bound_json(grant["source_access"])
    for item in grant["implementation_sources"]:
        if binding(item["path"]) != item:
            raise ValueError("declared implementation source changed: " + item["path"])
    verify_consumers(grant, sources)
    _context(grant, context)
    for key, cap in (("remaining_seconds", CAPS[kind]), ("first_block_seconds", 90)):
        if type(grant[key]) not in (int, float) or not 0 < grant[key] <= cap:
            raise ValueError("explicit bounded " + key + " required")
    if type(grant["output_max_bytes"]) is not int or not 0 < grant["output_max_bytes"] <= 1073741824:
        raise ValueError("invalid output byte cap")
    if type(grant["resume"]) is not bool:
        raise ValueError("explicit resume flag required")
    if not isinstance(grant["inputs"], dict) or not isinstance(grant["settings"], dict):
        raise ValueError("explicit producer inputs/settings required")
    recipe = {"kind": kind, "registration": grant["registration"],
              "source_access": grant["source_access"],
              "dependency_inventory": grant["dependency_inventory"],
              "implementation_sources": grant["implementation_sources"],
              "inputs": grant["inputs"], "settings": grant["settings"],
              "library_binding": grant.get("library_binding")}
    identity = hashlib.sha256(canonical(recipe)).hexdigest()
    if grant["resume"]:
        if grant["resume_identity"] != identity:
            raise ValueError("recovery source/input/recipe identity mismatch")
    elif grant["resume_identity"] is not None:
        raise ValueError("fresh execution cannot name a recovery identity")
    return grant, recipe, identity
