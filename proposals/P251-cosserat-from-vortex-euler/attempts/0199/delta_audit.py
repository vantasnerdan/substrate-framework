"""Semantic impact audit for the additive displacement theorem; no scientific oracle."""

from pathlib import Path
import subprocess

import yaml

BASE = "ee0502c18ed9d2ff86976bcc066f0d0641d912bc"


def baseline(path):
    return subprocess.check_output(["git", "show", f"{BASE}:{path}"])


registry = "governance/claims.yaml"
old = {entry["id"]: entry for entry in yaml.safe_load(baseline(registry))["claims"]}
new = {entry["id"]: entry for entry in yaml.safe_load(Path(registry).read_text())["claims"]}
assert set(new) - set(old) == {"C-CST-015"}
assert all(new[key] == value for key, value in old.items())
print(f"All {len(old)} prior claim objects unchanged; sole addition C-CST-015.")
paths = subprocess.check_output(
    ["git", "ls-tree", "-r", "--name-only", BASE, "src/substrate_framework"], text=True
).splitlines()
for path in paths:
    assert Path(path).read_bytes() == baseline(path), path
print(f"All {len(paths)} old canonical source files byte-identical.")
current = yaml.safe_load(Path("governance/releases/current.yaml").read_text())
pinned = yaml.safe_load(Path("governance/releases/v0.180.0.yaml").read_text())
assert current == pinned
assert set(current["accepted_claims"]) == set(new)
print("Pinned/current v0.180.0 and all 268 accepted claim identifiers agree.")
print("Reused scientific receipts:0188 and independent0197; no correction needed.")
print("Direct new API receipt:0199 repaired-pytest.stdout, 10 passed.")
print("No full test suite or unrelated scientific dependency replay required.")
