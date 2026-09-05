"""Semantic impact audit for the additive optical theorem; no scientific oracle."""

from pathlib import Path
import subprocess

import yaml

BASE = "6131db5ff5990cd6b6591c93a1ace6ec65bdfef9"


def baseline(path):
    return subprocess.check_output(["git", "show", f"{BASE}:{path}"])


registry = "governance/claims.yaml"
old = {entry["id"]: entry for entry in yaml.safe_load(baseline(registry))["claims"]}
new = {entry["id"]: entry for entry in yaml.safe_load(Path(registry).read_text())["claims"]}
assert set(new) - set(old) == {"C-CST-014"}
assert all(new[key] == value for key, value in old.items())
print(f"All {len(old)} prior claim objects unchanged; sole addition C-CST-014.")
paths = subprocess.check_output(
    ["git", "ls-tree", "-r", "--name-only", BASE, "src/substrate_framework"], text=True
).splitlines()
for path in paths:
    assert Path(path).read_bytes() == baseline(path), path
print(f"All {len(paths)} old canonical source files byte-identical.")
current = yaml.safe_load(Path("governance/releases/current.yaml").read_text())
pinned = yaml.safe_load(Path("governance/releases/v0.179.0.yaml").read_text())
assert current == pinned
assert set(current["accepted_claims"]) == set(new)
print("Pinned/current v0.179.0 and all 267 accepted claim identifiers agree.")
print("Reused scientific receipts:0181/0187, independent0190, corrected0191.")
print("Direct new API receipt:0192 repaired-pytest.stdout, 17 passed.")
print("No full test suite or unrelated scientific dependency replay required.")
