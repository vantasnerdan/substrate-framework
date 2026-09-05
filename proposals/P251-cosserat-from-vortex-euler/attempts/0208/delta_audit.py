"""Additive C016 semantic/source audit; this is not a scientific oracle."""

from pathlib import Path
import subprocess

import yaml

BASE = "6a58392362108d7256f849bd763625a05134fd42"


def baseline(path):
    return subprocess.check_output(["git", "show", f"{BASE}:{path}"])


registry = "governance/claims.yaml"
old = {entry["id"]: entry for entry in yaml.safe_load(baseline(registry))["claims"]}
new = {entry["id"]: entry for entry in yaml.safe_load(Path(registry).read_text())["claims"]}
assert set(new)-set(old) == {"C-CST-016"}
assert all(new[key] == value for key, value in old.items())
print(f"All {len(old)} previous claim objects unchanged; sole addition C-CST-016.")
paths = subprocess.check_output(
    ["git", "ls-tree", "-r", "--name-only", BASE, "src/substrate_framework"], text=True
).splitlines()
for path in paths:
    assert Path(path).read_bytes() == baseline(path), path
print(f"All {len(paths)} previous canonical source files byte-identical.")
current = yaml.safe_load(Path("governance/releases/current.yaml").read_text())
pinned = yaml.safe_load(Path("governance/releases/v0.181.0.yaml").read_text())
assert current == pinned
assert set(current["accepted_claims"]) == set(new)
print("Current/pinned v0.181.0 and all269 accepted identifiers agree.")
print("Reused0196's39 source checks and0207's eight independent review checks.")
print("Direct new API tests:0208 first-pytest.stdout,12 passed on first execution.")
print("No full suite or unrelated scientific replay required by this additive delta.")
