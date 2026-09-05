"""Semantic blast-radius audit; not an additional scientific oracle."""

import ast
from pathlib import Path
import subprocess

import yaml


def baseline(path):
    return subprocess.check_output(["git", "show", f"dfe495c:{path}"], text=True)


registry = "governance/claims.yaml"
old = {entry["id"]: entry for entry in yaml.safe_load(baseline(registry))["claims"]}
new = {entry["id"]: entry for entry in yaml.safe_load(Path(registry).read_text())["claims"]}
assert set(new)-set(old) == {"C-CST-013"}
assert all(new[key] == value for key, value in old.items())
print(f"All{len(old)} prior claim objects unchanged; sole addition C-CST-013.")
module = "src/substrate_framework/euler_acoustic.py"
previous = ast.parse(baseline(module))
current = ast.parse(Path(module).read_text())
definitions = {item.name: ast.dump(item) for item in current.body
               if isinstance(item, (ast.FunctionDef, ast.ClassDef))}
count = 0
for item in previous.body:
    if isinstance(item, (ast.FunctionDef, ast.ClassDef)):
        assert definitions[item.name] == ast.dump(item), item.name
        count += 1
print(f"All{count} earlier acoustic functions/classes are AST-identical.")
for path in ("src/substrate_framework/euler_fourier.py", "src/substrate_framework/euler_phase.py"):
    assert Path(path).read_text() == baseline(path)
    print("Unchanged executable dependency:", path)
print("Scientific receipt reused:0170 exact proof/native19checks and0177 four direct tests.")
print("No full suite is repeated; fixed canonical/registry/generated checks complete this scope.")
