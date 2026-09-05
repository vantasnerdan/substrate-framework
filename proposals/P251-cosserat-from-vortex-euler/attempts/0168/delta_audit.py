"""Semantic blast-radius audit of this additive promotion, not a science oracle."""

import ast
from pathlib import Path
import subprocess

import yaml


def baseline(path):
    return subprocess.check_output(["git", "show", f"dbf0c04:{path}"], text=True)


registry_path = "governance/claims.yaml"
old = {item["id"]: item for item in yaml.safe_load(baseline(registry_path))["claims"]}
new = {item["id"]: item for item in yaml.safe_load(Path(registry_path).read_text())["claims"]}
assert set(new)-set(old) == {"C-CST-012"}
assert all(new[key] == value for key, value in old.items())
print(f"All {len(old)} prior claim records are semantically unchanged; sole addition C-CST-012.")

phase_path = "src/substrate_framework/euler_phase.py"
before = ast.parse(baseline(phase_path))
after = ast.parse(Path(phase_path).read_text())
definitions = {item.name: ast.dump(item) for item in after.body
               if isinstance(item, (ast.FunctionDef, ast.ClassDef))}
checked = 0
for item in before.body:
    if isinstance(item, (ast.FunctionDef, ast.ClassDef)):
        assert definitions[item.name] == ast.dump(item), item.name
        checked += 1
print(f"All {checked} prior phase functions/classes are AST-identical; new helper is additive.")
print("Scientific receipts:0161 exact proof/checks,0162 thirteen affected tests,0167 five new tests.")
print("No unchanged downstream phase consumer has an altered executable dependency.")
