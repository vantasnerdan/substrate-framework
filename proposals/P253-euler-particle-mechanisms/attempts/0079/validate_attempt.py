"""Bounded artifact-consistency checks for P253/0079."""

from __future__ import annotations

import hashlib
from pathlib import Path

import yaml


HERE = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


assert (HERE / "activation-schema.exit").read_text() == "0\n"
assert sha256(HERE / "README.md") == (
    "9a61ba2d5df5677540924392453149eadab760b02eb75584563c24fe6d210942"
)

result = yaml.safe_load((HERE / "result.yaml").read_text())
assert result["attempt"] == "P253/0079"
assert result["parent_campaign_state"] == "active"
assert result["exhaustion_claimed"] is False
assert {entry["verdict"] for entry in result["route_verdicts"].values()} <= {
    "established",
    "refuted",
    "blocked",
}
assert result["route_verdicts"]["route_A_symmetry_enforced_multiplicity"][
    "verdict"
] == "refuted"
assert result["route_verdicts"]["route_B_distinct_sector_crossing_and_controls"][
    "verdict"
] == "blocked"
assert "o(1/N)" in result["exact_results"]["rational_ray_crossing_reduction"][
    "statement"
]
assert "return" in result["exact_results"]["coupled_gate_ledger"]["statement"]
assert "absolute" in result["exact_results"]["coupled_gate_ledger"][
    "statement"
]
assert (HERE / "symbolic-v4.exit").read_text() == "0\n"

print("PASS frozen README and activation receipt")
print("PASS route verdict vocabulary and active parent scope")
print("PASS result records C0 rational-ray bracket and two-sided return kernel")
print("PASS corrected symbolic oracle receipt")
