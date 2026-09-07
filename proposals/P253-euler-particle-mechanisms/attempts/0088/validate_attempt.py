#!/usr/bin/env python3
"""Static claim-boundary validation for P253/0088."""

from hashlib import sha256
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parent


def digest(name: str) -> str:
    return sha256((ROOT / name).read_bytes()).hexdigest()


def main() -> None:
    assert digest("README.md") == (
        "6b01cbe08fc0321610b64315725e46c6c582be7592939d58b9795a44c5a65393"
    )
    assert (ROOT / "activation-schema.exit").read_text().strip() == "0"
    assert (ROOT / "hessian-oracle-dual-riesz.exit").read_text().strip() == "0"
    assert not (ROOT / "hessian-oracle-dual-riesz.stderr.txt").read_text()

    derivation = (ROOT / "derivation.md").read_text()
    required = (
        "rho_0 integral zeta r^2",
        "rho_0 integral r zeta",
        "L_1^k=C+2k^2 s c",
        "P_row^* C_0^*g",
        "first permitted response order is `O(delta)`",
        "source-derived smooth-core",
        "gamma_12^col>0",
        "ell_1=L_St",
        "A_1^sharp=rho_0 A_1!=0",
        "No density of compact",
        "M_12(h)=integral w dot F_12^v",
        "M_(12,N)(h_N)=M_12^col(h_col)+o(1)",
        "No shrinking-",
        "Im(c_1 conjugate(c_2))!=0",
    )
    for phrase in required:
        assert phrase in derivation, phrase
    forbidden = (
        "response starts at order `delta`",
        "with exact leading loss `delta^1`",
        "P_S C_0^*g",
        "P_row^*F_12^v",
        "D_curv[q]!=0",
        "author-stage 0083",
        "0083_path_pending_review",
    )
    for phrase in forbidden:
        assert phrase not in derivation, phrase

    result = yaml.safe_load((ROOT / "result.yaml").read_text())
    assert result["route_verdicts"]["route_B_Hessian_Lagrange_row"]["verdict"] == "established"
    assert result["route_verdicts"]["route_B_zero_zero_diagonal_seed"]["verdict"] == "refuted"
    assert result["route_verdicts"]["route_B_curvature_mixed_diagonal_seed"]["verdict"] == "blocked"
    assert result["route_verdicts"]["route_B_offdiagonal_response"]["verdict"] == "established"
    assert result["parent_campaign_state"] == "active"
    assert result["exhaustion_claimed"] is False

    bad = []
    for path in ROOT.iterdir():
        if not path.is_file():
            continue
        for line_number, line in enumerate(
            path.read_text(errors="ignore").splitlines(), start=1
        ):
            if line.rstrip() != line:
                bad.append((path.name, line_number))
    assert not bad, bad

    print("PASS replayed README hash and activation exit")
    print("PASS physical r^2 KKS and residual fixed-fiber r row")
    print("PASS fixed-k massive curvature and graph-C1 claim boundary")
    print("PASS typed finite-row displacement and stabilizer verdict")
    print("PASS dual-Riesz/Sturm closure and fixed-column near-axis response")
    print("PASS reviewed 0083/0089 source-specific Cao transfer boundary")
    print("PASS route verdicts remain response-scoped and campaign active")
    print("PASS owned-file trailing-whitespace scan")


if __name__ == "__main__":
    main()
