"""Regression checks for public contribution and workflow safety surfaces."""

from __future__ import annotations

from pathlib import Path
import re
import tomllib

import yaml


ROOT = Path(__file__).resolve().parents[1]


def test_required_community_files_exist_and_are_linked() -> None:
    required = (
        "CODE_OF_CONDUCT.md",
        "CONTRIBUTING.md",
        "SECURITY.md",
        "SUPPORT.md",
        "LICENSE",
        "THIRD_PARTY_NOTICES.md",
        ".github/CODEOWNERS",
        ".github/dependabot.yml",
        ".github/workflows/validate.yml",
    )
    for relative_path in required:
        assert (ROOT / relative_path).is_file(), relative_path

    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    for relative_path in (
        "CODE_OF_CONDUCT.md",
        "CONTRIBUTING.md",
        "SECURITY.md",
        "LICENSE",
        "THIRD_PARTY_NOTICES.md",
    ):
        assert f"]({relative_path})" in readme


def test_license_metadata_and_exceptions_are_explicit() -> None:
    metadata = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    project = metadata["project"]
    assert project["license"] == "Apache-2.0"
    assert project["license-files"] == ["LICENSE", "LICENSES/*"]

    license_text = (ROOT / "LICENSE").read_text(encoding="utf-8")
    assert "Apache License" in license_text
    assert "Version 2.0, January 2004" in license_text

    notices = (ROOT / "THIRD_PARTY_NOTICES.md").read_text(encoding="utf-8")
    assert "incoming/einbein_1plus1D_tutorial.pdf" in notices
    assert "docs/tutorials/einbein_2plus1D/" in notices
    assert "docs/tutorials/einbein_3plus1D/" in notices
    assert "GPL-3.0-or-later" in notices
    assert "CODE_OF_CONDUCT.md" in notices
    assert "CC-BY-SA-4.0" in notices


def test_all_files_require_designated_code_owner_approval() -> None:
    codeowners = (ROOT / ".github/CODEOWNERS").read_text(encoding="utf-8")
    active_lines = [
        line for line in codeowners.splitlines() if line and not line.startswith("#")
    ]
    assert active_lines == ["* @vantasnerdan @axis-marbell @mlops-kelvin"]

    contributing = (ROOT / "CONTRIBUTING.md").read_text(encoding="utf-8")
    for login in ("@vantasnerdan", "@axis-marbell", "@mlops-kelvin"):
        assert login in contributing
    assert "do not receive merge access" in contributing
    assert "automatically after merge" in contributing


def test_pull_request_workflow_is_least_privilege_and_sha_pinned() -> None:
    workflow = (ROOT / ".github/workflows/validate.yml").read_text(encoding="utf-8")
    assert "pull_request_target:" not in workflow
    assert "secrets." not in workflow
    assert re.search(r"(?m)^permissions:\n  contents: read$", workflow)
    assert "persist-credentials: false" in workflow
    assert "scripts/validate.sh --full" in workflow
    assert "timeout-minutes:" in workflow

    uses_lines = [line.strip() for line in workflow.splitlines() if "uses:" in line]
    assert uses_lines
    for line in uses_lines:
        reference = line.split("#", 1)[0].strip()
        assert re.fullmatch(r"uses: [A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+@[0-9a-f]{40}", reference)


def test_dependabot_covers_python_and_actions() -> None:
    config = yaml.safe_load((ROOT / ".github/dependabot.yml").read_text(encoding="utf-8"))
    assert config["version"] == 2
    ecosystems = {entry["package-ecosystem"] for entry in config["updates"]}
    assert ecosystems == {"github-actions", "pip"}
    for entry in config["updates"]:
        assert entry["directory"] == "/"
        assert entry["schedule"]["interval"] == "weekly"


def test_issue_forms_preserve_issue_first_and_rights_boundaries() -> None:
    issue_directory = ROOT / ".github/ISSUE_TEMPLATE"
    config = yaml.safe_load((issue_directory / "config.yml").read_text(encoding="utf-8"))
    assert config["blank_issues_enabled"] is False

    for name in ("agent_task.yml", "bug_report.yml", "documentation.yml"):
        form = yaml.safe_load((issue_directory / name).read_text(encoding="utf-8"))
        assert form["name"]
        assert form["description"]
        serialized = (issue_directory / name).read_text(encoding="utf-8")
        assert "pull request" in serialized.lower()
        assert "issue" in serialized.lower()


def test_contribution_policy_blocks_sensitive_and_unlicensed_sources() -> None:
    contributing = (ROOT / "CONTRIBUTING.md").read_text(encoding="utf-8").lower()
    assert "credentials" in contributing
    assert "paywalled" in contributing
    assert "redistribution" in contributing
    assert "must not merge" in contributing


def test_rights_restricted_preparata_sources_are_not_distributed() -> None:
    source_directory = (
        ROOT / "proposals/P229-preparata-qcd-vacuum-audit/sources"
    )
    restricted = (
        "preparata1986-nuovo-cim-a96-366.pdf",
        "preparata1986-extracted.txt",
    )
    for name in restricted:
        assert not (source_directory / name).exists()

    ignore_rules = (ROOT / ".gitignore").read_text(encoding="utf-8")
    for name in restricted:
        path = f"/proposals/P229-preparata-qcd-vacuum-audit/sources/{name}"
        assert path in ignore_rules

    source_notice = (source_directory / "README.md").read_text(encoding="utf-8")
    assert "10.1007/BF02833896" in source_notice
    assert "d712d4a085582ac390701e0a1f43bf21796f0971eda9f355c48da8b0d98b5b8a" in source_notice
    assert "2cb0a6ca8b929d3e1513a27351173cbe05b354c1bcbc878ccaee77092c2dc961" in source_notice
