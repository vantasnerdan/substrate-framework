"""Deterministic lexical consumer audits for hash-pinned source trees.

Matches produced here are source-navigation evidence. A lexical occurrence in
code, a comment, or a docstring does not by itself establish a semantic or
scientific dependency.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
import hashlib
from pathlib import Path, PurePosixPath
import re
from typing import Mapping


@dataclass(frozen=True)
class SourceMatch:
    """One matched source file with content provenance."""

    path: str
    sha256: str
    groups: tuple[str, ...]


@dataclass(frozen=True)
class SourceAudit:
    """Deterministic result of one explicitly scoped lexical scan."""

    include_pattern: str
    exclusions: tuple[str, ...]
    token_patterns: tuple[tuple[str, str], ...]
    scanned_file_count: int
    matches: tuple[SourceMatch, ...]

    def paths_for(self, group: str) -> tuple[str, ...]:
        """Return sorted relative paths matching one declared token group."""

        declared = {name for name, _ in self.token_patterns}
        if group not in declared:
            raise KeyError(f"unknown token group: {group}")
        return tuple(match.path for match in self.matches if group in match.groups)

    def paths_with_all(self, *groups: str) -> tuple[str, ...]:
        """Return paths that match every named group."""

        if not groups:
            raise ValueError("at least one token group is required")
        declared = {name for name, _ in self.token_patterns}
        unknown = set(groups) - declared
        if unknown:
            raise KeyError(f"unknown token groups: {sorted(unknown)}")
        required = set(groups)
        return tuple(
            match.path
            for match in self.matches
            if required <= set(match.groups)
        )


@dataclass(frozen=True)
class NumpyTrapezoidCompatibility:
    """AST inventory of current and legacy NumPy trapezoidal-integration names."""

    numpy_aliases: tuple[str, ...]
    direct_legacy_attributes: int
    dynamic_legacy_getattrs: int
    imported_legacy_names: int
    direct_current_attributes: int
    dynamic_current_getattrs: int
    imported_current_names: int
    eager_legacy_default_fallbacks: int

    @property
    def legacy_references(self) -> int:
        """Return all executable references to the removed legacy name."""

        return (
            self.direct_legacy_attributes
            + self.dynamic_legacy_getattrs
            + self.imported_legacy_names
        )

    @property
    def current_references(self) -> int:
        """Return all executable references to the current NumPy name."""

        return (
            self.direct_current_attributes
            + self.dynamic_current_getattrs
            + self.imported_current_names
        )

    @property
    def requires_legacy_alias(self) -> bool:
        """Whether immutable execution needs an explicit legacy-name alias."""

        return self.legacy_references > 0


def audit_numpy_trapezoid_compatibility(
    source: str,
    *,
    filename: str = "<unknown>",
) -> NumpyTrapezoidCompatibility:
    """Detect direct, imported, and dynamic NumPy integration-name access.

    The AST audit catches both ``np.trapz`` and ``getattr(np, "trapz")``.
    It also identifies the eager-default bug in
    ``getattr(np, "trapezoid", getattr(np, "trapz"))``: Python evaluates the
    default argument before calling the outer ``getattr``, so that expression
    still aborts under NumPy versions where the legacy name was removed.
    Comments and docstrings do not count as executable references.
    """

    tree = ast.parse(source, filename=filename)
    aliases = {
        alias.asname or alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
        if alias.name == "numpy"
    }

    def dynamic_name(node: ast.AST, name: str) -> bool:
        return (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "getattr"
            and len(node.args) >= 2
            and isinstance(node.args[0], ast.Name)
            and node.args[0].id in aliases
            and isinstance(node.args[1], ast.Constant)
            and node.args[1].value == name
        )

    def direct_name(node: ast.AST, name: str) -> bool:
        return (
            isinstance(node, ast.Attribute)
            and isinstance(node.value, ast.Name)
            and node.value.id in aliases
            and node.attr == name
        )

    def imported_name(node: ast.AST, name: str) -> int:
        if not isinstance(node, ast.ImportFrom) or node.module != "numpy":
            return 0
        return sum(alias.name == name for alias in node.names)

    def contains_legacy(node: ast.AST) -> bool:
        return any(
            direct_name(descendant, "trapz")
            or dynamic_name(descendant, "trapz")
            for descendant in ast.walk(node)
        )

    nodes = tuple(ast.walk(tree))
    eager_fallbacks = sum(
        dynamic_name(node, "trapezoid")
        and len(node.args) >= 3
        and contains_legacy(node.args[2])
        for node in nodes
    )
    return NumpyTrapezoidCompatibility(
        numpy_aliases=tuple(sorted(aliases)),
        direct_legacy_attributes=sum(direct_name(node, "trapz") for node in nodes),
        dynamic_legacy_getattrs=sum(dynamic_name(node, "trapz") for node in nodes),
        imported_legacy_names=sum(imported_name(node, "trapz") for node in nodes),
        direct_current_attributes=sum(direct_name(node, "trapezoid") for node in nodes),
        dynamic_current_getattrs=sum(dynamic_name(node, "trapezoid") for node in nodes),
        imported_current_names=sum(imported_name(node, "trapezoid") for node in nodes),
        eager_legacy_default_fallbacks=eager_fallbacks,
    )


def _normalized_exclusions(exclusions: tuple[str, ...]) -> tuple[str, ...]:
    normalized: list[str] = []
    for exclusion in exclusions:
        path = PurePosixPath(exclusion)
        if path.is_absolute() or ".." in path.parts or str(path) in {"", "."}:
            raise ValueError(f"invalid relative exclusion: {exclusion!r}")
        normalized.append(path.as_posix().rstrip("/"))
    return tuple(sorted(set(normalized)))


def _is_excluded(relative_path: str, exclusions: tuple[str, ...]) -> bool:
    path = PurePosixPath(relative_path)
    return any(
        path == PurePosixPath(exclusion)
        or PurePosixPath(exclusion) in path.parents
        for exclusion in exclusions
    )


def audit_source_tokens(
    root: Path,
    token_patterns: Mapping[str, str],
    *,
    include_pattern: str = "**/*.py",
    exclusions: tuple[str, ...] = (),
    flags: int = re.IGNORECASE,
) -> SourceAudit:
    """Scan an explicit source root and return sorted, hashed lexical matches.

    ``token_patterns`` maps stable group names to regular-expression strings.
    Exclusions are component-aware relative paths beneath ``root``. Text is
    decoded as UTF-8 with replacement so undecodable bytes cannot silently
    remove a file from the scan.
    """

    source_root = Path(root)
    if not source_root.is_dir():
        raise ValueError(f"source root is not a directory: {source_root}")
    if not include_pattern:
        raise ValueError("include_pattern must be non-empty")
    if not token_patterns:
        raise ValueError("at least one token pattern is required")

    compiled: dict[str, re.Pattern[str]] = {}
    for name, pattern in token_patterns.items():
        if not isinstance(name, str) or not name:
            raise ValueError("token group names must be non-empty strings")
        if not isinstance(pattern, str) or not pattern:
            raise ValueError(f"token pattern for {name!r} must be non-empty")
        compiled[name] = re.compile(pattern, flags)

    normalized_exclusions = _normalized_exclusions(exclusions)
    files: list[tuple[str, Path]] = []
    for path in source_root.glob(include_pattern):
        if not path.is_file():
            continue
        relative = path.relative_to(source_root).as_posix()
        if not _is_excluded(relative, normalized_exclusions):
            files.append((relative, path))
    files.sort(key=lambda item: item[0])

    matches: list[SourceMatch] = []
    for relative, path in files:
        payload = path.read_bytes()
        text = payload.decode("utf-8", errors="replace")
        matched_groups = tuple(
            sorted(name for name, regex in compiled.items() if regex.search(text))
        )
        if matched_groups:
            matches.append(
                SourceMatch(
                    path=relative,
                    sha256=hashlib.sha256(payload).hexdigest(),
                    groups=matched_groups,
                )
            )

    return SourceAudit(
        include_pattern=include_pattern,
        exclusions=normalized_exclusions,
        token_patterns=tuple(sorted(token_patterns.items())),
        scanned_file_count=len(files),
        matches=tuple(matches),
    )
