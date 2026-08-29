"""Thin exact-check aggregator for the P248 campaign."""

from __future__ import annotations

from exact_compatibility_checks import run as run_compatibility
from exact_metric_checks import run as run_metric
from exact_relaxed_checks import run as run_relaxed


def main() -> int:
    run_relaxed()
    run_metric()
    run_compatibility()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
