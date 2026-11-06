"""Shared verification primitives for campaign and framework code."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Iterable


class CheckFailure(AssertionError):
    """Raised when a verification condition fails."""


class SuccessfulCheckTally(int):
    """Backward-compatible verifier success token with a visible check count.

    Historical campaign entry points return ``ledger.finish()`` directly from
    ``main`` and pass that value to :class:`SystemExit`.  Returning the positive
    check count as an ordinary integer made every successful verifier exit with
    a nonzero process status.  This token has integer value zero for the OS,
    while its string/format representation remains the actual tally used by
    older campaign summary lines.  New code should inspect ``passed_count``
    rather than coercing the token arithmetically.
    """

    passed_count: int

    def __new__(cls, passed_count: int) -> "SuccessfulCheckTally":
        instance = int.__new__(cls, 0)
        instance.passed_count = int(passed_count)
        return instance

    def __str__(self) -> str:
        return str(self.passed_count)

    def __repr__(self) -> str:
        return f"SuccessfulCheckTally(passed_count={self.passed_count})"

    def __format__(self, format_spec: str) -> str:
        return format(self.passed_count, format_spec)


@dataclass
class CheckLedger:
    """Collect named checks and emit the repository's terminal tally."""

    claim_id: str
    passed: list[str] = field(default_factory=list)

    def check(self, name: str, condition: object, detail: str = "") -> None:
        if not bool(condition):
            suffix = f" -- {detail}" if detail else ""
            raise CheckFailure(f"FAIL [{self.claim_id}] {name}{suffix}")
        self.passed.append(name)
        print(f"PASS [{self.claim_id}] {name}")

    def mutation_sensitive(
        self,
        name: str,
        predicate: Callable[[object], object],
        baseline: object,
        mutations: Iterable[object],
    ) -> None:
        """Require a predicate to accept the derivation and reject every mutation."""

        mutated = list(mutations)
        self.check(f"{name}: baseline", predicate(baseline), "baseline was rejected")
        self.check(f"{name}: mutations supplied", bool(mutated), "no mutations were tested")
        for index, candidate in enumerate(mutated, start=1):
            self.check(
                f"{name}: mutation {index} rejected",
                not bool(predicate(candidate)),
                "the check is insensitive to a changed input",
            )

    def finish(self) -> SuccessfulCheckTally:
        """Print the terminal tally and return a status-zero visible tally."""

        if not self.passed:
            raise CheckFailure(f"FAIL [{self.claim_id}] no checks ran")
        print(f"ALL {len(self.passed)} CHECKS PASS [{self.claim_id}]")
        return SuccessfulCheckTally(len(self.passed))
