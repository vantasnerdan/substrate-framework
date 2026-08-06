"""Exact real-scalar input contracts shared by symbolic gravity modules."""

from __future__ import annotations

from typing import Any

import sympy as sp


def exact_real(value: Any, name: str) -> sp.Expr:
    """Return an exact provably real SymPy expression."""

    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact rather than floating")
    if expression.is_real is not True:
        raise ValueError(f"{name} must be provably real")
    return expression


def positive_exact(value: Any, name: str) -> sp.Expr:
    """Return an exact provably positive real SymPy expression."""

    expression = exact_real(value, name)
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be provably positive")
    return expression
