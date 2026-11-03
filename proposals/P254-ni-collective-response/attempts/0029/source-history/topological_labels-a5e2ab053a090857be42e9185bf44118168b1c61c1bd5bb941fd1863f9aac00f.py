"""Exact algebra for integer winding labels.

The parity label defined here is a group character. Calling it a fermion or
boson statistic requires an additional physical spin-statistics map.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Any, Sequence

import sympy as sp


def _integer(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_integer is not True:
        raise ValueError(f"{name} must be an integer")
    return expression


def winding_parity(winding: Any) -> sp.Expr:
    """Return the multiplicative label ``(-1)**winding`` for integer winding."""

    value = _integer(winding, "winding")
    return sp.simplify((-1) ** value)


def combined_winding_parity(anchor_winding: Any, dressing_winding: Any) -> sp.Expr:
    """Return parity after adding integer anchor and dressing windings."""

    anchor = _integer(anchor_winding, "anchor_winding")
    dressing = _integer(dressing_winding, "dressing_winding")
    return winding_parity(anchor + dressing)


def _positive_concrete_integer(value: Any, name: str) -> int:
    expression = _integer(value, name)
    if expression.is_number is not True or expression.is_positive is not True:
        raise ValueError(f"{name} must be a positive concrete integer")
    return int(expression)


@dataclass(frozen=True)
class CyclicSignCharacterLedger:
    """Complete homomorphism ledger from additive ``C_order`` to signs."""

    order: int
    generator_images: tuple[int, ...]
    nontrivial_exists: bool
    nontrivial_kernel: tuple[int, ...] | None
    nontrivial_quotient_order: int | None
    nontrivial_faithful: bool


@dataclass(frozen=True)
class BinaryProductCharacterLedger:
    """Character-count and faithfulness ledger for ``C_2**rank``."""

    rank: int
    character_count: int
    nontrivial_character_count: int
    faithful_character_exists: bool


def cyclic_sign_character_ledger(order: Any) -> CyclicSignCharacterLedger:
    """Classify every homomorphism ``C_order -> {+1, -1}``.

    A homomorphism is fixed by the image ``epsilon`` of the additive
    generator, subject to ``epsilon**order == 1``.  Thus ``-1`` is allowed for
    every even order, not only order two.  Its kernel is the even residues and
    the character is faithful only at order two.
    """

    modulus = _positive_concrete_integer(order, "order")
    nontrivial_exists = modulus % 2 == 0
    kernel = tuple(range(0, modulus, 2)) if nontrivial_exists else None
    return CyclicSignCharacterLedger(
        order=modulus,
        generator_images=(1, -1) if nontrivial_exists else (1,),
        nontrivial_exists=nontrivial_exists,
        nontrivial_kernel=kernel,
        nontrivial_quotient_order=2 if nontrivial_exists else None,
        nontrivial_faithful=nontrivial_exists and kernel == (0,),
    )


def cyclic_sign_character(
    order: Any,
    element: Any,
    *,
    generator_image: Any,
) -> sp.Expr:
    """Evaluate a declared sign character of a finite additive cyclic group."""

    modulus = _positive_concrete_integer(order, "order")
    residue = _integer(element, "element")
    image = sp.sympify(generator_image)
    if image not in (sp.Integer(1), sp.Integer(-1)):
        raise ValueError("generator_image must be +1 or -1")
    if image**modulus != 1:
        raise ValueError("generator_image does not satisfy the cyclic relation")
    return sp.simplify(image ** sp.Mod(residue, modulus))


def _binary_vector(values: Sequence[Any], name: str) -> tuple[int, ...]:
    result = tuple(values)
    if not result:
        raise ValueError(f"{name} must be non-empty")
    parsed: list[int] = []
    for value in result:
        expression = sp.sympify(value)
        if expression not in (sp.Integer(0), sp.Integer(1)):
            raise ValueError(f"{name} entries must be zero or one")
        parsed.append(int(expression))
    return tuple(parsed)


def binary_product_character(
    coordinates: Sequence[Any],
    selector: Sequence[Any],
) -> sp.Integer:
    """Evaluate ``(-1)**(selector dot coordinates)`` on a binary product."""

    point = _binary_vector(coordinates, "coordinates")
    character = _binary_vector(selector, "selector")
    if len(point) != len(character):
        raise ValueError("coordinates and selector must have equal length")
    exponent = sum(left * right for left, right in zip(point, character)) % 2
    return sp.Integer(-1 if exponent else 1)


def binary_product_character_kernel(
    selector: Sequence[Any],
) -> tuple[tuple[int, ...], ...]:
    """Enumerate the exact kernel of one character of ``C_2**rank``."""

    character = _binary_vector(selector, "selector")
    points = product((0, 1), repeat=len(character))
    return tuple(
        point for point in points if binary_product_character(point, character) == 1
    )


def binary_product_character_ledger(rank: Any) -> BinaryProductCharacterLedger:
    """Count all sign characters of ``C_2**rank`` and classify faithfulness."""

    dimension = _positive_concrete_integer(rank, "rank")
    count = 2**dimension
    return BinaryProductCharacterLedger(
        rank=dimension,
        character_count=count,
        nontrivial_character_count=count - 1,
        faithful_character_exists=dimension == 1,
    )
