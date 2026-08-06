from __future__ import annotations

import pytest
import sympy as sp

from substrate_framework.exact_symbolic import exact_real, positive_exact


def test_exact_real_accepts_exact_provably_real_values() -> None:
    real_symbol = sp.Symbol("x", real=True)
    assert exact_real(real_symbol, "x") == real_symbol
    assert exact_real(sp.Rational(-3, 5), "ratio") == sp.Rational(-3, 5)


def test_exact_real_rejects_floats_and_nonreal_values() -> None:
    with pytest.raises(ValueError, match="exact rather than floating"):
        exact_real(sp.Float("1.25"), "value")
    with pytest.raises(ValueError, match="provably real"):
        exact_real(sp.I, "value")


def test_positive_exact_requires_a_provably_positive_value() -> None:
    positive_symbol = sp.Symbol("a", positive=True)
    assert positive_exact(positive_symbol, "a") == positive_symbol
    with pytest.raises(ValueError, match="provably positive"):
        positive_exact(sp.Integer(0), "value")
