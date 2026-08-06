from __future__ import annotations

import sympy as sp
import pytest

import substrate_framework as framework
from substrate_framework.dimensional_sine_gordon import dimensional_sine_gordon_coefficients
from substrate_framework.sine_gordon_fiber_source import (
    CanonicalFiberStress,
    canonical_sine_gordon_fiber_stress,
)


def test_public_package_exports_fiber_source_api() -> None:
    assert framework.CanonicalFiberStress is CanonicalFiberStress
    assert (
        framework.canonical_sine_gordon_fiber_stress
        is canonical_sine_gordon_fiber_stress
    )


def test_canonical_fiber_is_an_isolated_conserved_four_source_on_shell() -> None:
    x, t, y, z = sp.symbols("x t y z", real=True)
    field = sp.Function("u")(x, t)
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    source = canonical_sine_gordon_fiber_stress(field, x, t, y, z, coefficients)
    transverse = sp.DiracDelta(y) * sp.DiracDelta(z)
    assert source.transverse_density == transverse
    assert source.stress_energy[:2, :2] == (
        source.sine_gordon_stress.contravariant * transverse
    )
    assert source.stress_energy[2:, :] == sp.zeros(2, 4)
    assert source.stress_energy[:, 2:] == sp.zeros(4, 2)
    assert source.divergence[:2, :] == source.sine_gordon_stress.divergence * transverse
    assert source.divergence[2:, :] == sp.zeros(2, 1)
    assert source.divergence.subs(
        source.sine_gordon_stress.field_equation_residual, 0
    ) == sp.zeros(4, 1)


def test_fiber_embedding_requires_explicit_coordinate_symbols() -> None:
    x, t, y = sp.symbols("x t y", real=True)
    with pytest.raises(ValueError, match="transverse_z"):
        canonical_sine_gordon_fiber_stress(
            sp.Function("u")(x, t),
            x,
            t,
            y,
            0,
            dimensional_sine_gordon_coefficients(1, 1, 1),
        )
