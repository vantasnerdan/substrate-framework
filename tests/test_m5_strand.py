"""Matrix-derived, boundary-sensitive oracles for conditional P254 strand claims."""

from __future__ import annotations

import math

import pytest
import sympy as sp

from substrate_framework.m5_strand import (
    pair_block,
    pair_curvature_xy,
    radial_tension_infimum,
    spectral_potential,
    spectral_strand_infimum,
    trace_director_first_variation,
    trace_potential,
    trace_strand_infimum,
    trace_strand_leading_bound,
)


def test_matrix_curvature_distinguishes_source_norms_and_winding() -> None:
    rho, phi, mean = sp.symbols("rho phi mean", positive=True)
    b = sp.Function("b")(rho)
    for m in (1, 2, 3):
        fxy = pair_curvature_xy(b, rho, phi, m, mean)
        expected = 2 * m * b * b.diff(rho) / rho
        assert sp.simplify(fxy[0, 1] - expected) == 0
        assert sp.simplify(fxy[1, 0] + expected) == 0
        field_norm = sp.trace(fxy * fxy.T)
        assert sp.simplify(2 * field_norm - 16 * m**2 * (b * b.diff(rho) / rho) ** 2) == 0
        assert sp.simplify(4 * field_norm - 32 * m**2 * (b * b.diff(rho) / rho) ** 2) == 0
    # Both conventions cannot use the same coefficient: that error changes T by sqrt(2).
    assert 2 * field_norm != 4 * field_norm


def test_trace_potential_comes_from_four_distinct_matrix_traces() -> None:
    q, q0, mean, weight, g = sp.symbols("q q0 mean weight g", positive=True)
    b, b0 = sp.symbols("b b0", real=True)
    director = sp.symbols("director", real=True)
    pair = pair_block(b, 0, 1, mean)
    matrix = sp.diag(-g, 0, 0, director)
    matrix[1:3, 1:3] = pair
    eta = sp.diag(-1, 1, 1, 1)
    mixed = matrix * eta
    vacuum = [g, mean + b0, mean - b0, 1]
    deltas = [
        sp.expand(sp.trace(mixed**p) - sum(value**p for value in vacuum))
        for p in range(1, 5)
    ]
    assert deltas[0] == director - 1
    source_potential = weight * sum(delta**2 for delta in deltas)
    expected = trace_potential(q, q0, mean, weight)
    frozen = source_potential.subs(director, 1).subs({b**2: q, b0**2: q0})
    assert sp.factor(frozen - expected) == 0
    derivative = sp.diff(source_potential, director).subs(director, 1)
    derivative = derivative.subs({b**2: q, b0**2: q0})
    assert sp.factor(
        derivative - trace_director_first_variation(q, q0, mean, weight)
    ) == 0
    assert trace_director_first_variation(
        sp.Rational(9, 800), sp.Rational(9, 400), sp.Rational(3, 20), 1
    ) < 0
    assert sp.simplify(expected.subs(q, q0)) == 0


def test_director_relaxation_has_no_spatial_curvature_cost() -> None:
    rho, phi, mean = sp.symbols("rho phi mean", positive=True)
    b = sp.Function("b")(rho)
    director = sp.Function("director")(rho)
    field = sp.diag(-8, 0, 0, director)
    field[1:3, 1:3] = pair_block(b, phi, 1, mean)
    eta = sp.diag(-1, 1, 1, 1)
    radial, angular = field.diff(rho), field.diff(phi)
    curvature = (radial * eta * angular - angular * eta * radial) / rho
    frozen = curvature.subs({director: 1, director.diff(rho): 0})
    assert all(sp.simplify(left - right) == 0 for left, right in zip(curvature, frozen))


def test_spectral_tension_uses_core_boundary_and_half_winding() -> None:
    b, b0, beta = sp.symbols("b b0 beta", positive=True)
    q = sp.Symbol("q", nonnegative=True)
    mean = sp.Symbol("mean", real=True)
    pair = pair_block(b, 0, 1, mean)
    exterior = pair_block(b0, 0, 1, mean)
    assert sp.simplify(
        spectral_potential(b**2, b0) - sp.trace((pair - exterior) ** 2)
    ) == 0
    # Independent variable substitution dq=2b db in the square-completion cross term.
    integral = sp.integrate(2 * b * sp.sqrt(2) * (b0 - b), (b, 0, b0))
    for winding in (1, 2):
        expected = 2 * sp.pi * sp.sqrt(16) * winding * integral
        result = spectral_strand_infimum(beta, winding)
        assert sp.simplify(result - expected.subs(b0, beta / 2)) == 0
    assert sp.simplify(spectral_potential(q, beta / 2).subs(q, 0) - beta**2 / 2) == 0
    assert spectral_strand_infimum(sp.Rational(3, 10), 2) == 2 * spectral_strand_infimum(
        sp.Rational(3, 10), 1
    )
    # An unpinned axis q(0)=q0 admits the constant vacuum with T=0 instead.
    regular = radial_tension_infimum(32, 1, (1 - q) ** 2, q, 1).doit()
    assert sp.simplify(regular - sp.pi * sp.sqrt(32)) == 0
    assert regular > 0


def test_trace_variable_stiffness_and_distinct_small_split_power() -> None:
    scipy = pytest.importorskip("scipy.integrate")
    relative_corrections = []
    for delta in (0.03, 0.1, 0.3):
        b0 = delta / 2
        q0 = b0**2
        mean = b0
        weight = 1.0
        stiff = lambda q: 4 + 36 * mean**2 + (12 * mean**2 + 2 * (q + q0)) ** 2
        independent = 2 * math.pi * math.sqrt(32 * weight) * scipy.quad(
            lambda q: (q0 - q) * math.sqrt(stiff(q)), 0.0, q0,
            epsabs=1e-14,
        )[0]
        actual = float(trace_strand_infimum(b0, mean, weight, 1).evalf())
        leading = float(trace_strand_leading_bound(b0, mean, weight, 1))
        vacuum_upper = math.pi * math.sqrt(32 * weight * stiff(q0)) * b0**4
        assert actual == pytest.approx(independent, rel=1e-11)
        assert leading < actual < vacuum_upper
        relative_corrections.append(actual / delta**4 / (math.pi * math.sqrt(2) / 2) - 1)
    assert 0 < relative_corrections[0] < relative_corrections[1] < relative_corrections[2]
    assert relative_corrections[0] < 0.002


def test_winding_two_equality_profile_is_not_differentiable_at_axis() -> None:
    r, c = sp.symbols("r c", positive=True)
    # Square completion yields q~c² rho², b~c rho for *every* winding.
    # At x=+r and x=-r, m=2 has the same diagonal displacement +c r:
    # its two directional slopes in x disagree, unlike the m=1 block.
    for m in (1, 2):
        right = pair_block(c * r, 0, m, 0)[0, 0]
        left = pair_block(c * r, sp.pi, m, 0)[0, 0]
        assert sp.simplify(right - c * r) == 0
        assert sp.simplify(left - (-1) ** m * c * r) == 0
        if m == 2:
            assert sp.simplify(right / r + left / r) == 2 * c
        else:
            assert sp.simplify(right / r + left / r) == 0
