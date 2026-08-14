"""Exact and mutation-sensitive checks for the conditional flowing metric."""

from __future__ import annotations

import pytest
import sympy as sp

from substrate_framework.optical_geometry import (
    optical_metric_1d,
    optical_ricci_scalar_1d,
)
from substrate_framework.substrate_metric import (
    FlowingSubstrateMetric,
    flowing_christoffel_lab_frame,
    flowing_metric_ricci_scalar,
    flowing_substrate_metric,
)


def _independent_christoffels(
    metric: sp.Matrix,
    coordinates: tuple[sp.Symbol, sp.Symbol],
) -> tuple[tuple[tuple[sp.Expr, ...], ...], ...]:
    """Reconstruct ``Gamma^a_bc`` without the package geometry helpers."""

    inverse = metric.inv()
    dimension = len(coordinates)
    return tuple(
        tuple(
            tuple(
                sp.simplify(
                    sum(
                        inverse[upper, inner]
                        * (
                            sp.diff(metric[inner, second], coordinates[first])
                            + sp.diff(metric[inner, first], coordinates[second])
                            - sp.diff(metric[first, second], coordinates[inner])
                        )
                        for inner in range(dimension)
                    )
                    / 2
                )
                for second in range(dimension)
            )
            for first in range(dimension)
        )
        for upper in range(dimension)
    )


def test_declared_metric_inverse_determinant_and_null_roots() -> None:
    n, c0 = sp.symbols("n c0", positive=True)
    velocity = sp.symbols("V", real=True)
    result = flowing_substrate_metric(n, velocity, c0)
    expected_inverse = sp.Matrix(
        [
            [-n / c0**2, n * velocity / c0**2],
            [n * velocity / c0**2, 1 / n - n * velocity**2 / c0**2],
        ]
    )
    assert isinstance(result, FlowingSubstrateMetric)
    assert result.inverse == expected_inverse
    assert (result.inverse * result.covariant).applyfunc(sp.simplify) == sp.eye(2)
    assert result.determinant == -c0**2

    local_speed = sp.symbols("v_local", real=True)
    null_form = sp.expand(
        result.covariant[0, 0]
        + 2 * result.covariant[0, 1] * local_speed
        + result.covariant[1, 1] * local_speed**2
    )
    independently_solved = tuple(
        sorted(
            (sp.simplify(root) for root in sp.solve(null_form, local_speed)),
            key=str,
        )
    )
    assert tuple(sorted(result.massless_speeds, key=str)) == independently_solved
    assert set(result.massless_speeds) == {
        -velocity - c0 / n,
        -velocity + c0 / n,
    }


def test_static_limit_is_constant_rescaling_of_accepted_metric() -> None:
    t, x = sp.symbols("t x", real=True)
    c0 = sp.symbols("c0", positive=True)
    n = sp.Function("n", positive=True)(x)
    result = flowing_substrate_metric(n, 0, c0)
    accepted = optical_metric_1d(n, c0)
    assert result.covariant == c0**2 * accepted

    source_curvature = flowing_metric_ricci_scalar(n, 0, c0, t, x)
    accepted_curvature = optical_ricci_scalar_1d(n, x, c0)
    assert sp.simplify(source_curvature - accepted_curvature / c0**2) == 0


def test_minkowski_limit() -> None:
    c0 = sp.symbols("c0", positive=True)
    result = flowing_substrate_metric(1, 0, c0)
    assert result.covariant == sp.diag(-c0**2, 1)
    assert result.inverse == sp.diag(-1 / c0**2, 1)
    assert result.massless_speeds == (-c0, c0)


def test_general_christoffels_match_independent_reconstruction() -> None:
    t, x = sp.symbols("t x", real=True)
    c0 = sp.symbols("c0", positive=True)
    n = sp.Function("n", positive=True)(t, x)
    velocity = sp.Function("V", real=True)(t, x)
    metric = flowing_substrate_metric(n, velocity, c0).covariant
    direct = _independent_christoffels(metric, (t, x))
    returned = flowing_christoffel_lab_frame(n, velocity, c0, t, x)
    expected = {
        "Gamma^t_tt": direct[0][0][0],
        "Gamma^t_tX": direct[0][0][1],
        "Gamma^t_XX": direct[0][1][1],
        "Gamma^X_tt": direct[1][0][0],
        "Gamma^X_tX": direct[1][0][1],
        "Gamma^X_XX": direct[1][1][1],
    }
    assert {
        name: sp.simplify(returned[name] - expression)
        for name, expression in expected.items()
    } == {name: 0 for name in expected}
    assert any(expression != 0 for expression in returned.values())


def test_static_connection_matches_accepted_metric_connection() -> None:
    t, x = sp.symbols("t x", real=True)
    c0 = sp.symbols("c0", positive=True)
    n = sp.Function("n", positive=True)(x)
    returned = flowing_christoffel_lab_frame(n, 0, c0, t, x)
    accepted = _independent_christoffels(optical_metric_1d(n, c0), (t, x))
    assert sp.simplify(returned["Gamma^t_tX"] - accepted[0][0][1]) == 0
    assert sp.simplify(returned["Gamma^X_tt"] - accepted[1][0][0]) == 0
    assert returned["Gamma^t_XX"] == returned["Gamma^X_tX"] == 0


def test_load_bearing_metric_mutations_break_distinct_oracles() -> None:
    n, c0 = sp.symbols("n c0", positive=True)
    velocity = sp.symbols("V", real=True)
    correct = flowing_substrate_metric(n, velocity, c0)

    missing_timelike_term = sp.Matrix(
        [[n * velocity**2, n * velocity], [n * velocity, n]]
    )
    assert missing_timelike_term.det() == 0
    assert missing_timelike_term.det() != correct.determinant

    wrong_cross_sign = sp.Matrix(
        [
            [-c0**2 / n + n * velocity**2, -n * velocity],
            [-n * velocity, n],
        ]
    )
    local_speed = sp.symbols("v_local", real=True)
    wrong_null = (
        wrong_cross_sign[0, 0]
        + 2 * wrong_cross_sign[0, 1] * local_speed
        + wrong_cross_sign[1, 1] * local_speed**2
    )
    wrong_roots = set(sp.solve(wrong_null, local_speed))
    assert wrong_roots == {velocity - c0 / n, velocity + c0 / n}
    assert wrong_roots != set(correct.massless_speeds)


@pytest.mark.parametrize(
    ("index", "velocity", "speed", "message"),
    [
        (0, 0, 1, "index"),
        (-1, 0, 1, "index"),
        (sp.Symbol("n"), 0, 1, "index"),
        (1.0, 0, 1, "index"),
        (1, sp.Symbol("V"), 1, "flow_velocity"),
        (1, 0.5, 1, "flow_velocity"),
        (1, 0, 0, "signal_speed"),
        (1, 0, -1, "signal_speed"),
        (1, 0, 1.0, "signal_speed"),
    ],
)
def test_input_guards(
    index: object,
    velocity: object,
    speed: object,
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        flowing_substrate_metric(index, velocity, speed)


def test_chart_helpers_reject_coordinate_dependent_signal_speed() -> None:
    t, x = sp.symbols("t x", real=True)
    n = sp.Function("n", positive=True)(t, x)
    velocity = sp.Function("V", real=True)(t, x)
    speed = sp.Function("c0", positive=True)(t, x)
    with pytest.raises(ValueError, match="constant"):
        flowing_christoffel_lab_frame(n, velocity, speed, t, x)
    with pytest.raises(ValueError, match="constant"):
        flowing_metric_ricci_scalar(n, velocity, speed, t, x)


def test_public_api_is_explicit_and_narrow() -> None:
    import substrate_framework.substrate_metric as module

    assert set(module.__all__) == {
        "FlowingSubstrateMetric",
        "flowing_christoffel_lab_frame",
        "flowing_metric_ricci_scalar",
        "flowing_substrate_metric",
    }
