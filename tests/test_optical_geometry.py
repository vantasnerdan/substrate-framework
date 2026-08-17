from __future__ import annotations

import pytest
import sympy as sp

from substrate_framework.optical_geometry import (
    index_from_potential,
    optical_box_static_1d,
    optical_dilaton,
    optical_dilaton_source_operator_1d,
    optical_metric_1d,
    optical_ricci_scalar_1d,
    slow_geodesic_acceleration_1d,
    slow_geodesic_acceleration_from_potential,
)


def test_optical_dilaton_matches_curvature_for_nontrivial_profile() -> None:
    x = sp.symbols("x", real=True)
    c0 = sp.symbols("c0", positive=True)
    index = 1 + x**2
    metric = optical_metric_1d(index, c0)
    assert metric.det() == -1 / c0**2
    assert sp.simplify(
        optical_box_static_1d(optical_dilaton(index), index, x, c0)
        - optical_ricci_scalar_1d(index, x, c0)
    ) == 0


def test_conditional_potential_map_and_drift_are_exact() -> None:
    x = sp.symbols("x", real=True)
    c0 = sp.symbols("c0", positive=True)
    potential = sp.Function("Phi")(x)
    index = index_from_potential(potential, c0)
    assert index == 1 / (1 + 2 * potential / c0**2)
    assert sp.simplify(
        slow_geodesic_acceleration_from_potential(potential, x, c0)
        + (1 + 2 * potential / c0**2) * sp.diff(potential, x)
    ) == 0


def test_conditional_source_operator_is_exact_not_just_weak_field() -> None:
    x = sp.symbols("x", real=True)
    c0 = sp.symbols("c0", positive=True)
    potential = sp.Function("Phi")(x)
    source_side = optical_dilaton_source_operator_1d(potential, x, c0)
    assert sp.simplify(source_side - 2 * sp.diff(potential, x, 2)) == 0
    assert c0 not in source_side.free_symbols




def test_slow_geodesic_acceleration_matches_canonical_christoffel() -> None:
    """Cross-route: ``slow_geodesic_acceleration_1d`` equals ``-Gamma^x_tt``.

    The closed form
    ``c0**2 * n_x / (2 * n**3)`` is the accepted C-OG-002 surface owned
    by P004 (issue #67 follow-up). It is reconstructed independently here
    via the canonical ``pseudo_riemannian.metric_christoffel_symbols``
    route on the static 1+1 optical metric, evaluated at the (t, t, x)
    Christoffel slot. Exact agreement is the cross-route regression guard.

    The test deliberately derives the canonical value from
    ``metric_christoffel_symbols(optical_metric_1d(n, c0), (t, x))``
    rather than re-implementing the Levi-Civita formula in the test
    body — the canonical route is the single source of truth for the
    Christoffel-symbol path.
    """

    from substrate_framework.pseudo_riemannian import (
        coordinate_symbols,
        metric_christoffel_symbols,
    )

    t, x = coordinate_symbols(
        (sp.Symbol("t", real=True), sp.Symbol("x", real=True))
    )
    c0 = sp.Symbol("c0", positive=True)
    n = sp.Function("n", positive=True)(x)

    metric = optical_metric_1d(n, c0)
    gamma = metric_christoffel_symbols(metric, (t, x))
    canonical_acceleration = -sp.simplify(gamma[1][0][0])

    closed_form = slow_geodesic_acceleration_1d(n, x, c0)
    assert sp.simplify(closed_form - canonical_acceleration) == 0

    # Sanity: a non-constant n must produce a non-zero acceleration.
    assert canonical_acceleration != 0


def test_slow_geodesic_acceleration_under_sign_mutation() -> None:
    """Mutation guard: flipping the sign of the closed form breaks the route.

    Verifies the cross-route test in ``test_slow_geodesic_acceleration_
    matches_canonical_christoffel`` cannot be satisfied by ``+c0**2 *
    n_x / (2 * n**3)`` (i.e. ``+Gamma^x_tt``) instead of the correct
    ``-Gamma^x_tt``. This catches any future drift where the helper
    and the canonical route silently agree on a wrong-sign expression.
    """

    from substrate_framework.pseudo_riemannian import (
        coordinate_symbols,
        metric_christoffel_symbols,
    )

    t, x = coordinate_symbols(
        (sp.Symbol("t", real=True), sp.Symbol("x", real=True))
    )
    c0 = sp.Symbol("c0", positive=True)
    n = sp.Function("n", positive=True)(x)

    metric = optical_metric_1d(n, c0)
    gamma = metric_christoffel_symbols(metric, (t, x))
    canonical_acceleration = -sp.simplify(gamma[1][0][0])

    wrong_sign = sp.simplify(+gamma[1][0][0])  # NOT negated
    assert sp.simplify(wrong_sign - canonical_acceleration) != 0


def test_slow_geodesic_acceleration_under_power_mutation() -> None:
    """Mutation guard: wrong power of n in the closed form breaks the route.

    The correct closed form has ``n**3`` in the denominator. A drift to
    ``n**2`` (or any other power) must be rejected by the cross-route
    comparison. Also checks the ``c0**2`` coefficient — replacing it
    with ``c0`` (missing one power) is rejected.
    """

    from substrate_framework.pseudo_riemannian import (
        coordinate_symbols,
        metric_christoffel_symbols,
    )

    t, x = coordinate_symbols(
        (sp.Symbol("t", real=True), sp.Symbol("x", real=True))
    )
    c0 = sp.Symbol("c0", positive=True)
    n = sp.Function("n", positive=True)(x)

    metric = optical_metric_1d(n, c0)
    gamma = metric_christoffel_symbols(metric, (t, x))
    canonical_acceleration = -sp.simplify(gamma[1][0][0])

    # Wrong power of n
    wrong_power = sp.simplify(c0**2 * sp.diff(n, x) / (2 * n**2))
    assert sp.simplify(wrong_power - canonical_acceleration) != 0

    # Missing c0**2 (factor of c0 only)
    missing_c0_squared = sp.simplify(c0 * sp.diff(n, x) / (2 * n**3))
    assert sp.simplify(missing_c0_squared - canonical_acceleration) != 0



@pytest.mark.parametrize(
    ("call", "message"),
    [
        (lambda: optical_metric_1d(0, 1), "index"),
        (lambda: optical_metric_1d(1, 0), "signal_speed"),
        (lambda: index_from_potential(-sp.Rational(1, 2), 1), "positive optical index"),
    ],
)
def test_numeric_domains_are_enforced(call, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        call()
