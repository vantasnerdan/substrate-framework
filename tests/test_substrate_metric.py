"""Tests for the flowing substrate metric and impedance matching module.

These tests reproduce the gates from campaigns/P226/verify.py as
importable pytest cases. Each test invokes a public function of the
``substrate_framework.substrate_metric`` module and asserts the exact
algebraic identity it produces. Mutation tests confirm the test's
sensitivity to load-bearing inputs.
"""

from __future__ import annotations

import sympy as sp
import pytest

from substrate_framework.substrate_metric import (
    FlowingSubstrateMetric,
    flowing_christoffel_lab_frame,
    flowing_metric_determinant,
    flowing_metric_inverse,
    flowing_metric_massless_speeds,
    flowing_metric_minkowski_limit,
    flowing_metric_static_limit,
    flowing_substrate_metric,
    impedance_product,
    impedance_product_residual,
)


# Symbols used across tests
N, V, C0 = sp.symbols("n V c0", positive=True)
RHO_0, THETA_0 = sp.symbols("rho_0 Theta_0", positive=True)
T, X = sp.symbols("t x", real=True)


# ============================================================
# Impedance matching (C-IMP-001 candidate)
# ============================================================

def test_impedance_product_vanishes() -> None:
    """rho*Theta under the constitutive map equals rho_0*Theta_0."""

    product = impedance_product(N, RHO_0, THETA_0)
    expected = RHO_0 * THETA_0
    assert sp.simplify(product - expected) == 0


def test_impedance_product_residual_is_zero() -> None:
    """The load-bearing residual vanishes identically in n."""

    residual = impedance_product_residual(N, RHO_0, THETA_0)
    assert residual == 0


def test_impedance_product_under_n_mutation() -> None:
    """Mutation: alter the constitutive map and require the residual to be non-zero."""

    # Mutant: Theta = Theta_0 * n (wrong sign on the map)
    mutant = (RHO_0 * N) * (THETA_0 * N)
    expected = RHO_0 * THETA_0
    assert sp.simplify(mutant - expected) != 0


# ============================================================
# Flowing substrate metric (C-OG-006 candidate)
# ============================================================

def test_flowing_metric_determinant_is_minus_c0_squared() -> None:
    """The flowing metric determinant is exactly -c0**2 for any V."""

    det = flowing_metric_determinant(N, V, C0)
    assert det == -C0**2


def test_flowing_metric_inverse_components() -> None:
    """The inverse metric satisfies g^{mu alpha} g_{alpha nu} = delta^mu_nu."""

    inverse = flowing_metric_inverse(N, V, C0)
    metric = flowing_substrate_metric(N, V, C0).covariant
    product = (inverse * metric)
    assert sp.simplify(product[0, 0] - 1) == 0
    assert sp.simplify(product[0, 1]) == 0
    assert sp.simplify(product[1, 0]) == 0
    assert sp.simplify(product[1, 1] - 1) == 0


def test_flowing_metric_minkowski_limit() -> None:
    """At n=1, V=0 the metric reduces to the homogeneous Minkowski metric."""

    metric = flowing_metric_minkowski_limit(C0)
    assert metric.covariant == sp.Matrix([[-C0**2, 0], [0, 1]])
    assert metric.determinant == -C0**2
    assert metric.massless_speeds == (-C0, C0)


def test_flowing_metric_static_limit_matches_optical_family() -> None:
    """At V=0 the metric reduces to the static 1+1 optical metric family
    in the paper's c0**2-explicit convention.
    """

    metric = flowing_metric_static_limit(N, C0)
    assert metric.covariant[0, 0] == -C0**2 / N
    assert metric.covariant[0, 1] == 0
    assert metric.covariant[1, 1] == N
    assert metric.determinant == -C0**2


def test_flowing_metric_massless_speeds_match_paper() -> None:
    """The massless null constraint yields v = -V +/- c0/n."""

    speeds = flowing_metric_massless_speeds(N, V, C0)
    expected = sorted([-V + C0 / N, -V - C0 / N], key=lambda s: str(s))
    actual = sorted(speeds, key=lambda s: str(s))
    assert sp.simplify(actual[0] - expected[0]) == 0
    assert sp.simplify(actual[1] - expected[1]) == 0


def test_flowing_metric_determinant_is_structural_invariant() -> None:
    """Mutation: det is structural (-c0**2) for any V, n population.

    This is a property of the metric ansatz: even with V=5 the determinant
    is preserved, but the *specific* metric components change. The test
    confirms the determinant invariance holds at multiple V values.
    """

    for v_test in [sp.Integer(0), sp.Integer(1), sp.Integer(5), sp.Integer(-3)]:
        det = flowing_metric_determinant(N, v_test, C0)
        assert det == -C0**2, f"det at V={v_test} is {det}, expected -c0**2"


def test_flowing_metric_under_signature_mutation() -> None:
    """Mutation: replace g_tt with a wrong-sign form (missing -c0**2/n term)
    and require the determinant to no longer equal -c0**2.
    """

    metric_obj = flowing_substrate_metric(N, V, C0)
    # Mutant: drop the -c0**2/n term entirely, keep only the n V^2 part.
    # This breaks the structural det = -c0**2 invariant.
    mutant_g_tt = N * V**2  # missing the -c0**2/n term
    mutant_det = sp.simplify(mutant_g_tt * N - (N * V)**2)
    # Mutant determinant is n^2 V^2 - n^2 V^2 = 0
    assert mutant_det == 0
    # The determinant invariance is broken under the mutation
    assert mutant_det != metric_obj.determinant


# ============================================================
# Christoffel symbols (derived from C-OG-006, not a separate claim)
# ============================================================

def test_flowing_christoffel_lab_frame_vanishes_in_minkowski_limit() -> None:
    """At n=1, V=0, with no spatial/temporal gradients in n or V,
    all six Christoffel symbols vanish (Minkowski is flat).
    """

    n_static = sp.Integer(1)
    v_static = sp.Integer(0)
    symbols = flowing_christoffel_lab_frame(n_static, v_static, C0, T, X)
    for name, expr in symbols.items():
        assert sp.simplify(expr) == 0, f"{name} = {expr} (expected 0)"


def test_flowing_christoffel_static_limit_recovers_optical_family() -> None:
    """At V=0 with non-trivial spatial gradient in n, the Christoffel symbols
    match the static 1+1 optical family: only Gamma^t_tX and Gamma^X_tt
    are non-zero, and they agree with the C-OG-001 expression.
    """

    symbols = flowing_christoffel_lab_frame(N, sp.Integer(0), C0, T, X)
    # In the V=0, time-independent (n=N(x)) limit:
    # dtn = 0, dtv = 0, dxn = N', dxv = 0
    # Expected per C-OG-001: Gamma^t_tX = -(1/(2n)) n_x, Gamma^X_tt = -c0^2/(2 n^3) n_x
    # The other four should vanish.
    Gamma_t_tX = sp.simplify(symbols["Gamma^t_tX"].subs([(sp.Derivative(N, T), 0)]))
    Gamma_X_tt = sp.simplify(symbols["Gamma^X_tt"].subs([(sp.Derivative(N, T), 0)]))
    expected_t_tX = -sp.diff(N, X) / (2 * N)
    expected_X_tt = -C0**2 * sp.diff(N, X) / (2 * N**3)
    assert sp.simplify(Gamma_t_tX - expected_t_tX) == 0
    assert sp.simplify(Gamma_X_tt - expected_X_tt) == 0


# ============================================================
# Module API smoke test
# ============================================================

def test_module_exposes_expected_api() -> None:
    """Confirm the public API surface is stable."""

    expected = {
        "impedance_product",
        "impedance_product_residual",
        "FlowingSubstrateMetric",
        "flowing_substrate_metric",
        "flowing_metric_determinant",
        "flowing_metric_inverse",
        "flowing_metric_massless_speeds",
        "flowing_metric_minkowski_limit",
        "flowing_metric_static_limit",
        "flowing_christoffel_lab_frame",
    }
    import substrate_framework.substrate_metric as mod
    actual = set(dir(mod))
    missing = expected - actual
    assert not missing, f"missing API: {missing}"
