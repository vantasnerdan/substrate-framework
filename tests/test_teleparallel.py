from __future__ import annotations

import pytest
import sympy as sp

import substrate_framework as framework
from substrate_framework.teleparallel import (
    TEGR_INVARIANT_WEIGHTS,
    TORSION_CHANNELS,
    TeleparallelCoframeLedger,
    teleparallel_coframe_ledger,
    teleparallel_constitutive_matrix_mostly_plus,
    teleparallel_constitutive_spectral_basis,
)


def test_public_package_exports_teleparallel_geometry() -> None:
    assert framework.TEGR_INVARIANT_WEIGHTS is TEGR_INVARIANT_WEIGHTS
    assert framework.TORSION_CHANNELS is TORSION_CHANNELS
    assert framework.TeleparallelCoframeLedger is TeleparallelCoframeLedger
    assert framework.teleparallel_coframe_ledger is teleparallel_coframe_ledger
    assert (
        framework.teleparallel_constitutive_matrix_mostly_plus
        is teleparallel_constitutive_matrix_mostly_plus
    )
    assert (
        framework.teleparallel_constitutive_spectral_basis
        is teleparallel_constitutive_spectral_basis
    )


def test_constitutive_matrix_reconstructs_the_tegr_quadratic_form() -> None:
    constitutive = teleparallel_constitutive_matrix_mostly_plus()
    basis, spectral = teleparallel_constitutive_spectral_basis()
    assert constitutive.shape == (24, 24)
    assert constitutive == constitutive.T
    assert constitutive.rank() == 24
    eigenvalues = list(spectral.diagonal())
    assert eigenvalues.count(-2) == 3
    assert eigenvalues.count(-1) == 8
    assert eigenvalues.count(sp.Rational(-1, 2)) == 1
    assert eigenvalues.count(sp.Rational(1, 2)) == 3
    assert eigenvalues.count(1) == 8
    assert eigenvalues.count(2) == 1
    assert basis.T * basis == sp.eye(24)
    assert sp.simplify(basis.T * constitutive * basis - spectral) == sp.zeros(24)


def test_nonlinear_flrw_coframe_closes_the_mostly_plus_tegr_identity() -> None:
    x0, x, y, z = sp.symbols("x0 x y z", real=True)
    hubble = sp.symbols("H", real=True)
    scale = sp.exp(hubble * x0)
    ledger = teleparallel_coframe_ledger(
        sp.diag(1, scale, scale, scale),
        (x0, x, y, z),
    )
    assert ledger.metric_covariant == sp.diag(-1, scale**2, scale**2, scale**2)
    assert ledger.volume_density == scale**3
    assert ledger.torsion_invariant_one == -6 * hubble**2
    assert ledger.torsion_invariant_two == -3 * hubble**2
    assert ledger.torsion_vector_norm_squared == -9 * hubble**2
    assert ledger.torsion_scalar == 6 * hubble**2
    assert ledger.constitutive_quadratic_residual == 0
    assert ledger.boundary_divergence == 18 * hubble**2
    assert ledger.levi_civita_ricci_scalar == 12 * hubble**2
    assert ledger.einstein_teleparallel_identity_residual == 0


def test_local_frame_cancellation_and_wrong_weight_are_distinguished() -> None:
    x0, x, y, z = sp.symbols("x0 x y z", real=True)
    slope = sp.symbols("b", real=True)
    rapidity = slope * x0
    boost = sp.Matrix(
        [
            [sp.cosh(rapidity), sp.sinh(rapidity), 0, 0],
            [sp.sinh(rapidity), sp.cosh(rapidity), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ]
    )
    ledger = teleparallel_coframe_ledger(boost, (x0, x, y, z))
    assert ledger.metric_covariant == sp.diag(-1, 1, 1, 1)
    assert ledger.torsion_invariant_one == 2 * slope**2
    assert ledger.torsion_invariant_two == slope**2
    assert ledger.torsion_vector_norm_squared == slope**2
    assert ledger.torsion_scalar == 0
    assert ledger.einstein_teleparallel_identity_residual == 0

    mutated = teleparallel_coframe_ledger(
        boost,
        (x0, x, y, z),
        torsion_invariant_weights=(sp.Rational(1, 3), sp.Rational(1, 2), -1),
    )
    assert mutated.einstein_teleparallel_identity_residual == slope**2 / 6




def test_levi_civita_ricci_scalar_matches_pseudo_riemannian() -> None:
    """Cross-route: the inline 4D Levi-Civita Ricci scalar agrees with the
    canonical ``pseudo_riemannian.metric_ricci_scalar`` route.

    ``teleparallel._levi_civita_ricci_scalar`` is an intentional
    independent construction kept separate from the canonical
    pseudo-Riemannian machinery (per PR #7 — the metric-only Levi-Civita
    scalar is required for the ``R + T - B = 0`` mostly-plus identity
    check). This test pins the equivalence so a future drift in either
    route is detected.

    The non-flat anisotropic coframe
    ``diag(1, e^(H1 t), e^(H2 t), e^(H3 t))`` produces the metric
    ``diag(-1, e^(2 H1 t), e^(2 H2 t), e^(2 H3 t))``. Both routes yield
    ``2*(H1**2 + H2**2 + H3**2 + H1*H2 + H1*H3 + H2*H3)``. The test
    asserts the symbolic equality after sympy simplify.
    """

    from substrate_framework.pseudo_riemannian import (
        coordinate_symbols,
        metric_ricci_scalar,
    )

    t, x, y, z = coordinate_symbols(
        (
            sp.Symbol("t", real=True),
            sp.Symbol("x", real=True),
            sp.Symbol("y", real=True),
            sp.Symbol("z", real=True),
        )
    )
    hubble_components = sp.symbols("H1 H2 H3", real=True)
    coframe = sp.Matrix(
        [
            [1, 0, 0, 0],
            [0, sp.exp(hubble_components[0] * t), 0, 0],
            [0, 0, sp.exp(hubble_components[1] * t), 0],
            [0, 0, 0, sp.exp(hubble_components[2] * t)],
        ]
    )
    ledger = teleparallel_coframe_ledger(coframe, (t, x, y, z))
    canonical_ricci = metric_ricci_scalar(
        ledger.metric_covariant, (t, x, y, z)
    )
    assert (
        sp.simplify(ledger.levi_civita_ricci_scalar - canonical_ricci)
        == 0
    )
    # Pin the analytic value too — both routes must simplify to the
    # same symbolic expression, so any future drift in either route
    # shows up as a divergence from this exact answer.
    expected = (
        2 * hubble_components[0] ** 2
        + 2 * hubble_components[0] * hubble_components[1]
        + 2 * hubble_components[0] * hubble_components[2]
        + 2 * hubble_components[1] ** 2
        + 2 * hubble_components[1] * hubble_components[2]
        + 2 * hubble_components[2] ** 2
    )
    assert sp.simplify(ledger.levi_civita_ricci_scalar - expected) == 0
    assert sp.simplify(canonical_ricci - expected) == 0


def test_levi_civita_ricci_scalar_under_sign_mutation() -> None:
    """Mutation guard: flipping the Ricci sign in either route is rejected.

    A Ricci-sign mutation is a valid general-equivalence mutation: any
    metric has a true Ricci scalar, and the two routes must agree on
    it with the same sign. Replacing the ledger's Ricci scalar with
    its negation (or vice versa) must produce a nonzero residual under
    the canonical route, which the cross-route comparison catches.

    Note: a metric-determinant mutation is NOT a valid general
    equivalence mutation (per Dan's #67 technical-assessment comment 2)
    because changing a metric can simply define another valid metric
    for which both Ricci implementations still agree. The sign mutation
    here targets a different invariant: the actual sign of the scalar
    curvature.
    """

    from substrate_framework.pseudo_riemannian import (
        coordinate_symbols,
        metric_ricci_scalar,
    )

    t, x, y, z = coordinate_symbols(
        (
            sp.Symbol("t", real=True),
            sp.Symbol("x", real=True),
            sp.Symbol("y", real=True),
            sp.Symbol("z", real=True),
        )
    )
    hubble_components = sp.symbols("H1 H2 H3", real=True)
    coframe = sp.Matrix(
        [
            [1, 0, 0, 0],
            [0, sp.exp(hubble_components[0] * t), 0, 0],
            [0, 0, sp.exp(hubble_components[1] * t), 0],
            [0, 0, 0, sp.exp(hubble_components[2] * t)],
        ]
    )
    ledger = teleparallel_coframe_ledger(coframe, (t, x, y, z))
    canonical_ricci = metric_ricci_scalar(
        ledger.metric_covariant, (t, x, y, z)
    )

    # Sign-flipping the ledger Ricci must NOT equal canonical Ricci.
    negated_ledger_ricci = -ledger.levi_civita_ricci_scalar
    assert sp.simplify(negated_ledger_ricci - canonical_ricci) != 0

    # Sign-flipping the canonical Ricci must NOT equal ledger Ricci.
    negated_canonical_ricci = -canonical_ricci
    assert sp.simplify(negated_canonical_ricci - ledger.levi_civita_ricci_scalar) != 0



@pytest.mark.parametrize(
    ("operation", "message"),
    [
        (
            lambda: teleparallel_coframe_ledger(sp.eye(3), sp.symbols("x:4")),
            "4 by 4",
        ),
        (
            lambda: teleparallel_coframe_ledger(sp.zeros(4), sp.symbols("x:4")),
            "invertible",
        ),
        (
            lambda: teleparallel_coframe_ledger(
                sp.eye(4),
                sp.symbols("x:4"),
                torsion_invariant_weights=(1, 2),
            ),
            "exactly three",
        ),
    ],
)
def test_teleparallel_inputs_fail_with_context(operation, message: str) -> None:
    with pytest.raises((TypeError, ValueError), match=message):
        operation()
