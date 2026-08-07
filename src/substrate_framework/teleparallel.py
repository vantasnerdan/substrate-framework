"""Exact mostly-plus teleparallel geometry and constitutive ledgers.

This module evaluates the standard TEGR torsion scalar and the identity
``R + T - B = 0`` for a supplied nondegenerate coframe in Weitzenbock gauge.
It also exposes the exact 24-component constitutive quadratic form.  It does
not derive a coframe, a microscopic action, Newton's constant, matter
coupling, or a compact lattice realization.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Sequence

import sympy as sp

from .exact_symbolic import exact_real as _exact_real


TEGR_INVARIANT_WEIGHTS = (
    sp.Rational(1, 4),
    sp.Rational(1, 2),
    sp.Integer(-1),
)
TORSION_CHANNELS = tuple(
    (internal, first, second)
    for internal in range(4)
    for first in range(4)
    for second in range(first + 1, 4)
)


@dataclass(frozen=True)
class TeleparallelCoframeLedger:
    """Exact tensors and the TEGR/Einstein identity for one coframe."""

    coordinates: tuple[sp.Symbol, sp.Symbol, sp.Symbol, sp.Symbol]
    coframe: sp.ImmutableMatrix
    inverse_coframe: sp.ImmutableMatrix
    metric_covariant: sp.ImmutableMatrix
    metric_contravariant: sp.ImmutableMatrix
    volume_density: sp.Expr
    torsion_invariant_weights: tuple[sp.Expr, sp.Expr, sp.Expr]
    torsion_invariant_one: sp.Expr
    torsion_invariant_two: sp.Expr
    torsion_vector_norm_squared: sp.Expr
    torsion_scalar: sp.Expr
    torsion_vector_covariant: sp.ImmutableMatrix
    torsion_vector_contravariant: sp.ImmutableMatrix
    local_torsion_channels: sp.ImmutableMatrix
    constitutive_quadratic_residual: sp.Expr
    boundary_divergence: sp.Expr
    levi_civita_ricci_scalar: sp.Expr
    einstein_teleparallel_identity_residual: sp.Expr


def _quadratic_weights(values: Sequence[Any]) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    if len(values) != 3:
        raise ValueError("torsion invariant weights must contain exactly three values")
    weights = tuple(_exact_real(value, "torsion invariant weight") for value in values)
    if all(value == 0 for value in weights):
        raise ValueError("torsion invariant weights cannot all vanish")
    return weights  # type: ignore[return-value]


def _local_torsion_quadratic(channels: Sequence[Any]) -> sp.Expr:
    if len(channels) != len(TORSION_CHANNELS):
        raise ValueError("local torsion channels must contain exactly 24 values")
    signs = (-1, 1, 1, 1)
    torsion = [
        [[sp.Integer(0) for _ in range(4)] for _ in range(4)]
        for _ in range(4)
    ]
    for value, (internal, first, second) in zip(channels, TORSION_CHANNELS):
        expression = sp.sympify(value)
        torsion[internal][first][second] = expression
        torsion[internal][second][first] = -expression
    torsion_down = [
        [
            [
                signs[internal] * torsion[internal][first][second]
                for second in range(4)
            ]
            for first in range(4)
        ]
        for internal in range(4)
    ]
    torsion_up = [
        [
            [
                signs[internal]
                * signs[first]
                * signs[second]
                * torsion_down[internal][first][second]
                for second in range(4)
            ]
            for first in range(4)
        ]
        for internal in range(4)
    ]
    invariant_one = sum(
        torsion_up[internal][first][second]
        * torsion_down[internal][first][second]
        for internal in range(4)
        for first in range(4)
        for second in range(4)
    )
    invariant_two = sum(
        torsion_up[internal][first][second]
        * torsion_down[second][first][internal]
        for internal in range(4)
        for first in range(4)
        for second in range(4)
    )
    torsion_vector = [
        sum(torsion[internal][internal][first] for internal in range(4))
        for first in range(4)
    ]
    vector_norm = sum(
        signs[first] * torsion_vector[first] ** 2 for first in range(4)
    )
    return sp.expand(invariant_one / 4 + invariant_two / 2 - vector_norm)


@lru_cache(maxsize=1)
def teleparallel_constitutive_matrix_mostly_plus() -> sp.ImmutableMatrix:
    r"""Return ``K`` with ``T=tau.T*K*tau`` in the canonical 24 channels.

    Channel order is ``(a,b,c)`` with ``a=0..3`` and ``0<=b<c<=3`` for the
    local coframe components ``T**a_bc``.
    """

    symbols = sp.symbols(f"tau0:{len(TORSION_CHANNELS)}", real=True)
    quadratic = _local_torsion_quadratic(symbols)
    matrix = sp.hessian(quadratic, symbols) / 2
    if matrix != matrix.T:
        raise AssertionError("teleparallel constitutive matrix must be symmetric")
    reconstructed = sp.expand(
        (sp.Matrix(symbols).T * matrix * sp.Matrix(symbols))[0]
    )
    if sp.expand(reconstructed - quadratic) != 0:
        raise AssertionError("constitutive matrix did not reconstruct torsion scalar")
    return sp.ImmutableMatrix(matrix)


@lru_cache(maxsize=1)
def teleparallel_constitutive_spectral_basis(
) -> tuple[sp.ImmutableMatrix, sp.ImmutableMatrix]:
    """Return an exact orthonormal eigenbasis and diagonal weights for ``K``."""

    constitutive = sp.Matrix(teleparallel_constitutive_matrix_mostly_plus())
    columns: list[sp.Matrix] = []
    eigenvalues: list[sp.Expr] = []
    for eigenvalue, multiplicity, vectors in sorted(
        constitutive.eigenvects(),
        key=lambda item: float(item[0]),
    ):
        orthonormal = sp.GramSchmidt(vectors, orthonormal=True)
        if len(orthonormal) != multiplicity:
            raise AssertionError("constitutive eigenspace lost multiplicity")
        columns.extend(orthonormal)
        eigenvalues.extend([eigenvalue] * multiplicity)
    basis = sp.Matrix.hstack(*columns)
    diagonal = sp.diag(*eigenvalues)
    if basis.T * basis != sp.eye(len(TORSION_CHANNELS)):
        raise AssertionError("constitutive spectral basis must be orthonormal")
    if sp.simplify(basis.T * constitutive * basis - diagonal) != sp.zeros(
        len(TORSION_CHANNELS)
    ):
        raise AssertionError("constitutive spectral basis did not diagonalize K")
    return sp.ImmutableMatrix(basis), sp.ImmutableMatrix(diagonal)


def _validated_coframe(
    coframe: Any,
    coordinates: Sequence[sp.Symbol],
) -> tuple[sp.Matrix, tuple[sp.Symbol, sp.Symbol, sp.Symbol, sp.Symbol]]:
    frame = sp.Matrix(coframe)
    if frame.shape != (4, 4):
        raise ValueError("coframe must have shape 4 by 4")
    if len(coordinates) != 4 or any(
        not isinstance(coordinate, sp.Symbol) for coordinate in coordinates
    ):
        raise ValueError("coordinates must contain exactly four SymPy Symbols")
    if any(entry.has(sp.Float) for entry in frame):
        raise ValueError("coframe must be exact rather than floating")
    determinant = sp.simplify(frame.det())
    if determinant.is_zero is True:
        raise ValueError("coframe must be invertible")
    return frame, tuple(coordinates)  # type: ignore[return-value]


def _coordinate_torsion(
    frame: sp.Matrix,
    inverse_frame: sp.Matrix,
    coordinates: tuple[sp.Symbol, sp.Symbol, sp.Symbol, sp.Symbol],
) -> list[list[list[sp.Expr]]]:
    """Return ``T**rho_mu_nu`` from the coframe curl."""

    return [
        [
            [
                sp.simplify(
                    sum(
                        inverse_frame[rho, internal]
                        * (
                            sp.diff(frame[internal, nu], coordinates[mu])
                            - sp.diff(frame[internal, mu], coordinates[nu])
                        )
                        for internal in range(4)
                    )
                )
                for nu in range(4)
            ]
            for mu in range(4)
        ]
        for rho in range(4)
    ]


def _torsion_index_forms(
    torsion: list[list[list[sp.Expr]]],
    metric: sp.Matrix,
    inverse_metric: sp.Matrix,
) -> tuple[list[list[list[sp.Expr]]], list[list[list[sp.Expr]]]]:
    """Return the all-lowered and all-raised torsion components."""

    lowered = [
        [
            [
                sp.simplify(
                    sum(
                        metric[rho, sigma] * torsion[sigma][mu][nu]
                        for sigma in range(4)
                    )
                )
                for nu in range(4)
            ]
            for mu in range(4)
        ]
        for rho in range(4)
    ]
    raised = [
        [
            [
                sp.simplify(
                    sum(
                        inverse_metric[rho, alpha]
                        * inverse_metric[mu, beta]
                        * inverse_metric[nu, gamma]
                        * lowered[alpha][beta][gamma]
                        for alpha in range(4)
                        for beta in range(4)
                        for gamma in range(4)
                    )
                )
                for nu in range(4)
            ]
            for mu in range(4)
        ]
        for rho in range(4)
    ]
    return lowered, raised


def _torsion_contractions(
    torsion: list[list[list[sp.Expr]]],
    lowered: list[list[list[sp.Expr]]],
    raised: list[list[list[sp.Expr]]],
    inverse_metric: sp.Matrix,
) -> tuple[sp.Expr, sp.Expr, list[sp.Expr], list[sp.Expr], sp.Expr]:
    """Return the two tensor contractions, vector, and vector norm."""

    invariant_one = sp.simplify(
        sum(
            raised[rho][mu][nu] * lowered[rho][mu][nu]
            for rho in range(4)
            for mu in range(4)
            for nu in range(4)
        )
    )
    invariant_two = sp.simplify(
        sum(
            raised[rho][mu][nu] * lowered[nu][mu][rho]
            for rho in range(4)
            for mu in range(4)
            for nu in range(4)
        )
    )
    vector_down = [
        sp.simplify(sum(torsion[nu][nu][mu] for nu in range(4)))
        for mu in range(4)
    ]
    vector_up = [
        sp.simplify(
            sum(inverse_metric[mu, nu] * vector_down[nu] for nu in range(4))
        )
        for mu in range(4)
    ]
    vector_norm = sp.simplify(
        sum(vector_down[mu] * vector_up[mu] for mu in range(4))
    )
    return invariant_one, invariant_two, vector_down, vector_up, vector_norm


def _local_torsion_channel_vector(
    frame: sp.Matrix,
    inverse_frame: sp.Matrix,
    coordinates: tuple[sp.Symbol, sp.Symbol, sp.Symbol, sp.Symbol],
) -> sp.ImmutableMatrix:
    """Return the 24 local ``T**a_bc`` components in canonical order."""

    return sp.ImmutableMatrix(
        [
            sp.simplify(
                sum(
                    inverse_frame[mu, first]
                    * inverse_frame[nu, second]
                    * (
                        sp.diff(frame[internal, nu], coordinates[mu])
                        - sp.diff(frame[internal, mu], coordinates[nu])
                    )
                    for mu in range(4)
                    for nu in range(4)
                )
            )
            for internal, first, second in TORSION_CHANNELS
        ]
    )


def _levi_civita_ricci_scalar(
    metric: sp.Matrix,
    inverse_metric: sp.Matrix,
    coordinates: tuple[sp.Symbol, sp.Symbol, sp.Symbol, sp.Symbol],
) -> sp.Expr:
    """Construct the Levi-Civita scalar independently of torsion."""

    christoffel = [
        [
            [
                sp.simplify(
                    sum(
                        inverse_metric[rho, sigma]
                        * (
                            sp.diff(metric[sigma, nu], coordinates[mu])
                            + sp.diff(metric[sigma, mu], coordinates[nu])
                            - sp.diff(metric[mu, nu], coordinates[sigma])
                        )
                        / 2
                        for sigma in range(4)
                    )
                )
                for nu in range(4)
            ]
            for mu in range(4)
        ]
        for rho in range(4)
    ]
    ricci = [
        [
            sp.simplify(
                sum(
                    sp.diff(christoffel[rho][mu][nu], coordinates[rho])
                    - sp.diff(christoffel[rho][mu][rho], coordinates[nu])
                    + sum(
                        christoffel[rho][rho][sigma]
                        * christoffel[sigma][mu][nu]
                        - christoffel[rho][nu][sigma]
                        * christoffel[sigma][mu][rho]
                        for sigma in range(4)
                    )
                    for rho in range(4)
                )
            )
            for nu in range(4)
        ]
        for mu in range(4)
    ]
    return sp.simplify(
        sum(
            inverse_metric[mu, nu] * ricci[mu][nu]
            for mu in range(4)
            for nu in range(4)
        )
    )


def teleparallel_coframe_ledger(
    coframe: Any,
    coordinates: Sequence[sp.Symbol],
    *,
    torsion_invariant_weights: Sequence[Any] = TEGR_INVARIANT_WEIGHTS,
) -> TeleparallelCoframeLedger:
    r"""Construct ``R``, ``T``, and ``B`` independently from one coframe.

    For ``e**a_mu`` and ``eta_ab=diag(-1,1,1,1)``, this function uses
    ``T**rho_mn=e_a**rho*(partial_m e**a_n-partial_n e**a_m)``,
    ``v_m=T**n_nm``, and ``T=w1*I1+w2*I2+w3*v_m*v**m``.  At the TEGR
    weights, the independently constructed Levi-Civita scalar obeys
    ``R+T-B=0`` with ``B=2/e*partial_m(e*v**m)``.

    Non-TEGR weights are accepted deliberately so mutation tests can observe
    a nonzero residual.  The function uses a zero inertial spin connection;
    callers are responsible for supplying the corresponding coframe gauge.
    """

    frame, coordinate_tuple = _validated_coframe(coframe, coordinates)
    weights = _quadratic_weights(torsion_invariant_weights)
    eta = sp.diag(-1, 1, 1, 1)
    inverse_frame = frame.inv().applyfunc(sp.simplify)
    metric = (frame.T * eta * frame).applyfunc(sp.simplify)
    inverse_metric = metric.inv().applyfunc(sp.simplify)
    volume = sp.simplify(sp.sqrt(-metric.det()))

    torsion = _coordinate_torsion(frame, inverse_frame, coordinate_tuple)
    torsion_down, torsion_up = _torsion_index_forms(
        torsion,
        metric,
        inverse_metric,
    )
    (
        invariant_one,
        invariant_two,
        torsion_vector_down,
        torsion_vector_up,
        vector_norm,
    ) = _torsion_contractions(
        torsion,
        torsion_down,
        torsion_up,
        inverse_metric,
    )
    torsion_scalar = sp.simplify(
        weights[0] * invariant_one
        + weights[1] * invariant_two
        + weights[2] * vector_norm
    )
    local_torsion_channels = _local_torsion_channel_vector(
        frame,
        inverse_frame,
        coordinate_tuple,
    )
    constitutive_quadratic = sp.simplify(
        (
            local_torsion_channels.T
            * teleparallel_constitutive_matrix_mostly_plus()
            * local_torsion_channels
        )[0]
    )
    constitutive_residual = sp.simplify(
        torsion_scalar - constitutive_quadratic
    )
    if weights == TEGR_INVARIANT_WEIGHTS and constitutive_residual != 0:
        raise AssertionError(
            "local 24-channel constitutive map did not reproduce torsion scalar"
        )
    boundary = sp.simplify(
        2
        / volume
        * sum(
            sp.diff(volume * torsion_vector_up[mu], coordinate_tuple[mu])
            for mu in range(4)
        )
    )
    ricci_scalar = _levi_civita_ricci_scalar(
        metric,
        inverse_metric,
        coordinate_tuple,
    )
    identity_residual = sp.simplify(ricci_scalar + torsion_scalar - boundary)
    if weights == TEGR_INVARIANT_WEIGHTS and identity_residual != 0:
        raise AssertionError(
            "mostly-plus teleparallel identity failed: "
            f"R+T-B={identity_residual!s}"
        )
    return TeleparallelCoframeLedger(
        coordinates=coordinate_tuple,
        coframe=sp.ImmutableMatrix(frame),
        inverse_coframe=sp.ImmutableMatrix(inverse_frame),
        metric_covariant=sp.ImmutableMatrix(metric),
        metric_contravariant=sp.ImmutableMatrix(inverse_metric),
        volume_density=volume,
        torsion_invariant_weights=weights,
        torsion_invariant_one=invariant_one,
        torsion_invariant_two=invariant_two,
        torsion_vector_norm_squared=vector_norm,
        torsion_scalar=torsion_scalar,
        torsion_vector_covariant=sp.ImmutableMatrix(torsion_vector_down),
        torsion_vector_contravariant=sp.ImmutableMatrix(torsion_vector_up),
        local_torsion_channels=local_torsion_channels,
        constitutive_quadratic_residual=constitutive_residual,
        boundary_divergence=boundary,
        levi_civita_ricci_scalar=ricci_scalar,
        einstein_teleparallel_identity_residual=identity_residual,
    )
