"""Sine-Gordon plaquette normalization and teleparallel geometry ledgers.

The microscopic input is a compact plaquette phase ``q`` on a network with
one sine-Gordon channel per transverse area ``ell**2``.  The channel uses the
same positive onsite line density ``mu`` as the declared dimensional
sine-Gordon action.  On the principal compact branch, the chord map

``tau = 2*sin(q/2)/ell``

turns its lifted cosine density into ``-mu*tau**2/2`` exactly.  Local Lorentz
invariance selects the unique three-invariant TEGR quadratic form up to an
overall multiplier; the channel density fixes that remaining multiplier.

With the mostly-plus torsion convention used below,
``R(Levi-Civita) + T - B = 0``.  The collective action is therefore
``-mu*T/2 = mu*R/2`` up to the displayed boundary density.  Consequently
``kappa=1/mu`` and ``G=c**4/(8*pi*mu)``.  No observed Newton constant, bare
Einstein term, determinant regulator, or fitted coefficient enters this
module.

The local-Lorentz uniqueness input is the coframe constitutive theorem in
Itin, Hehl, and Obukhov, arXiv:1611.05759.  The code independently constructs
both sides of the nonlinear teleparallel identity for supplied coframes.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Sequence

import sympy as sp

from .dimensional_sine_gordon import (
    DimensionalSineGordonCoefficients,
    DimensionalSineGordonScales,
    dimensional_sine_gordon_coefficients,
    dimensional_sine_gordon_scales,
)
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
class CompactTorsionChannelLedger:
    """Exact compact-phase to collective-torsion normalization for one channel."""

    phase: sp.Expr
    cell_length: sp.Expr
    onsite_line_density: sp.Expr
    transverse_line_density: sp.Expr
    constitutive_weight: sp.Expr
    chord_torsion: sp.Expr
    geometric_small_phase_torsion: sp.Expr
    microscopic_lagrangian_density: sp.Expr
    quadratic_collective_density: sp.Expr
    cosine_chord_identity_residual: sp.Expr


@dataclass(frozen=True)
class SineGordonLinkCoframe:
    """Collective coframe and metric of four microscopic cell-link covectors."""

    microscopic_link_covectors: sp.ImmutableMatrix
    cell_length: sp.Expr
    collective_coframe: sp.ImmutableMatrix
    metric_covariant: sp.ImmutableMatrix
    metric_contravariant: sp.ImmutableMatrix
    volume_density: sp.Expr


@dataclass(frozen=True)
class CompactTeleparallelActionLedger:
    """All 24 compact plaquette channels and their exact TEGR quadratic map."""

    local_channel_order: tuple[tuple[int, int, int], ...]
    phases: sp.ImmutableMatrix
    cell_length: sp.Expr
    constitutive_matrix: sp.ImmutableMatrix
    spectral_basis: sp.ImmutableMatrix
    spectral_weights: sp.ImmutableMatrix
    spectral_chord_torsion_channels: sp.ImmutableMatrix
    local_chord_torsion_channels: sp.ImmutableMatrix
    microscopic_lagrangian_density: sp.Expr
    collective_teleparallel_density: sp.Expr
    compact_collective_identity_residual: sp.Expr


@dataclass(frozen=True)
class SineGordonTeleparallelGravity:
    """Upstream scale and coefficient closure of the SG-to-TEGR construction."""

    coefficients: DimensionalSineGordonCoefficients
    scales: DimensionalSineGordonScales
    cell_length: sp.Expr
    transverse_line_density: sp.Expr
    torsion_invariant_weights: tuple[sp.Expr, sp.Expr, sp.Expr]
    quadratic_coframe_parameter_count: int
    local_lorentz_ratio_dimension: int
    teleparallel_action_coefficient: sp.Expr
    einstein_hilbert_coefficient: sp.Expr
    einstein_stress_coupling: sp.Expr
    newton_constant: sp.Expr
    radiative_tensor_polarizations: int
    independent_gravitational_parameters: int


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


def _validated_coefficients(
    coefficients: DimensionalSineGordonCoefficients,
) -> DimensionalSineGordonCoefficients:
    if not isinstance(coefficients, DimensionalSineGordonCoefficients):
        raise TypeError("coefficients must be a DimensionalSineGordonCoefficients record")
    return dimensional_sine_gordon_coefficients(
        coefficients.inertia,
        coefficients.gradient,
        coefficients.onsite,
    )


def _quadratic_weights(values: Sequence[Any]) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    if len(values) != 3:
        raise ValueError("torsion invariant weights must contain exactly three values")
    weights = tuple(_exact_real(value, "torsion invariant weight") for value in values)
    if all(value == 0 for value in weights):
        raise ValueError("torsion invariant weights cannot all vanish")
    return weights  # type: ignore[return-value]


def compact_torsion_channel_ledger(
    phase: Any,
    coefficients: DimensionalSineGordonCoefficients,
    *,
    constitutive_weight: Any = 1,
) -> CompactTorsionChannelLedger:
    r"""Return the exact lifted-cosine/chord-torsion identity for one channel.

    The physical cell length is not a new parameter: it is the sine-Gordon
    length ``ell=sqrt(T/mu)``.  A single line channel per area ``ell**2``
    lifts ``-mu*(1-cos(q))`` to a volume density.  Defining the signed chord
    torsion on the principal branch by ``tau=2*sin(q/2)/ell`` gives

    ``-mu*(1-cos(q))/ell**2 = -mu*tau**2/2``.

    ``constitutive_weight`` is a component of the Lorentz-selected quadratic
    form, not an additional overall coupling.
    """

    q = _exact_real(phase, "phase")
    weight = _exact_real(constitutive_weight, "constitutive_weight")
    validated = _validated_coefficients(coefficients)
    scales = dimensional_sine_gordon_scales(validated)
    length = scales.length
    transverse_density = sp.simplify(1 / length**2)
    chord_torsion = sp.simplify(2 * sp.sin(q / 2) / length)
    microscopic = sp.simplify(
        -validated.onsite
        * transverse_density
        * weight
        * (1 - sp.cos(q))
    )
    quadratic = sp.simplify(-validated.onsite * weight * chord_torsion**2 / 2)
    residual = sp.trigsimp(microscopic - quadratic)
    if residual != 0:
        raise AssertionError("compact cosine did not reduce to chord torsion exactly")
    return CompactTorsionChannelLedger(
        phase=q,
        cell_length=length,
        onsite_line_density=validated.onsite,
        transverse_line_density=transverse_density,
        constitutive_weight=weight,
        chord_torsion=chord_torsion,
        geometric_small_phase_torsion=sp.simplify(q / length),
        microscopic_lagrangian_density=microscopic,
        quadratic_collective_density=quadratic,
        cosine_chord_identity_residual=residual,
    )


def collective_sine_gordon_link_coframe(
    microscopic_link_covectors: Any,
    coefficients: DimensionalSineGordonCoefficients,
) -> SineGordonLinkCoframe:
    r"""Map four physical cell-link covectors to a Lorentzian geometry.

    Rows of ``L**a_mu`` are one clock-like and three spatial microscopic link
    covectors with length dimension.  Dividing by the derived SG cell length
    gives ``e**a_mu=L**a_mu/ell``.  The collective geometry is then the typed
    contraction ``g=e.T*diag(-1,1,1,1)*e``; invertibility of the link matrix
    is exactly the nondegeneracy condition.  No refractive-index analogy or
    independently supplied metric enters the map.
    """

    links = sp.Matrix(microscopic_link_covectors)
    if links.shape != (4, 4):
        raise ValueError("microscopic_link_covectors must have shape 4 by 4")
    if any(entry.has(sp.Float) for entry in links):
        raise ValueError("microscopic_link_covectors must be exact rather than floating")
    if sp.simplify(links.det()).is_zero is True:
        raise ValueError("microscopic_link_covectors must be invertible")
    validated = _validated_coefficients(coefficients)
    length = dimensional_sine_gordon_scales(validated).length
    coframe = (links / length).applyfunc(sp.simplify)
    metric = (coframe.T * sp.diag(-1, 1, 1, 1) * coframe).applyfunc(sp.simplify)
    inverse_metric = metric.inv().applyfunc(sp.simplify)
    return SineGordonLinkCoframe(
        microscopic_link_covectors=sp.ImmutableMatrix(links),
        cell_length=length,
        collective_coframe=sp.ImmutableMatrix(coframe),
        metric_covariant=sp.ImmutableMatrix(metric),
        metric_contravariant=sp.ImmutableMatrix(inverse_metric),
        volume_density=sp.simplify(sp.sqrt(-metric.det())),
    )


def _local_torsion_quadratic(channels: Sequence[Any]) -> sp.Expr:
    if len(channels) != len(TORSION_CHANNELS):
        raise ValueError("local torsion channels must contain exactly 24 values")
    signs = (-1, 1, 1, 1)
    torsion = [[[sp.Integer(0) for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for value, (internal, first, second) in zip(channels, TORSION_CHANNELS):
        expression = sp.sympify(value)
        torsion[internal][first][second] = expression
        torsion[internal][second][first] = -expression
    torsion_down = [
        [
            [signs[internal] * torsion[internal][first][second] for second in range(4)]
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
    r"""Return the explicit 24-channel matrix ``K`` with ``T=tau.T*K*tau``.

    Channel order is ``(a,b,c)`` with ``a=0..3`` and ``0<=b<c<=3`` for the
    local coframe components ``T**a_bc``.  Building ``K`` from the three
    invariant contractions makes the compact-phase-to-coframe map executable;
    no phrase such as "torsion-like" is used as an identification.
    """

    symbols = sp.symbols(f"tau0:{len(TORSION_CHANNELS)}", real=True)
    quadratic = _local_torsion_quadratic(symbols)
    matrix = sp.hessian(quadratic, symbols) / 2
    if matrix != matrix.T:
        raise AssertionError("teleparallel constitutive matrix must be symmetric")
    reconstructed = sp.expand((sp.Matrix(symbols).T * matrix * sp.Matrix(symbols))[0])
    if sp.expand(reconstructed - quadratic) != 0:
        raise AssertionError("constitutive matrix did not reconstruct the torsion scalar")
    return sp.ImmutableMatrix(matrix)


@lru_cache(maxsize=1)
def teleparallel_constitutive_spectral_basis() -> tuple[sp.ImmutableMatrix, sp.ImmutableMatrix]:
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


def compact_teleparallel_action_ledger(
    phases: Sequence[Any],
    coefficients: DimensionalSineGordonCoefficients,
) -> CompactTeleparallelActionLedger:
    r"""Map all compact coframe plaquettes to the TEGR density exactly.

    ``q_I`` are 24 compact phases in an exact orthonormal eigenbasis of the
    Lorentz-selected constitutive matrix.  Their chord torsions are
    ``s_I=2*sin(q_I/2)/ell`` and their fixed spectral weights are ``d_I``.
    The microscopic density is the ordinary signed cosine sum

    ``-mu/ell**2 * sum_I d_I*(1-cos(q_I))``.

    It is identically ``-mu*s.T*D*s/2``.  The typed orthogonal map
    ``tau=U*s`` then gives ``-mu*tau.T*K*tau/2=-mu*T/2`` in local coframe
    channels.  Each microscopic phase remains individually ``2*pi`` periodic.
    """

    if len(phases) != len(TORSION_CHANNELS):
        raise ValueError("phases must contain exactly 24 compact plaquette values")
    phase_vector = sp.ImmutableMatrix(
        [_exact_real(phase, "compact plaquette phase") for phase in phases]
    )
    validated = _validated_coefficients(coefficients)
    scales = dimensional_sine_gordon_scales(validated)
    constitutive = teleparallel_constitutive_matrix_mostly_plus()
    spectral_basis, spectral_weights = teleparallel_constitutive_spectral_basis()
    spectral_torsion = sp.ImmutableMatrix(
        2 * phase_vector.applyfunc(lambda phase: sp.sin(phase / 2)) / scales.length
    )
    local_torsion = sp.ImmutableMatrix(spectral_basis * spectral_torsion)
    microscopic = sp.trigsimp(
        -validated.onsite
        / scales.length**2
        * sum(
            spectral_weights[index, index] * (1 - sp.cos(phase_vector[index]))
            for index in range(len(TORSION_CHANNELS))
        )
    )
    collective = sp.trigsimp(
        -validated.onsite
        / 2
        * (local_torsion.T * constitutive * local_torsion)[0]
    )
    residual = sp.trigsimp(microscopic - collective)
    if residual != 0:
        raise AssertionError("compact plaquette action did not reduce to TEGR")
    return CompactTeleparallelActionLedger(
        local_channel_order=TORSION_CHANNELS,
        phases=phase_vector,
        cell_length=scales.length,
        constitutive_matrix=constitutive,
        spectral_basis=spectral_basis,
        spectral_weights=spectral_weights,
        spectral_chord_torsion_channels=spectral_torsion,
        local_chord_torsion_channels=local_torsion,
        microscopic_lagrangian_density=microscopic,
        collective_teleparallel_density=collective,
        compact_collective_identity_residual=residual,
    )


def sine_gordon_teleparallel_gravity(
    coefficients: DimensionalSineGordonCoefficients,
) -> SineGordonTeleparallelGravity:
    r"""Derive the nonlinear Einstein coefficient from SG cell parameters.

    The general parity-even local quadratic torsion action has three invariant
    coefficients.  The imported local-Lorentz theorem leaves one common
    multiplier and fixes the conventional TEGR ratios ``(1/4, 1/2, -1)``.
    The exact channel normalization fixes the collective action coefficient
    to ``-mu/2``.  Since ``R=-T+B``, the Einstein-Hilbert coefficient is
    ``mu/2=1/(2*kappa)``.
    """

    validated = _validated_coefficients(coefficients)
    scales = dimensional_sine_gordon_scales(validated)
    onsite = validated.onsite
    teleparallel_coefficient = sp.simplify(-onsite / 2)
    einstein_coefficient = sp.simplify(onsite / 2)
    einstein_coupling = sp.simplify(1 / onsite)
    newton = sp.simplify(scales.signal_speed**4 / (8 * sp.pi * onsite))
    if sp.simplify(einstein_coefficient - 1 / (2 * einstein_coupling)) != 0:
        raise AssertionError("Einstein-Hilbert and stress-coupling coefficients disagree")
    if sp.simplify(
        einstein_coupling - 8 * sp.pi * newton / scales.signal_speed**4
    ) != 0:
        raise AssertionError("Newton and Einstein stress couplings disagree")
    return SineGordonTeleparallelGravity(
        coefficients=validated,
        scales=scales,
        cell_length=scales.length,
        transverse_line_density=sp.simplify(1 / scales.length**2),
        torsion_invariant_weights=TEGR_INVARIANT_WEIGHTS,
        quadratic_coframe_parameter_count=3,
        local_lorentz_ratio_dimension=1,
        teleparallel_action_coefficient=teleparallel_coefficient,
        einstein_hilbert_coefficient=einstein_coefficient,
        einstein_stress_coupling=einstein_coupling,
        newton_constant=newton,
        radiative_tensor_polarizations=2,
        independent_gravitational_parameters=0,
    )


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
                    sum(metric[rho, sigma] * torsion[sigma][mu][nu] for sigma in range(4))
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
    ``v_m=T**n_nm``, and
    ``T=w1*I1+w2*I2+w3*v_m*v**m``.

    At the TEGR weights ``(1/4,1/2,-1)``, the independently constructed
    Levi-Civita scalar obeys ``R+T-B=0`` with
    ``B=2/e*partial_m(e*v**m)``.  Non-TEGR weights are accepted deliberately
    so mutation tests can observe a nonzero residual.
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
    constitutive_residual = sp.simplify(torsion_scalar - constitutive_quadratic)
    if weights == TEGR_INVARIANT_WEIGHTS and constitutive_residual != 0:
        raise AssertionError(
            "local 24-channel constitutive map did not reproduce the torsion scalar"
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
