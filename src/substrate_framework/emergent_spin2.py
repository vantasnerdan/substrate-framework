"""Rank-audited collective spin-2 completion of the sine-Gordon wall frame.

The wall construction supplies six symmetric spatial-tensor components and a
quadratic normalization.  This module keeps the logically separate steps
visible:

* local volume preservation and force balance define the admissible
  trace-free/transverse wall configurations;
* their constraint matrix has nullity two for every provably nonzero spatial
  wavevector;
* linearized relabeling invariance fixes the four two-derivative covariant
  tensor invariants to the unique Fierz--Pauli ray; and
* the wall normalization fixes the coefficient on that ray and therefore the
  source coupling.

The volume/force-balance restriction and the emergent relabeling symmetry are
declared premises of this collective model.  The functions derive all ranks,
coefficient ratios, dispersions, and coupling factors that follow from them;
they do not present those premises as consequences of a scalar ``ell=2`` mode.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Iterable

import sympy as sp

from .dimensional_sine_gordon import DimensionalSineGordonCoefficients
from .exact_symbolic import exact_real as _exact_real
from .exact_symbolic import positive_exact as _positive_exact
from .sine_gordon_wall_network import (
    SineGordonWallNetworkScales,
    sine_gordon_wall_network_scales,
)


_SPATIAL_COMPONENT_ORDER = ("xx", "yy", "zz", "xy", "xz", "yz")
_SPACETIME_COMPONENTS = (
    (0, 0),
    (0, 1),
    (0, 2),
    (0, 3),
    (1, 1),
    (1, 2),
    (1, 3),
    (2, 2),
    (2, 3),
    (3, 3),
)
_INVARIANT_ORDER = (
    "d_h_d_h",
    "div_h_div_h",
    "div_h_d_trace",
    "d_trace_d_trace",
)


@dataclass(frozen=True)
class FierzPauliGaugeLedger:
    """Exact coefficient solve for the gauge-invariant kinetic operator."""

    invariant_order: tuple[str, ...]
    sampled_constraint_rank: int
    allowed_coefficient_dimension: int
    normalized_coefficient_ray: tuple[sp.Expr, ...]
    symbolic_gauge_residual: sp.ImmutableMatrix


@dataclass(frozen=True)
class SpinTwoModeLedger:
    """Constraint, spectrum, and finite-spacing data at one wavevector."""

    component_order: tuple[str, ...]
    constraint_names: tuple[str, ...]
    wavevector: sp.ImmutableMatrix
    constraint_matrix: sp.ImmutableMatrix
    constraint_rank: int
    admissible_basis: sp.ImmutableMatrix
    physical_mode_count: int
    projected_frobenius_metric: sp.ImmutableMatrix
    projected_kinetic_rank: int
    signal_speed: sp.Expr
    continuum_angular_frequency_squared: sp.Expr
    lattice_angular_frequency_squared: sp.Expr
    leading_lattice_correction: sp.Expr
    relative_lattice_correction: sp.Expr


@dataclass(frozen=True)
class SineGordonSpinTwoCouplingLedger:
    """Gauge-completed wall normalization and its source coupling."""

    wall_scales: SineGordonWallNetworkScales
    gauge: FierzPauliGaugeLedger
    normalized_fierz_pauli_coefficients: tuple[sp.Expr, ...]
    trace_reversed_source_coefficient: sp.Expr
    einstein_coupling: sp.Expr
    newton_constant: sp.Expr


@dataclass(frozen=True)
class SineGordonSpinTwoLedger:
    """Composed coupling and physical-mode evidence at one wavevector."""

    coupling: SineGordonSpinTwoCouplingLedger
    modes: SpinTwoModeLedger


def _exact_tuple(values: Iterable[Any], length: int, name: str) -> tuple[sp.Expr, ...]:
    raw_values = tuple(values)
    if len(raw_values) != length:
        raise ValueError(
            f"{name} must contain exactly {length} entries; got {len(raw_values)}"
        )
    result: list[sp.Expr] = []
    for index, value in enumerate(raw_values):
        try:
            result.append(_exact_real(value, f"{name}[{index}]"))
        except ValueError as error:
            raise ValueError(f"{error}; got {value!r}") from error
    return tuple(result)


def _spacetime_tensor_variables() -> tuple[tuple[sp.Symbol, ...], sp.ImmutableMatrix]:
    variables = sp.symbols("h00 h01 h02 h03 h11 h12 h13 h22 h23 h33", real=True)
    entries: list[list[sp.Expr]] = [[sp.Integer(0) for _ in range(4)] for _ in range(4)]
    for variable, (mu, nu) in zip(variables, _SPACETIME_COMPONENTS, strict=True):
        entries[mu][nu] = variable
        entries[nu][mu] = variable
    return variables, sp.ImmutableMatrix(entries)


@lru_cache(maxsize=32)
def _kinetic_hessian(
    momentum: tuple[sp.Expr, ...],
    coefficients: tuple[sp.Expr, ...],
) -> sp.ImmutableMatrix:
    """Return the Hessian of the four-invariant quadratic density."""

    variables, h_covariant = _spacetime_tensor_variables()
    metric = sp.diag(-1, 1, 1, 1)
    h_contravariant = metric * h_covariant * metric
    k_covariant = sp.ImmutableMatrix(momentum)
    k_contravariant = metric * k_covariant
    momentum_squared = (k_covariant.T * k_contravariant)[0]
    trace = sp.trace(metric * h_covariant)
    tensor_norm = sum(
        h_covariant[mu, nu] * h_contravariant[mu, nu]
        for mu in range(4)
        for nu in range(4)
    )
    divergence_contravariant = sp.ImmutableMatrix(
        [
            sum(k_covariant[mu] * h_contravariant[mu, nu] for mu in range(4))
            for nu in range(4)
        ]
    )
    divergence_covariant = metric * divergence_contravariant
    invariants = (
        momentum_squared * tensor_norm,
        (divergence_contravariant.T * divergence_covariant)[0],
        sum(
            divergence_contravariant[nu] * k_covariant[nu] * trace
            for nu in range(4)
        ),
        momentum_squared * trace**2,
    )
    density = sum(
        coefficient * invariant
        for coefficient, invariant in zip(coefficients, invariants, strict=True)
    )
    return sp.ImmutableMatrix(sp.hessian(sp.expand(density), variables))


@lru_cache(maxsize=32)
def _gauge_generator(momentum: tuple[sp.Expr, ...]) -> sp.ImmutableMatrix:
    """Return ``delta h_mu_nu=k_mu*xi_nu+k_nu*xi_mu`` in component order."""

    return sp.ImmutableMatrix(
        [
            [
                (momentum[mu] if nu == alpha else 0)
                + (momentum[nu] if mu == alpha else 0)
                for alpha in range(4)
            ]
            for mu, nu in _SPACETIME_COMPONENTS
        ]
    )


def fierz_pauli_gauge_residual(
    coefficients: Iterable[Any],
    momentum: Iterable[Any],
) -> sp.ImmutableMatrix:
    """Return the exact kinetic-operator residual on a gauge variation."""

    coefficient_tuple = _exact_tuple(coefficients, 4, "coefficient")
    momentum_tuple = _exact_tuple(momentum, 4, "momentum")
    return sp.ImmutableMatrix(
        sp.simplify(
            _kinetic_hessian(momentum_tuple, coefficient_tuple)
            * _gauge_generator(momentum_tuple)
        )
    )


@lru_cache(maxsize=1)
def fierz_pauli_gauge_ledger() -> FierzPauliGaugeLedger:
    r"""Solve the most-general two-derivative quadratic action exactly.

    Gauge invariance under ``delta h_mu_nu=partial_mu xi_nu+partial_nu xi_mu``
    leaves the single ray ``(1,-2,2,-1)`` for the ordered invariants.  The
    solve stacks exact non-null and non-collinear sample momenta to determine
    the coefficient-space rank, then verifies the resulting ray at a fully
    symbolic four-momentum.
    """

    basis = tuple(
        tuple(sp.Integer(int(index == column)) for index in range(4))
        for column in range(4)
    )
    samples = (
        (sp.Integer(1), sp.Integer(2), sp.Integer(3), sp.Integer(5)),
        (sp.Integer(2), sp.Integer(-1), sp.Integer(1), sp.Integer(3)),
        (sp.Integer(3), sp.Integer(1), sp.Integer(-2), sp.Integer(4)),
    )
    constraint_rows: list[list[sp.Expr]] = []
    for momentum in samples:
        basis_residuals = [
            _kinetic_hessian(momentum, coefficient) * _gauge_generator(momentum)
            for coefficient in basis
        ]
        for row in range(10):
            for column in range(4):
                constraint_rows.append(
                    [residual[row, column] for residual in basis_residuals]
                )
    constraint_matrix = sp.ImmutableMatrix(constraint_rows)
    nullspace = constraint_matrix.nullspace()
    if len(nullspace) != 1:
        raise AssertionError("gauge solve did not select a unique coefficient ray")
    ray_vector = nullspace[0]
    ray_vector = sp.simplify(ray_vector / ray_vector[0])
    ray = tuple(ray_vector)

    k0, k1, k2, k3 = sp.symbols("k_0 k_1 k_2 k_3", real=True)
    symbolic_residual = fierz_pauli_gauge_residual(ray, (k0, k1, k2, k3))
    if symbolic_residual != sp.zeros(10, 4):
        raise AssertionError("sampled gauge ray failed the symbolic momentum check")

    return FierzPauliGaugeLedger(
        invariant_order=_INVARIANT_ORDER,
        sampled_constraint_rank=constraint_matrix.rank(),
        allowed_coefficient_dimension=len(nullspace),
        normalized_coefficient_ray=ray,
        symbolic_gauge_residual=symbolic_residual,
    )


def collective_tensor_constraint_matrix(
    wavevector: Iterable[Any],
    *,
    include_volume_constraint: bool = True,
) -> sp.ImmutableMatrix:
    """Return force-balance rows ``k_i h_ij=0`` and optional ``tr(h)=0``."""

    kx, ky, kz = _exact_tuple(wavevector, 3, "wavevector")
    rows = [
        [kx, 0, 0, ky, kz, 0],
        [0, ky, 0, kx, 0, kz],
        [0, 0, kz, 0, kx, ky],
    ]
    if include_volume_constraint:
        rows.append([1, 1, 1, 0, 0, 0])
    return sp.ImmutableMatrix(rows)


def tensor_mode_count_from_constraints(constraints: Any) -> int:
    """Derive the symmetric-tensor nullity from an exact constraint matrix."""

    matrix = sp.ImmutableMatrix(constraints)
    if matrix.cols != 6:
        raise ValueError("constraints must act on six symmetric tensor components")
    return matrix.cols - matrix.rank()


def spin_two_mode_ledger(
    wavevector: Iterable[Any],
    signal_speed: Any,
    lattice_spacing: Any,
) -> SpinTwoModeLedger:
    r"""Derive the two TT modes and their continuum/lattice dispersion.

    The lattice control is the local nearest-neighbor symbol
    ``4*c**2*sum_i sin(k_i*ell/2)**2/ell**2``.  Its leading departure from the
    common relativistic cone is
    ``-c**2*ell**2*sum_i k_i**4/12``; this exposes the preferred-frame cutoff
    rather than claiming exact microscopic Lorentz invariance.
    """

    components = _exact_tuple(wavevector, 3, "wavevector")
    if not any(component.is_zero is False for component in components):
        raise ValueError("wavevector must be provably nonzero")
    speed = _positive_exact(signal_speed, "signal_speed")
    spacing = _positive_exact(lattice_spacing, "lattice_spacing")
    constraints = collective_tensor_constraint_matrix(components)
    rank = constraints.rank()
    basis_vectors = constraints.nullspace()
    if len(basis_vectors) != tensor_mode_count_from_constraints(constraints):
        raise AssertionError("constraint nullspace dimension is inconsistent with rank")
    basis = sp.ImmutableMatrix.hstack(*basis_vectors)
    frobenius_metric = sp.diag(1, 1, 1, 2, 2, 2)
    projected_metric = sp.ImmutableMatrix(sp.simplify(basis.T * frobenius_metric * basis))
    if (
        projected_metric.rank() != len(basis_vectors)
        or projected_metric.is_positive_definite is not True
    ):
        raise AssertionError("admissible wall modes do not have positive kinetic norm")

    norm_squared = sp.simplify(sum(component**2 for component in components))
    continuum = sp.simplify(speed**2 * norm_squared)
    lattice = sp.simplify(
        4
        * speed**2
        / spacing**2
        * sum(sp.sin(component * spacing / 2) ** 2 for component in components)
    )
    leading_correction = sp.simplify(
        -speed**2 * spacing**2 * sum(component**4 for component in components) / 12
    )
    relative_correction = sp.simplify(leading_correction / continuum)

    return SpinTwoModeLedger(
        component_order=_SPATIAL_COMPONENT_ORDER,
        constraint_names=("force_x", "force_y", "force_z", "volume"),
        wavevector=sp.ImmutableMatrix(components),
        constraint_matrix=constraints,
        constraint_rank=rank,
        admissible_basis=basis,
        physical_mode_count=len(basis_vectors),
        projected_frobenius_metric=projected_metric,
        projected_kinetic_rank=projected_metric.rank(),
        signal_speed=speed,
        continuum_angular_frequency_squared=continuum,
        lattice_angular_frequency_squared=lattice,
        leading_lattice_correction=leading_correction,
        relative_lattice_correction=relative_correction,
    )


def sine_gordon_spin_two_coupling_ledger(
    coefficients: DimensionalSineGordonCoefficients,
) -> SineGordonSpinTwoCouplingLedger:
    """Complete the wall normalization on the derived gauge-invariant ray."""

    wall = sine_gordon_wall_network_scales(coefficients)
    gauge = fierz_pauli_gauge_ledger()
    normalization = sp.simplify(-wall.spin2_spacetime_normalization / 2)
    completed_coefficients = tuple(
        sp.simplify(normalization * coefficient)
        for coefficient in gauge.normalized_coefficient_ray
    )
    coupling = sp.simplify(1 / (4 * wall.spin2_spacetime_normalization))
    newton = sp.simplify(wall.signal_speed**4 * coupling / (8 * sp.pi))
    source = sp.simplify(-2 * coupling)
    return SineGordonSpinTwoCouplingLedger(
        wall_scales=wall,
        gauge=gauge,
        normalized_fierz_pauli_coefficients=completed_coefficients,
        trace_reversed_source_coefficient=source,
        einstein_coupling=coupling,
        newton_constant=newton,
    )


def sine_gordon_spin_two_ledger(
    coefficients: DimensionalSineGordonCoefficients,
    wavevector: Iterable[Any],
) -> SineGordonSpinTwoLedger:
    """Compose the gauge-completed coupling with the physical-mode audit."""

    coupling = sine_gordon_spin_two_coupling_ledger(coefficients)
    modes = spin_two_mode_ledger(
        wavevector,
        coupling.wall_scales.signal_speed,
        coupling.wall_scales.profile_length,
    )
    return SineGordonSpinTwoLedger(
        coupling=coupling,
        modes=modes,
    )
