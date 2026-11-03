"""Exact conditional algebra for a free symmetric tensor.

This module provides two deliberately separate utilities.  First, it solves
the coefficient kernel obtained after *supplying* linearized relabeling
invariance for the most-general Lorentz-covariant quadratic two-derivative
density.  Second, it computes the nullspace obtained after *supplying* spatial
trace and transversality constraints.

The utilities do not show that either structure emerges from a microscopic
model.  In particular, a positive full-rank microscopic kinetic Hessian is not
turned into a gauge theory by evaluating it only on the returned TT subspace.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Iterable

import sympy as sp

from .exact_symbolic import exact_real as _exact_real


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
    """Exact coefficient solve conditional on the supplied gauge variation."""

    invariant_order: tuple[str, ...]
    sampled_constraint_rank: int
    allowed_coefficient_dimension: int
    normalized_coefficient_ray: tuple[sp.Expr, ...]
    symbolic_gauge_residual: sp.ImmutableMatrix


@dataclass(frozen=True)
class SuppliedTTConstraintLedger:
    """Nullspace and norm after supplied trace/transversality conditions."""

    component_order: tuple[str, ...]
    constraint_names: tuple[str, ...]
    wavevector: sp.ImmutableMatrix
    constraint_matrix: sp.ImmutableMatrix
    constraint_rank: int
    admissible_basis: sp.ImmutableMatrix
    remaining_amplitude_count: int
    projected_frobenius_metric: sp.ImmutableMatrix


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
    entries: list[list[sp.Expr]] = [
        [sp.Integer(0) for _ in range(4)] for _ in range(4)
    ]
    for variable, (mu, nu) in zip(variables, _SPACETIME_COMPONENTS, strict=True):
        entries[mu][nu] = variable
        entries[nu][mu] = variable
    return variables, sp.ImmutableMatrix(entries)


@lru_cache(maxsize=32)
def _kinetic_hessian(
    momentum: tuple[sp.Expr, ...],
    coefficients: tuple[sp.Expr, ...],
) -> sp.ImmutableMatrix:
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
    """Return the kinetic residual on a supplied linearized gauge variation."""

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
    r"""Solve the coefficient ray conditional on linearized gauge invariance.

    For the ordered invariant basis, requiring
    ``delta h_mu_nu=k_mu*xi_nu+k_nu*xi_mu`` to lie in the quadratic kernel
    leaves the unique normalized ray ``(1,-2,2,-1)``.  Exact samples determine
    the coefficient rank and a symbolic four-momentum verifies the result.
    """

    coefficient_basis = tuple(
        tuple(sp.Integer(int(index == column)) for index in range(4))
        for column in range(4)
    )
    samples = (
        (sp.Integer(1), sp.Integer(2), sp.Integer(3), sp.Integer(5)),
        (sp.Integer(2), sp.Integer(-1), sp.Integer(1), sp.Integer(3)),
        (sp.Integer(3), sp.Integer(1), sp.Integer(-2), sp.Integer(4)),
    )
    rows: list[list[sp.Expr]] = []
    for momentum in samples:
        residuals = [
            _kinetic_hessian(momentum, coefficient)
            * _gauge_generator(momentum)
            for coefficient in coefficient_basis
        ]
        for row in range(10):
            for column in range(4):
                rows.append([residual[row, column] for residual in residuals])
    constraint_matrix = sp.ImmutableMatrix(rows)
    nullspace = constraint_matrix.nullspace()
    if len(nullspace) != 1:
        raise AssertionError("gauge solve did not select a unique coefficient ray")
    ray_vector = sp.simplify(nullspace[0] / nullspace[0][0])
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


def supplied_tt_constraint_matrix(
    wavevector: Iterable[Any],
    *,
    include_trace_constraint: bool = True,
) -> sp.ImmutableMatrix:
    """Return supplied rows ``k_i h_ij=0`` and optional ``tr(h)=0``."""

    kx, ky, kz = _exact_tuple(wavevector, 3, "wavevector")
    rows = [
        [kx, 0, 0, ky, kz, 0],
        [0, ky, 0, kx, 0, kz],
        [0, 0, kz, 0, kx, ky],
    ]
    if include_trace_constraint:
        rows.append([1, 1, 1, 0, 0, 0])
    return sp.ImmutableMatrix(rows)


def symmetric_tensor_mode_count(constraints: Any) -> int:
    """Return the nullity of supplied constraints on six symmetric entries."""

    matrix = sp.ImmutableMatrix(constraints)
    if matrix.cols != 6:
        raise ValueError("constraints must act on six symmetric tensor components")
    return matrix.cols - matrix.rank()


def supplied_tt_constraint_ledger(
    wavevector: Iterable[Any],
) -> SuppliedTTConstraintLedger:
    """Audit the TT nullspace at one provably nonzero spatial wavevector."""

    components = _exact_tuple(wavevector, 3, "wavevector")
    if not any(component.is_zero is False for component in components):
        raise ValueError("wavevector must be provably nonzero")
    constraints = supplied_tt_constraint_matrix(components)
    basis_vectors = constraints.nullspace()
    basis = sp.ImmutableMatrix.hstack(*basis_vectors)
    frobenius_metric = sp.diag(1, 1, 1, 2, 2, 2)
    projected_metric = sp.ImmutableMatrix(
        sp.simplify(basis.T * frobenius_metric * basis)
    )
    if (
        projected_metric.rank() != len(basis_vectors)
        or projected_metric.is_positive_definite is not True
    ):
        raise AssertionError("supplied TT subspace lost its positive norm")
    return SuppliedTTConstraintLedger(
        component_order=_SPATIAL_COMPONENT_ORDER,
        constraint_names=("transverse_x", "transverse_y", "transverse_z", "trace"),
        wavevector=sp.ImmutableMatrix(components),
        constraint_matrix=constraints,
        constraint_rank=constraints.rank(),
        admissible_basis=basis,
        remaining_amplitude_count=len(basis_vectors),
        projected_frobenius_metric=projected_metric,
    )
