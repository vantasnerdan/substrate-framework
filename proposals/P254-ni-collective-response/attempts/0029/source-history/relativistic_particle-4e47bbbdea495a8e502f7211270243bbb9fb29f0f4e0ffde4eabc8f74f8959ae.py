"""Exact conditional algebra for a relativistic point-particle worldline.

Authority status: implements accepted claims C-WLN-001 through C-WLN-003 in
release v0.160.0.
The module uses the mostly-plus convention and keeps the positive-mass
auxiliary-field elimination separate from the zero-mass specialization. It
does not derive a background metric, quantize a particle, identify a physical
species, or imply field-theory conformal symmetry or helicity content.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .exact_symbolic import exact_real as _exact_real
from .exact_symbolic import positive_exact as _positive_exact
from .pseudo_riemannian import (
    exact_metric_matrix,
    metric_christoffel_from_derivatives,
    metric_inverse,
)


def _nonnegative_exact(value: Any, name: str) -> sp.Expr:
    expression = _exact_real(value, name)
    if expression.is_nonnegative is not True:
        raise ValueError(f"{name} must be provably nonnegative")
    return expression


def _exact_column(values: Any, name: str) -> sp.ImmutableMatrix:
    matrix = sp.ImmutableMatrix(values)
    if matrix.cols != 1:
        raise ValueError(f"{name} must be a column vector")
    entries = tuple(matrix)
    for index, value in enumerate(entries):
        _exact_real(value, f"{name}[{index}]")
    return matrix


@dataclass(frozen=True)
class MassiveEinbeinLedger:
    """Positive-mass auxiliary-field elimination data."""

    velocity_norm: sp.Expr
    einbein: sp.Expr
    mass: sp.Expr
    signal_speed: sp.Expr
    lagrangian: sp.Expr
    einbein_euler_derivative: sp.Expr
    mass_shell_constraint: sp.Expr
    positive_einbein_root: sp.Expr
    recovered_square_root_lagrangian: sp.Expr


@dataclass(frozen=True)
class MasslessEinbeinLedger:
    """Zero-mass action and auxiliary-field constraint data."""

    velocity_norm: sp.Expr
    einbein: sp.Expr
    lagrangian: sp.Expr
    einbein_euler_derivative: sp.Expr
    null_constraint: sp.Expr


@dataclass(frozen=True)
class EinbeinHamiltonianLedger:
    """Covector momentum, inverse Legendre map, and pure constraint."""

    metric_covariant: sp.ImmutableMatrix
    metric_contravariant: sp.ImmutableMatrix
    covariant_momentum: sp.ImmutableMatrix
    velocity: sp.ImmutableMatrix
    lagrangian: sp.Expr
    legendre_transform: sp.Expr
    hamiltonian: sp.Expr
    mass_shell_constraint: sp.Expr


def einbein_lagrangian(
    velocity_norm: Any,
    einbein: Any,
    mass: Any,
    signal_speed: Any,
) -> sp.Expr:
    r"""Return ``sigma/(2e) - e (m c0)^2/2`` for exact supplied scalars."""

    sigma = _exact_real(velocity_norm, "velocity_norm")
    e_value = _positive_exact(einbein, "einbein")
    mass_value = _nonnegative_exact(mass, "mass")
    speed = _positive_exact(signal_speed, "signal_speed")
    return sp.simplify(
        sigma / (2 * e_value) - e_value * (mass_value * speed) ** 2 / 2
    )


def massive_einbein_ledger(
    velocity_norm: Any,
    einbein: Any,
    mass: Any,
    signal_speed: Any,
) -> MassiveEinbeinLedger:
    """Derive the positive auxiliary-field root and square-root action."""

    sigma = _exact_real(velocity_norm, "velocity_norm")
    if sigma.is_negative is not True:
        raise ValueError("velocity_norm must be provably negative")
    e_value = _positive_exact(einbein, "einbein")
    mass_value = _positive_exact(mass, "mass")
    speed = _positive_exact(signal_speed, "signal_speed")
    lagrangian = einbein_lagrangian(sigma, e_value, mass_value, speed)
    euler = sp.simplify(
        -sigma / (2 * e_value**2) - (mass_value * speed) ** 2 / 2
    )
    constraint = sp.simplify(sigma + e_value**2 * (mass_value * speed) ** 2)
    root = sp.simplify(sp.sqrt(-sigma) / (mass_value * speed))
    recovered = einbein_lagrangian(sigma, root, mass_value, speed)
    expected = sp.simplify(-mass_value * speed * sp.sqrt(-sigma))
    if sp.simplify(recovered - expected) != 0:
        raise AssertionError("positive einbein elimination did not recover the root action")
    return MassiveEinbeinLedger(
        velocity_norm=sigma,
        einbein=e_value,
        mass=mass_value,
        signal_speed=speed,
        lagrangian=lagrangian,
        einbein_euler_derivative=euler,
        mass_shell_constraint=constraint,
        positive_einbein_root=root,
        recovered_square_root_lagrangian=recovered,
    )


def massless_einbein_ledger(
    velocity_norm: Any,
    einbein: Any,
) -> MasslessEinbeinLedger:
    """Return the zero-mass action whose auxiliary equation is ``sigma=0``."""

    sigma = _exact_real(velocity_norm, "velocity_norm")
    e_value = _positive_exact(einbein, "einbein")
    lagrangian = sp.simplify(sigma / (2 * e_value))
    return MasslessEinbeinLedger(
        velocity_norm=sigma,
        einbein=e_value,
        lagrangian=lagrangian,
        einbein_euler_derivative=sp.simplify(-sigma / (2 * e_value**2)),
        null_constraint=sigma,
    )


def einbein_hamiltonian_ledger(
    metric: Any,
    covariant_momentum: Any,
    einbein: Any,
    mass: Any,
    signal_speed: Any,
) -> EinbeinHamiltonianLedger:
    """Derive the pure-constraint Hamiltonian from canonical covector data."""

    metric_matrix = exact_metric_matrix(metric)
    inverse = metric_inverse(metric_matrix)
    momentum = _exact_column(covariant_momentum, "covariant_momentum")
    if momentum.rows != metric_matrix.rows:
        raise ValueError("covariant_momentum dimension must match metric")
    e_value = _positive_exact(einbein, "einbein")
    mass_value = _nonnegative_exact(mass, "mass")
    speed = _positive_exact(signal_speed, "signal_speed")
    velocity = sp.ImmutableMatrix(e_value * inverse * momentum)
    velocity_norm = sp.simplify((velocity.T * metric_matrix * velocity)[0])
    lagrangian = einbein_lagrangian(
        velocity_norm, e_value, mass_value, speed
    )
    legendre = sp.simplify((momentum.T * velocity)[0] - lagrangian)
    constraint = sp.simplify(
        (momentum.T * inverse * momentum)[0] + (mass_value * speed) ** 2
    )
    hamiltonian = sp.simplify(e_value * constraint / 2)
    if sp.simplify(legendre - hamiltonian) != 0:
        raise AssertionError("Legendre transform did not close on the constraint")
    return EinbeinHamiltonianLedger(
        metric_covariant=sp.ImmutableMatrix(metric_matrix),
        metric_contravariant=sp.ImmutableMatrix(inverse),
        covariant_momentum=momentum,
        velocity=velocity,
        lagrangian=lagrangian,
        legendre_transform=legendre,
        hamiltonian=hamiltonian,
        mass_shell_constraint=constraint,
    )


def worldline_reparametrization_residual(
    velocity_norm: Any,
    einbein: Any,
    mass: Any,
    signal_speed: Any,
    parameter_rate: Any,
) -> sp.Expr:
    """Return the transformed density-times-Jacobian minus the original."""

    sigma = _exact_real(velocity_norm, "velocity_norm")
    e_value = _positive_exact(einbein, "einbein")
    mass_value = _nonnegative_exact(mass, "mass")
    speed = _positive_exact(signal_speed, "signal_speed")
    rate = _positive_exact(parameter_rate, "parameter_rate")
    original = einbein_lagrangian(sigma, e_value, mass_value, speed)
    transformed = einbein_lagrangian(
        sigma / rate**2,
        e_value / rate,
        mass_value,
        speed,
    )
    return sp.simplify(rate * transformed - original)


def constant_e_gauge_rate(einbein: Any, target_einbein: Any) -> sp.Expr:
    """Return ``df/dtau=e/e_target`` that makes ``e/(df/dtau)`` constant."""

    e_value = _positive_exact(einbein, "einbein")
    target = _positive_exact(target_einbein, "target_einbein")
    return sp.simplify(e_value / target)


def massless_worldline_weyl_residual(
    velocity_norm: Any,
    einbein: Any,
    conformal_factor: Any,
) -> sp.Expr:
    """Return the massless density change under ``g,e -> Omega^2 g,e``."""

    sigma = _exact_real(velocity_norm, "velocity_norm")
    e_value = _positive_exact(einbein, "einbein")
    factor = _positive_exact(conformal_factor, "conformal_factor")
    original = massless_einbein_ledger(sigma, e_value).lagrangian
    transformed = massless_einbein_ledger(
        factor**2 * sigma,
        factor**2 * e_value,
    ).lagrangian
    return sp.simplify(transformed - original)


def massive_mass_term_weyl_change(
    einbein: Any,
    mass: Any,
    signal_speed: Any,
    conformal_factor: Any,
) -> sp.Expr:
    """Return the massive-term change under ``e -> Omega^2 e``."""

    e_value = _positive_exact(einbein, "einbein")
    mass_value = _positive_exact(mass, "mass")
    speed = _positive_exact(signal_speed, "signal_speed")
    factor = _positive_exact(conformal_factor, "conformal_factor")
    original = -e_value * (mass_value * speed) ** 2 / 2
    transformed = -factor**2 * e_value * (mass_value * speed) ** 2 / 2
    return sp.simplify(transformed - original)


def einbein_geodesic_acceleration(
    metric: Any,
    metric_derivatives: Any,
    velocity: Any,
    einbein: Any,
    einbein_derivative: Any,
) -> sp.ImmutableMatrix:
    r"""Return ``-Gamma^a_bc v^b v^c + (edot/e) v^a`` from local data."""

    metric_matrix = exact_metric_matrix(metric)
    vector = _exact_column(velocity, "velocity")
    if vector.rows != metric_matrix.rows:
        raise ValueError("velocity dimension must match metric")
    e_value = _positive_exact(einbein, "einbein")
    e_dot = _exact_real(einbein_derivative, "einbein_derivative")
    gamma = metric_christoffel_from_derivatives(
        metric_inverse(metric_matrix), metric_derivatives
    )
    dimension = metric_matrix.rows
    return sp.ImmutableMatrix(
        [
            sp.simplify(
                -sum(
                    gamma[upper, rho, nu] * vector[rho] * vector[nu]
                    for rho in range(dimension)
                    for nu in range(dimension)
                )
                + e_dot * vector[upper] / e_value
            )
            for upper in range(dimension)
        ]
    )
