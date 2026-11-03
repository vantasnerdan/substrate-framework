"""Exact local-U(1) gauge algebra for the accepted complex-scalar convention.

This module declares a connection and covariant derivative but no Maxwell
kinetic term, gauge-field equation, propagating photon, or physical charge map.
"""

from __future__ import annotations

from typing import Any

import sympy as sp


def _positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    return expression


def u1_covariant_derivative(
    field: Any,
    connection_component: Any,
    coordinate: sp.Symbol,
    coupling: Any,
) -> sp.Expr:
    """Return ``(partial_mu-i*e*A_mu) field`` for charge convention ``+e``."""

    e = _positive(coupling, "coupling")
    psi = sp.sympify(field)
    return sp.simplify(
        sp.diff(psi, coordinate)
        - sp.I * e * sp.sympify(connection_component) * psi
    )


def local_u1_transform(
    field: Any,
    connections: tuple[Any, Any],
    gauge_parameter: Any,
    coordinates: tuple[sp.Symbol, sp.Symbol],
    coupling: Any,
) -> tuple[sp.Expr, tuple[sp.Expr, sp.Expr]]:
    """Transform a charge-``+e`` field and two covariant connection components."""

    e = _positive(coupling, "coupling")
    chi = sp.sympify(gauge_parameter)
    phase = sp.exp(sp.I * e * chi)
    transformed_connections = tuple(
        sp.sympify(connection) + sp.diff(chi, coordinate)
        for connection, coordinate in zip(connections, coordinates, strict=True)
    )
    return phase * sp.sympify(field), transformed_connections


def gauged_scalar_kinetic_density(
    field: Any,
    conjugate_field: Any,
    connections: tuple[Any, Any],
    coordinates: tuple[sp.Symbol, sp.Symbol],
    coupling: Any,
) -> sp.Expr:
    """Return ``(D_0 Psi)* D_0 Psi-(D_1 Psi)* D_1 Psi`` in signature ``(+,-)``."""

    e = _positive(coupling, "coupling")
    psi = sp.sympify(field)
    psi_conjugate = sp.sympify(conjugate_field)
    derivatives: list[tuple[sp.Expr, sp.Expr]] = []
    for connection, coordinate in zip(connections, coordinates, strict=True):
        gauge = sp.sympify(connection)
        derivative = sp.diff(psi, coordinate) - sp.I * e * gauge * psi
        conjugate_derivative = (
            sp.diff(psi_conjugate, coordinate)
            + sp.I * e * gauge * psi_conjugate
        )
        derivatives.append((derivative, conjugate_derivative))
    return sp.simplify(
        derivatives[0][1] * derivatives[0][0]
        - derivatives[1][1] * derivatives[1][0]
    )


def u1_field_strength(
    connection_mu: Any,
    connection_nu: Any,
    coordinate_mu: sp.Symbol,
    coordinate_nu: sp.Symbol,
) -> sp.Expr:
    """Return ``F_mu_nu=partial_mu A_nu-partial_nu A_mu``."""

    return sp.simplify(
        sp.diff(sp.sympify(connection_nu), coordinate_mu)
        - sp.diff(sp.sympify(connection_mu), coordinate_nu)
    )


def u1_covariant_commutator(
    field: Any,
    connection_mu: Any,
    connection_nu: Any,
    coordinate_mu: sp.Symbol,
    coordinate_nu: sp.Symbol,
    coupling: Any,
) -> sp.Expr:
    """Return ``[D_mu,D_nu] field`` derived from the declared connection."""

    e = _positive(coupling, "coupling")
    psi = sp.sympify(field)
    a_mu, a_nu = sp.sympify(connection_mu), sp.sympify(connection_nu)
    d_mu_d_nu = u1_covariant_derivative(
        u1_covariant_derivative(psi, a_nu, coordinate_nu, e),
        a_mu,
        coordinate_mu,
        e,
    )
    d_nu_d_mu = u1_covariant_derivative(
        u1_covariant_derivative(psi, a_mu, coordinate_mu, e),
        a_nu,
        coordinate_nu,
        e,
    )
    return sp.simplify(d_mu_d_nu - d_nu_d_mu)


def finite_energy_winding_flux(winding: Any, coupling: Any) -> sp.Expr:
    """Return ``2*pi*N/e`` for an explicitly integer finite-energy winding."""

    integer = sp.sympify(winding)
    if integer.is_integer is not True:
        raise ValueError("winding must be an integer")
    e = _positive(coupling, "coupling")
    return sp.simplify(2 * sp.pi * integer / e)


def u1_holonomy(flux: Any, coupling: Any) -> sp.Expr:
    """Return the charge-``e`` loop phase ``exp(i*e*flux)``."""

    e = _positive(coupling, "coupling")
    return sp.simplify(sp.exp(sp.I * e * sp.sympify(flux)))
