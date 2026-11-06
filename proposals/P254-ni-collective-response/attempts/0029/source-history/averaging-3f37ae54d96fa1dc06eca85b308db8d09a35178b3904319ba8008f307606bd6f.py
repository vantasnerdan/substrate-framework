"""Exact coarse-graining machinery for constant-density Euler balance.

Conventions:

- a spatial filter is convolution with a normalized, symmetric kernel
  ``g`` of width ``Delta``: ``(S q)(x) = int g_Delta(y) q(x+y) dy``;
- symmetric normalized kernels have vanishing odd moments and even
  moments ``m_2k = int eta**(2k) g(eta) d eta``;
- on polynomials the moment series terminates, so every identity in this
  module is exact on a declared function class and linear in its input;
- the sub-filter momentum flux is ``Pi_ij = S(v_i v_j) - ubar_i ubar_j``
  with ``ubar = S v``; no closure is introduced here.

Importing this module executes no simulation.
"""

from __future__ import annotations

from typing import Any

import sympy as sp


def _top_hat_moment(order: int) -> sp.Rational:
    """Even moment of the unit top-hat shape on [-1/2, 1/2]."""

    return sp.Rational(1, (2 * order + 1) * 2 ** (2 * order))


def kernel_even_moments(kernel: str, max_order: int) -> list[sp.Expr]:
    """Even moments m_0, m_2, ..., m_{2*max_order} for a declared kernel."""

    if kernel not in {"tophat", "gaussian"}:
        raise ValueError("kernel must be 'tophat' or 'gaussian'")
    moments: list[sp.Expr] = []
    for k in range(max_order + 1):
        if kernel == "tophat":
            moments.append(_top_hat_moment(k))
        else:
            moment = sp.Rational(1)
            for j in range(1, k + 1):
                moment *= sp.Integer(2 * j - 1)
            moments.append(moment)
    return moments


def filter_polynomial(
    polynomial: sp.Expr,
    variable: sp.Symbol,
    width: Any,
    kernel: str = "tophat",
) -> sp.Expr:
    """Exact S[q] for polynomial q via the terminating moment series.

    For a polynomial the Taylor series of q(x+y) terminates, so
    ``S[q] = sum_k m_{2k} Delta^{2k} q^{(2k)} / (2k)!`` is exact, not an
    approximation. Direct integration for the top-hat kernel agrees and is
    checked in tests.
    """

    width = sp.sympify(width)
    poly = sp.Poly(sp.expand(polynomial), variable)
    degree = poly.degree()
    if degree < 0:  # sympy yields -oo for constant polynomials
        degree = 0
    half = max(int((degree + 1) // 2), 0)
    moments = kernel_even_moments(kernel, half)
    filtered = sp.Integer(0)
    for k in range(half + 1):
        derivative = sp.diff(polynomial, variable, 2 * k)
        if derivative != 0:
            filtered += moments[k] * width ** (2 * k) / sp.factorial(2 * k) * derivative
    return sp.expand(filtered)


def filter_direct_tophat(
    polynomial: sp.Expr, variable: sp.Symbol, width: Any
) -> sp.Expr:
    """Exact top-hat convolution over the cell [x - w/2, x + w/2].

    The support half-width is ``width/2`` so the shape lives on
    [-1/2, 1/2] and its second moment is 1/12, matching the series
    convention exactly.
    """

    width = sp.sympify(width)
    y = sp.Symbol("_filter_y", real=True)
    integrand = polynomial.subs(variable, variable + y) / width
    return sp.expand(sp.integrate(integrand, (y, -width / 2, width / 2)))


def commutation_residual(
    polynomial: sp.Expr,
    variable: sp.Symbol,
    width: Any,
    kernel: str = "tophat",
) -> sp.Expr:
    """Residual of [S, d/dx] q; identically zero for convolution filters."""

    filtered = filter_polynomial(polynomial, variable, width, kernel)
    filtered_derivative = filter_polynomial(
        sp.diff(polynomial, variable), variable, width, kernel
    )
    return sp.simplify(sp.diff(filtered, variable) - filtered_derivative)


def subfilter_flux(
    velocity_pair: tuple[sp.Expr, sp.Expr],
    coordinates: tuple[sp.Symbol, ...],
    time: sp.Symbol,
    width: Any,
    kernel: str = "tophat",
) -> sp.ImmutableMatrix:
    """Exact Pi_ij = S(v_i v_j) - (S v_i)(S v_j) for polynomial fields.

    The input pair is (v_i, v_j) as expressions in coordinates and time;
    filtering is applied per coordinate (product kernel).
    """

    width = sp.sympify(width)
    v_a, v_b = velocity_pair
    filtered_aa = _filter_multivariate(v_a * v_a, coordinates, width, kernel)
    filtered_ab = _filter_multivariate(v_a * v_b, coordinates, width, kernel)
    filtered_bb = _filter_multivariate(v_b * v_b, coordinates, width, kernel)
    mean_a = _filter_multivariate(v_a, coordinates, width, kernel)
    mean_b = _filter_multivariate(v_b, coordinates, width, kernel)
    flux = sp.Matrix(
        [
            [filtered_aa - mean_a * mean_a, filtered_ab - mean_a * mean_b],
            [filtered_ab - mean_a * mean_b, filtered_bb - mean_b * mean_b],
        ]
    )
    return sp.ImmutableMatrix(sp.expand(flux))


def _filter_multivariate(
    expression: Any,
    coordinates: tuple[sp.Symbol, ...],
    width: Any,
    kernel: str,
) -> sp.Expr:
    """Filter an expression coordinate-wise (time dependence untouched)."""

    filtered = sp.sympify(expression)
    for coordinate in coordinates:
        if filtered.has(coordinate):
            filtered = filter_polynomial(filtered, coordinate, width, kernel)
    return filtered


def filtered_balance_residual(
    velocity: tuple[sp.Expr, sp.Expr],
    pressure: sp.Expr,
    body_force: tuple[sp.Expr, sp.Expr],
    density: Any,
    coordinates: tuple[sp.Symbol, ...],
    time: sp.Symbol,
    width: Any,
    kernel: str = "tophat",
) -> list[sp.Expr]:
    """Residual of the exact filtered balance on a polynomial solution class.

    Given microscopic fields satisfying constant-density incompressible
    Euler, the residual of

        d_t ubar_i + d_j(ubar_i ubar_j + Pi_ij) + (1/rho) d_i pbar - fbar_i

    is computed exactly. For convolution filters S commutes with every
    coordinate and time derivative, so the residual is identically zero;
    the returned expressions are the symbolic proof objects.
    """

    if len(velocity) != 2 or len(body_force) != 2 or len(coordinates) != 2:
        raise ValueError("this helper verifies the declared 2-D balance class")

    def mean(value: Any) -> sp.Expr:
        return _filter_multivariate(value, coordinates, width, kernel)

    flux = subfilter_flux(velocity, coordinates, time, width, kernel)
    residual: list[sp.Expr] = []
    for i in range(2):
        inertial = sp.diff(mean(velocity[i]), time)
        transport = sum(
            sp.diff(mean(velocity[i]) * mean(velocity[j]) + flux[i, j], coordinates[j])
            for j in range(2)
        )
        pressure_gradient = sp.diff(mean(pressure), coordinates[i]) / density
        forcing = mean(body_force[i])
        residual.append(
            sp.simplify(inertial + transport + pressure_gradient - forcing)
        )
    return residual


def microscopic_balance_residual(
    velocity: tuple[sp.Expr, sp.Expr],
    pressure: sp.Expr,
    body_force: tuple[sp.Expr, sp.Expr],
    density: Any,
    coordinates: tuple[sp.Symbol, ...],
    time: sp.Symbol,
) -> list[sp.Expr]:
    """Residual of constant-density incompressible Euler for declared fields."""

    residual: list[sp.Expr] = []
    divergence = sum(
        sp.diff(velocity[a], coordinates[a]) for a in range(2)
    )
    if sp.simplify(divergence) != 0:
        raise ValueError("declared velocity field must be divergence free")
    for i in range(2):
        convective = sum(
            velocity[j] * sp.diff(velocity[i], coordinates[j]) for j in range(2)
        )
        residual.append(
            sp.simplify(
                sp.diff(velocity[i], time)
                + convective
                + sp.diff(pressure, coordinates[i]) / density
                - body_force[i]
            )
        )
    return residual


def barotropic_balance_residual(
    velocity: tuple[sp.Expr, sp.Expr],
    density: Any,
    pressure: Any,
    body_force: tuple[sp.Expr, sp.Expr],
    coordinates: tuple[sp.Symbol, ...],
    time: sp.Symbol,
) -> list[sp.Expr]:
    """Residual of compressible barotropic Euler for declared fields.

    The closure ``pressure = p(rho)`` is a *declared premise*: this helper
    never assumes a particular equation of state, it only evaluates the
    momentum balance
    ``d_t v_i + v_j d_j v_i + d_i p(rho)/rho - f_i`` together with the mass
    continuity equation ``d_t rho + d_j(rho v_j)``. Unlike the
    incompressible residual above there is no solenoidality requirement;
    pressure is here a state function of density, not a constraint
    multiplier.
    """

    if len(velocity) != 2 or len(body_force) != 2 or len(coordinates) != 2:
        raise ValueError("this helper verifies the declared 2-D balance class")
    rho = sp.sympify(density)
    p = sp.sympify(pressure)
    residual: list[sp.Expr] = []
    for i in range(2):
        convective = sum(
            velocity[j] * sp.diff(velocity[i], coordinates[j]) for j in range(2)
        )
        residual.append(
            sp.simplify(
                sp.diff(velocity[i], time)
                + convective
                + sp.diff(p, coordinates[i]) / rho
                - body_force[i]
            )
        )
    divergence = sum(sp.diff(rho * velocity[a], coordinates[a]) for a in range(2))
    residual.append(sp.simplify(sp.diff(rho, time) + divergence))
    return residual


def leonard_expansion_residual(
    expression: sp.Expr,
    variable: sp.Symbol,
    width: Any,
    truncation_order: int,
    kernel: str = "tophat",
) -> sp.Expr:
    """Residual after truncating the exact moment series at declared order.

    For polynomial input the series terminates, so at or beyond the
    polynomial degree the residual is identically zero; below it, the
    residual is the exact discarded tail, demonstrating the O(Delta^2)
    structure of the first neglected term.
    """

    moments = kernel_even_moments(kernel, truncation_order)
    series = sp.Integer(0)
    for k in range(truncation_order + 1):
        derivative = sp.diff(expression, variable, 2 * k)
        if derivative != 0:
            series += (
                moments[k]
                * sp.sympify(width) ** (2 * k)
                / sp.factorial(2 * k)
                * derivative
            )
    exact = filter_polynomial(expression, variable, width, kernel)
    return sp.simplify(exact - series)
