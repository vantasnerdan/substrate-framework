"""Exact mixed-coordinate coefficients of a declared cosine potential.

These helpers perform classical local calculus on
``amplitude*(1-cos(background+high_scale*H+low_scale*L))``.  The coordinate
names do not supply a frequency split, quantization, phonon dictionary,
matrix element, kinematic channel, transition rate, or material realization.
"""

from __future__ import annotations

from typing import Any

import sympy as sp


def _nonnegative_integer(value: Any, *, name: str) -> int:
    expression = sp.sympify(value)
    if (
        expression.is_number is not True
        or expression.is_integer is not True
        or expression.is_nonnegative is not True
    ):
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(expression)


def cosine_mixed_coefficient(
    high_order: Any,
    low_order: Any,
    *,
    amplitude: Any = 1,
    high_scale: Any = 1,
    low_scale: Any = 1,
    background: Any = 0,
) -> sp.Expr:
    """Return the coefficient of ``H**high_order*L**low_order``.

    The expanded function is
    ``amplitude*(1-cos(background+high_scale*H+low_scale*L))``.
    Polynomial coefficients equal origin mixed derivatives divided by
    ``high_order!*low_order!``.  All normalization factors and the expansion
    background remain explicit.
    """

    j = _nonnegative_integer(high_order, name="high_order")
    k = _nonnegative_integer(low_order, name="low_order")
    amplitude_expression = sp.sympify(amplitude)
    high_normalization = sp.sympify(high_scale)
    low_normalization = sp.sympify(low_scale)
    background_expression = sp.sympify(background)
    argument = sp.Symbol("_cosine_vertex_argument", real=True)
    potential = amplitude_expression * (1 - sp.cos(argument))
    derivative = sp.diff(potential, argument, j + k).subs(
        argument, background_expression
    )
    return sp.simplify(
        derivative
        * high_normalization**j
        * low_normalization**k
        / (sp.factorial(j) * sp.factorial(k))
    )


def cosine_mixed_derivative(
    high_order: Any,
    low_order: Any,
    **kwargs: Any,
) -> sp.Expr:
    """Return the raw mixed derivative for the declared cosine split.

    This differs from :func:`cosine_mixed_coefficient` by the explicit
    factorial ``high_order!*low_order!``.
    """

    j = _nonnegative_integer(high_order, name="high_order")
    k = _nonnegative_integer(low_order, name="low_order")
    return sp.simplify(
        sp.factorial(j)
        * sp.factorial(k)
        * cosine_mixed_coefficient(j, k, **kwargs)
    )


def vacuum_cosine_mixed_coefficient(
    high_order: Any,
    low_order: Any,
    *,
    amplitude: Any = 1,
    high_scale: Any = 1,
    low_scale: Any = 1,
) -> sp.Expr:
    """Return the exact zero-background mixed coefficient.

    The constant coefficient is zero.  For positive total order ``m=j+k``,
    odd ``m`` gives zero, while even ``m`` gives
    ``amplitude*(-1)**(m/2+1)*high_scale**j*low_scale**k/(j!*k!)``.
    """

    return cosine_mixed_coefficient(
        high_order,
        low_order,
        amplitude=amplitude,
        high_scale=high_scale,
        low_scale=low_scale,
        background=0,
    )


def vacuum_one_high_coefficient(
    low_order: Any,
    *,
    amplitude: Any = 1,
    high_scale: Any = 1,
    low_scale: Any = 1,
) -> sp.Expr:
    """Return the zero-background coefficient of ``H*L**low_order``.

    It vanishes for even ``low_order``.  For odd ``n`` it is
    ``amplitude*(-1)**((n-1)/2)*high_scale*low_scale**n/n!``.  A nonzero
    classical coefficient is not by itself a quantum transition amplitude.
    """

    n = _nonnegative_integer(low_order, name="low_order")
    return vacuum_cosine_mixed_coefficient(
        1,
        n,
        amplitude=amplitude,
        high_scale=high_scale,
        low_scale=low_scale,
    )


def cosine_mixed_taylor_polynomial(
    high_field: Any,
    low_field: Any,
    total_order: Any,
    *,
    amplitude: Any = 1,
    high_scale: Any = 1,
    low_scale: Any = 1,
    background: Any = 0,
) -> sp.Expr:
    """Return the total-degree Taylor polynomial through ``total_order``.

    The returned finite polynomial is not the full cosine.  Callers making an
    approximation claim must control the separately omitted remainder on their
    declared field domain.
    """

    order = _nonnegative_integer(total_order, name="total_order")
    high_expression = sp.sympify(high_field)
    low_expression = sp.sympify(low_field)
    polynomial = sp.Integer(0)
    for total in range(order + 1):
        for high_power in range(total + 1):
            low_power = total - high_power
            coefficient = cosine_mixed_coefficient(
                high_power,
                low_power,
                amplitude=amplitude,
                high_scale=high_scale,
                low_scale=low_scale,
                background=background,
            )
            polynomial += (
                coefficient
                * high_expression**high_power
                * low_expression**low_power
            )
    return sp.expand(polynomial)


def cosine_quadratic_gap(argument: Any) -> sp.Expr:
    """Return the exact gap between the quadratic and cosine potentials.

    For a real coordinate ``x`` this is
    ``x**2/2 - (1-cos(x))``.  It is globally nonnegative and no larger than
    :func:`cosine_quadratic_gap_bound`.  The function returns the exact
    expression; it does not silently choose an approximation domain.
    """

    coordinate = sp.sympify(argument)
    return sp.simplify(coordinate**2 / 2 - (1 - sp.cos(coordinate)))


def cosine_quadratic_gap_bound(argument: Any) -> sp.Expr:
    """Return the global fourth-order upper bound for the quadratic gap.

    For every real ``x``,
    ``0 <= cosine_quadratic_gap(x) <= x**4/24``.
    """

    coordinate = sp.sympify(argument)
    return coordinate**4 / sp.factorial(4)


def harmonic_cycle_mean_square(peak_amplitude: Any) -> sp.Expr:
    """Return the full-cycle mean square for a harmonic phase.

    If ``phi(t)=P*cos(omega*t+delta)`` with nonzero real ``omega`` and real
    peak amplitude ``P``, then the mean of ``phi(t)**2`` over one complete
    period is ``P**2/2``.  Calling ``P`` an RMS amplitude would therefore
    change the convention by a factor of ``sqrt(2)``.
    """

    peak = sp.sympify(peak_amplitude)
    return sp.simplify(peak**2 / 2)


def harmonic_rms_from_peak(peak_amplitude: Any) -> sp.Expr:
    """Return the nonnegative RMS amplitude of a real harmonic phase."""

    peak = sp.sympify(peak_amplitude)
    return sp.Abs(peak) / sp.sqrt(2)


def sufficient_cosine_quadratic_domain(relative_tolerance: Any) -> sp.Expr:
    """Return a sufficient symmetric domain radius for a relative error.

    The relative error is measured against the quadratic potential.  For
    ``x != 0`` it is at most ``x**2/12``; therefore a positive tolerance
    ``epsilon`` is guaranteed on ``|x| <= sqrt(12*epsilon)``.  A concrete
    nonpositive tolerance is rejected.  Symbolic callers retain the explicit
    obligation that their tolerance be positive.
    """

    tolerance = sp.sympify(relative_tolerance)
    if tolerance.is_number and tolerance.is_positive is not True:
        raise ValueError("relative_tolerance must be positive")
    return sp.sqrt(12 * tolerance)
