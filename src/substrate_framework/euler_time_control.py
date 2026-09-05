"""Finite frequency algebra for the proposed fixed-window Euler response bridge.

The frequency weights approximate a polynomial times a carrier on compact
time intervals. They do not construct Euler quasimodes or an observation
gain. Using them for a physical response requires the same actual generator,
detector and positive material law for all sources, with full norm costs.
"""

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class BandPolynomial:
    """Combined coherent amplitudes, with a conservative finite norm bound."""

    frequencies: tuple[sp.Expr, ...]
    weights: tuple[sp.Expr, ...]
    coefficient_norm_bound: sp.Expr


def finite_band_polynomial(coefficients, carrier, step):
    """Return frequencies/weights for exp(i carrier t) P((exp(i step t)-1)/(i step)).

    Coefficients are ordered from constant term upward. A real nonzero step
    tending to zero gives exp(i carrier t) P(t) in every fixed C^m time norm.
    Choose the step after the finite polynomial degree to retain all returned
    frequencies inside the intended open band. The coefficient norm can grow
    like an inverse power of the step and stays explicit in the result.
    """
    coefficients = tuple(map(sp.sympify, coefficients))
    if not coefficients:
        raise ValueError("at least one polynomial coefficient is required")
    carrier, step = map(sp.sympify, (carrier, step))
    if step.is_zero is True:
        raise ValueError("frequency step must be nonzero")
    if step.is_real is False or carrier.is_real is False:
        raise ValueError("carrier and frequency step must be real")
    count = len(coefficients)
    weights = tuple(sp.simplify(sum(
        coefficients[n]*sp.binomial(n, j)*(-1)**(n-j)/(sp.I*step)**n
        for n in range(j, count))) for j in range(count))
    bound = sum(sp.Abs(value)*(2/sp.Abs(step))**n
                for n, value in enumerate(coefficients))
    return BandPolynomial(tuple(carrier+j*step for j in range(count)),
                          weights, bound)
