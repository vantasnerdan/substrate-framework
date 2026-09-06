"""Exact symbol tests for an Euler-to-Maxwell emergence proposal.

The helpers expose a route boundary.  They do not assert that Maxwell theory,
charge, or a particle emerges from Euler.
"""

from dataclasses import dataclass
from typing import Any, Sequence

import sympy as sp


@dataclass(frozen=True)
class UniformSymbolLedger:
    """Characteristic data for physical transverse Euler and Maxwell modes."""

    euler_characteristic: sp.Expr
    maxwell_characteristic: sp.Expr
    euler_temporal_frequencies: tuple[sp.Expr, sp.Expr]
    maxwell_temporal_frequencies: tuple[sp.Expr, sp.Expr]


def uniform_symbol_ledger(
    spectral_parameter: Any,
    wave_number: Any,
    advection_frequency: Any,
    wave_speed: Any,
) -> UniformSymbolLedger:
    """Return the two-polarization characteristic polynomials.

    ``spectral_parameter`` is the generator variable ``lambda``.  The Euler
    generator on either physical polarization is ``-i*advection_frequency``;
    the Maxwell generator has roots ``+/- i*wave_speed*wave_number``.
    """

    lam = sp.sympify(spectral_parameter)
    k = sp.sympify(wave_number)
    adv = sp.sympify(advection_frequency)
    speed = sp.sympify(wave_speed)
    if speed.is_positive is False:
        raise ValueError("wave_speed must be positive")
    euler = sp.expand((lam + sp.I * adv) ** 2)
    maxwell = sp.expand((lam**2 + speed**2 * k**2) ** 2)
    return UniformSymbolLedger(
        euler_characteristic=euler,
        maxwell_characteristic=maxwell,
        euler_temporal_frequencies=(adv, adv),
        maxwell_temporal_frequencies=(-speed * k, speed * k),
    )


def transverse_scalar_source(wavevector: Sequence[Any], scalar_amplitude: Any):
    """Project the local isotropic scalar source ``i*k*q`` transversely."""

    k = sp.ImmutableMatrix(wavevector)
    if k.shape != (3, 1):
        raise ValueError("wavevector must have three components")
    norm_squared = sp.simplify(k.dot(k))
    if norm_squared == 0:
        raise ValueError("wavevector must be nonzero")
    projector = sp.eye(3) - (k * k.T) / norm_squared
    source = sp.I * k * sp.sympify(scalar_amplitude)
    return sp.ImmutableMatrix(projector.applyfunc(sp.simplify)), sp.ImmutableMatrix(
        (projector * source).applyfunc(sp.simplify)
    )


def clebsch_shift_kernel(
    beta_gradient: Sequence[Any], function_derivative: Any
) -> sp.ImmutableMatrix:
    """Return ``delta(dphi+alpha dbeta)`` for a Clebsch gauge shift.

    The shift is ``delta phi=-F(beta)``, ``delta alpha=F'(beta)``, and
    ``delta beta=0``.  Its physical velocity variation vanishes identically.
    """

    grad_beta = sp.ImmutableMatrix(beta_gradient)
    if grad_beta.shape != (3, 1):
        raise ValueError("beta_gradient must have three components")
    derivative = sp.sympify(function_derivative)
    variation = -derivative * grad_beta + derivative * grad_beta
    return sp.ImmutableMatrix(variation.applyfunc(sp.simplify))
