"""Exact conditional phase-ensemble and activation-threshold algebra.

These utilities distinguish an expected directional quadratic observable from
total energy and a declared Brownian phase law from a physical noise model.
They do not derive a physical emitter population, phase-diffusion coefficient,
normalization, effective temperature, stochastic escape process, material
barrier, event channel, or event payload.
"""

from __future__ import annotations

from typing import Any

import sympy as sp


def _exact_real(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact rather than floating")
    if expression.is_real is not True:
        raise ValueError(f"{name} must be provably real")
    return expression


def _exact_positive(value: Any, name: str) -> sp.Expr:
    expression = _exact_real(value, name)
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be provably positive")
    return expression


def _exact_nonnegative(value: Any, name: str) -> sp.Expr:
    expression = _exact_real(value, name)
    if expression.is_nonnegative is not True:
        raise ValueError(f"{name} must be provably nonnegative")
    return expression


def _exact_integer(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact rather than floating")
    if expression.is_integer is not True:
        raise ValueError(f"{name} must be a provable integer")
    return expression


def _unit_interval(value: Any, name: str) -> sp.Expr:
    expression = _exact_real(value, name)
    if expression.is_number and (
        bool(expression < 0) or bool(expression > 1)
    ):
        raise ValueError(f"{name} must lie in the closed unit interval")
    return expression


def iid_equal_amplitude_expected_intensity(
    count: Any,
    per_source_intensity: Any,
    pair_coherence: Any,
) -> sp.Expr:
    r"""Return the expected directional intensity of iid unit phasors.

    For positive integer ``N``, per-source intensity ``I1>0``, iid unit
    phasors ``z_j`` with common mean ``mu``, and
    ``V=abs(mu)**2`` supplied through ``pair_coherence``, independence gives

    ``I1*E[abs(sum_j z_j)**2] = I1*(N + N*(N-1)*V)``.

    The result is a directional quadratic observable. Integrating over
    direction, changing the per-source normalization, or interpreting it as
    thermal or population energy requires a separate model.
    """

    population = sp.sympify(count)
    if population.has(sp.Float):
        raise ValueError("count must be exact rather than floating")
    if population.is_integer is not True or population.is_positive is not True:
        raise ValueError("count must be a provably positive integer")
    intensity = _exact_positive(per_source_intensity, "per_source_intensity")
    coherence = _unit_interval(pair_coherence, "pair_coherence")
    return sp.simplify(
        intensity
        * (population + population * (population - 1) * coherence)
    )


def gaussian_phase_pair_coherence(phase_variance: Any) -> sp.Expr:
    r"""Return ``exp(-variance)`` for iid centered Gaussian phase noise.

    A centered Gaussian phase with variance ``sigma**2`` has mean phasor
    ``exp(-sigma**2/2)`` and hence pair coherence
    ``abs(E[exp(i*phase)])**2=exp(-sigma**2)``. The Gaussian phase law is a
    declared ensemble premise, not a dynamical decoherence derivation.
    """

    variance = _exact_real(phase_variance, "phase_variance")
    if variance.is_nonnegative is not True:
        raise ValueError("phase_variance must be provably nonnegative")
    return sp.exp(-variance)


def brownian_phase_characteristic(
    diffusion: Any,
    elapsed_time: Any,
    harmonic: Any = 1,
) -> sp.Expr:
    r"""Return the characteristic mode of a declared Brownian phase.

    For ``delta_t=sqrt(2*D)*W_t`` with ``D>=0``, ``t>=0``, and integer
    harmonic ``n``, the exact characteristic function is
    ``E[exp(i*n*delta_t)]=exp(-n**2*D*t)``. This conditional stochastic law
    does not derive ``D`` or map the phase to a physical oscillator.
    """

    phase_diffusion = _exact_nonnegative(diffusion, "diffusion")
    duration = _exact_nonnegative(elapsed_time, "elapsed_time")
    mode = _exact_integer(harmonic, "harmonic")
    return sp.exp(-(mode**2) * phase_diffusion * duration)


def brownian_phase_pair_coherence(diffusion: Any, elapsed_time: Any) -> sp.Expr:
    r"""Return iid pair coherence ``exp(-2*D*t)`` at a declared time.

    One Brownian phase has variance ``2*D*t`` and mean phasor
    ``exp(-D*t)``. The iid equal-time pair coherence is the squared mean
    phasor, ``exp(-2*D*t)``. Neither factor is a survival probability.
    """

    phase_diffusion = _exact_nonnegative(diffusion, "diffusion")
    duration = _exact_nonnegative(elapsed_time, "elapsed_time")
    return gaussian_phase_pair_coherence(2 * phase_diffusion * duration)


def _uniform_exponential_window_average(rate: sp.Expr, window: sp.Expr) -> sp.Expr:
    if rate.is_zero is True:
        return sp.Integer(1)
    average = (1 - sp.exp(-rate * window)) / (rate * window)
    if rate.is_positive is True:
        return average
    return sp.Piecewise((1, sp.Eq(rate, 0)), (average, True))


def brownian_mean_phasor_window_average(diffusion: Any, window: Any) -> sp.Expr:
    r"""Return the uniform-window average of ``exp(-D*t)``.

    For ``T>0`` this is ``(1-exp(-D*T))/(D*T)`` when ``D>0`` and one when
    ``D=0``. It differs from the endpoint mean phasor ``exp(-D*T)``.
    """

    phase_diffusion = _exact_nonnegative(diffusion, "diffusion")
    duration = _exact_positive(window, "window")
    return _uniform_exponential_window_average(phase_diffusion, duration)


def brownian_pair_coherence_window_average(diffusion: Any, window: Any) -> sp.Expr:
    r"""Return the uniform-window average of pair coherence ``exp(-2*D*t)``."""

    phase_diffusion = _exact_nonnegative(diffusion, "diffusion")
    duration = _exact_positive(window, "window")
    return _uniform_exponential_window_average(2 * phase_diffusion, duration)


def damped_brownian_coherent_mean_factor(
    damping: Any,
    diffusion: Any,
    elapsed_time: Any,
) -> sp.Expr:
    r"""Return ``exp(-(Gamma/2+D)*t)`` for declared independent factors.

    The first factor is a deterministic coordinate-amplitude envelope and the
    second is a Brownian mean phasor. Their product is a coherent mean-field
    factor, not a survival or population fraction.
    """

    loss = _exact_nonnegative(damping, "damping")
    phase_diffusion = _exact_nonnegative(diffusion, "diffusion")
    duration = _exact_nonnegative(elapsed_time, "elapsed_time")
    return sp.exp(-(loss / 2 + phase_diffusion) * duration)


def damped_brownian_coherent_quadratic_factor(
    damping: Any,
    diffusion: Any,
    elapsed_time: Any,
) -> sp.Expr:
    r"""Return the quadratic factor ``exp(-(Gamma+2*D)*t)``.

    This is the square of :func:`damped_brownian_coherent_mean_factor` under
    the same declared independent deterministic-envelope and Brownian-phase
    assumptions. It is not a total energy or survival probability.
    """

    mean_factor = damped_brownian_coherent_mean_factor(
        damping,
        diffusion,
        elapsed_time,
    )
    return sp.simplify(mean_factor**2)


def population_activation_scale(
    population: Any,
    per_source_scale: Any,
    coherence: Any,
) -> sp.Expr:
    r"""Return the declared interpolation ``theta*n*(1+(n-1)*V)``.

    Here ``n>0``, ``theta>0``, and ``V`` lies in ``[0,1]``. The function is
    the continuous-coordinate extension of the iid equal-amplitude expected
    intensity. Calling it energy, temperature, or an activation scale is a
    separate physical premise.
    """

    coordinate = _exact_positive(population, "population")
    scale = _exact_positive(per_source_scale, "per_source_scale")
    visibility = _unit_interval(coherence, "coherence")
    return sp.simplify(
        scale * coordinate * (1 + (coordinate - 1) * visibility)
    )


def continuous_population_threshold(
    barrier: Any,
    per_source_scale: Any,
    coherence: Any,
) -> sp.Expr:
    r"""Return the positive continuous population solving ``Theta=barrier``.

    With ``x=barrier/per_source_scale``, the incoherent branch ``V=0`` is
    ``n=x``. For ``V>0`` the unique positive root is

    ``(sqrt((1-V)**2+4*V*x) - (1-V))/(2*V)``.

    A physical integer count requires taking the ceiling only after a model
    declares the continuous interpolation and its threshold convention.
    """

    activation_barrier = _exact_positive(barrier, "barrier")
    scale = _exact_positive(per_source_scale, "per_source_scale")
    visibility = _unit_interval(coherence, "coherence")
    ratio = sp.simplify(activation_barrier / scale)
    if visibility.is_zero is True:
        return ratio
    positive_branch = sp.simplify(
        (
            sp.sqrt((1 - visibility) ** 2 + 4 * visibility * ratio)
            - (1 - visibility)
        )
        / (2 * visibility)
    )
    if visibility.is_zero is False:
        return positive_branch
    return sp.Piecewise(
        (ratio, sp.Eq(visibility, 0)),
        (positive_branch, True),
    )


def activated_relative_response(barrier: Any, activation_scale: Any) -> sp.Expr:
    r"""Return the dimensionless activated factor ``exp(-E/Theta)``.

    This factor is not a rate. A rate additionally needs a dimensionful
    prefactor and a justified stochastic or kinetic model.
    """

    activation_barrier = _exact_positive(barrier, "barrier")
    scale = _exact_positive(activation_scale, "activation_scale")
    return sp.exp(-activation_barrier / scale)
