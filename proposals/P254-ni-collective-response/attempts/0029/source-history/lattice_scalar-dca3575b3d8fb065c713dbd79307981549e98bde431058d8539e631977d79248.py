"""Exact nearest-neighbour scalar-lattice and continuum-limit ledgers.

The spatial lattice is one-dimensional, uniform, and periodic, with positive
spacing ``a``.  In normalized sine-Gordon units the continuum target is

``phi_tt - phi_xx + m**2*sin(phi) = 0``.

This module keeps four statements separate: the exact finite-spacing stencil,
its local smooth-field Taylor expansion, the Riemann-normalized discrete
action, and modewise long-wave convergence.  It does not derive a lattice,
choose its spacing, or establish convergence of nonlinear solutions.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import sympy as sp


def _positive_exact(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact rather than floating")
    if expression.is_real is not True or expression.is_positive is not True:
        raise ValueError(f"{name} must be provably positive and real")
    return expression


def _positive(value: Any, name: str) -> sp.Expr:
    result = sp.sympify(value)
    if result.is_number and (result.is_real is not True or result.is_positive is not True):
        raise ValueError(f"{name} must be positive and real")
    return result


def _nonnegative(value: Any, name: str) -> sp.Expr:
    result = sp.sympify(value)
    if result.is_number and (
        result.is_real is not True or result.is_nonnegative is not True
    ):
        raise ValueError(f"{name} must be nonnegative and real")
    return result


def _periodic_state(values: Sequence[Any], name: str) -> tuple[sp.Expr, ...]:
    result = tuple(sp.sympify(value) for value in values)
    if len(result) < 2:
        raise ValueError(f"{name} must contain at least two periodic sites")
    return result


def forward_difference(left: Any, right: Any, spacing: Any) -> sp.Expr:
    """Return the oriented nearest-neighbour difference ``(right-left)/a``."""

    a = _positive(spacing, "spacing")
    return sp.simplify((sp.sympify(right) - sp.sympify(left)) / a)


def centered_second_difference(
    left: Any,
    center: Any,
    right: Any,
    spacing: Any,
) -> sp.Expr:
    """Return ``(right - 2*center + left)/a**2``."""

    a = _positive(spacing, "spacing")
    return sp.simplify(
        (sp.sympify(right) - 2 * sp.sympify(center) + sp.sympify(left)) / a**2
    )


def centered_taylor_laplacian(
    derivatives: Sequence[Any],
    spacing: Any,
) -> sp.Expr:
    """Derive the centered stencil from a finite Taylor jet.

    ``derivatives[r]`` denotes the derivative of order ``r`` at the center.
    The highest supplied order must be even and at least two.  Both neighbour
    series are constructed explicitly, so cancellation of odd derivatives and
    every surviving factorial coefficient are derived rather than inserted.
    """

    jet = tuple(sp.sympify(value) for value in derivatives)
    highest = len(jet) - 1
    if highest < 2 or highest % 2:
        raise ValueError("Taylor jet must end at an even order of at least two")
    a = _positive(spacing, "spacing")
    plus = sp.Add(
        *(jet[order] * a**order / sp.factorial(order) for order in range(highest + 1))
    )
    minus = sp.Add(
        *(
            jet[order] * (-a) ** order / sp.factorial(order)
            for order in range(highest + 1)
        )
    )
    return sp.expand((plus - 2 * jet[0] + minus) / a**2)


def centered_taylor_remainder_bound(
    spacing: Any,
    next_even_derivative_bound: Any,
    *,
    retained_derivative_order: int = 6,
) -> sp.Expr:
    """Bound the centered-stencil remainder after an even derivative order.

    If ``phi`` has derivative order ``r+2`` bounded by ``M`` on
    ``[x-a,x+a]`` and the modified equation retains even derivatives through
    order ``r``, the two Lagrange remainders give

    ``abs(remainder) <= 2*M*a**r/(r+2)!``.
    """

    if (
        isinstance(retained_derivative_order, bool)
        or retained_derivative_order < 2
        or retained_derivative_order % 2
    ):
        raise ValueError("retained_derivative_order must be an even integer at least two")
    a = _positive(spacing, "spacing")
    bound = _nonnegative(next_even_derivative_bound, "derivative bound")
    return sp.simplify(
        2
        * bound
        * a**retained_derivative_order
        / sp.factorial(retained_derivative_order + 2)
    )


def lattice_laplacian_symbol(wavenumber: Any, spacing: Any) -> sp.Expr:
    """Return the exact centered-Laplacian Fourier symbol.

    For the mode ``exp(I*k*j*a)`` the symbol is
    ``-4*sin(k*a/2)**2/a**2``.  It is even and periodic under
    ``k -> k + 2*pi/a``; a unique physical representative therefore requires
    a declared first-Brillouin-zone convention.
    """

    k = sp.sympify(wavenumber)
    a = _positive(spacing, "spacing")
    return -4 * sp.sin(k * a / 2) ** 2 / a**2


def lattice_spatial_frequency_squared(wavenumber: Any, spacing: Any) -> sp.Expr:
    """Return the nonnegative exact lattice quantity ``-Delta_a(k)``."""

    return -lattice_laplacian_symbol(wavenumber, spacing)


def linearized_lattice_dispersion_squared(
    wavenumber: Any,
    spacing: Any,
    mass: Any = 1,
) -> sp.Expr:
    """Return ``omega**2=m**2+4*sin(k*a/2)**2/a**2``.

    This is the exact plane-wave dispersion of the linearization about a
    cosine-potential minimum.  It is not the nonlinear dispersion of a finite
    amplitude sine-Gordon solution.
    """

    m = _nonnegative(mass, "mass")
    return sp.simplify(m**2 + lattice_spatial_frequency_squared(wavenumber, spacing))


def lattice_mode_relative_deficit(wavenumber: Any, spacing: Any) -> sp.Expr:
    """Return ``1-kappa_a**2/k**2`` for a nonzero continuum wave number.

    The removable ``k=0`` value is zero.  A concrete zero input is rejected so
    callers cannot silently divide by zero; symbolic callers may take the
    limit explicitly.
    """

    k = sp.sympify(wavenumber)
    if k.is_number and sp.simplify(k) == 0:
        raise ValueError("wavenumber must be nonzero for a relative deficit")
    return sp.simplify(1 - lattice_spatial_frequency_squared(k, spacing) / k**2)


def periodic_lattice_lagrangian(
    field: Sequence[Any],
    velocity: Sequence[Any],
    spacing: Any,
    mass: Any = 1,
) -> sp.Expr:
    """Return the Riemann-normalized periodic lattice Lagrangian.

    For ``N`` sites this is

    ``a*sum_j[dot(phi_j)**2/2 - ((phi_(j+1)-phi_j)/a)**2/2
              - m**2*(1-cos(phi_j))]``.

    The wraparound bond from site ``N-1`` to site zero is included exactly.
    """

    phi = _periodic_state(field, "field")
    speed = _periodic_state(velocity, "velocity")
    if len(phi) != len(speed):
        raise ValueError("field and velocity must have the same number of sites")
    a = _positive(spacing, "spacing")
    m = _nonnegative(mass, "mass")
    density = []
    for index, value in enumerate(phi):
        gradient = forward_difference(value, phi[(index + 1) % len(phi)], a)
        density.append(
            speed[index] ** 2 / 2
            - gradient**2 / 2
            - m**2 * (1 - sp.cos(value))
        )
    return sp.simplify(a * sp.Add(*density))


def periodic_lattice_eom_residual(
    field: Sequence[Any],
    acceleration: Sequence[Any],
    spacing: Any,
    mass: Any = 1,
) -> tuple[sp.Expr, ...]:
    """Return every exact periodic discrete sine-Gordon EOM residual.

    Each entry is ``ddot(phi_j)-Delta_a(phi)_j+m**2*sin(phi_j)`` and follows
    from :func:`periodic_lattice_lagrangian` by sitewise variation.
    """

    phi = _periodic_state(field, "field")
    accel = _periodic_state(acceleration, "acceleration")
    if len(phi) != len(accel):
        raise ValueError("field and acceleration must have the same number of sites")
    a = _positive(spacing, "spacing")
    m = _nonnegative(mass, "mass")
    return tuple(
        sp.simplify(
            accel[index]
            - centered_second_difference(
                phi[(index - 1) % len(phi)],
                value,
                phi[(index + 1) % len(phi)],
                a,
            )
            + m**2 * sp.sin(value)
        )
        for index, value in enumerate(phi)
    )


def periodic_action_error_bound(
    length: Any,
    time_duration: Any,
    spacing: Any,
    max_abs_phi_x: Any,
    max_abs_phi_xx: Any,
    max_abs_phi_t: Any,
    max_abs_phi_tx: Any,
    mass: Any = 1,
) -> sp.Expr:
    """Return a sufficient sampled-action error bound for a smooth field.

    Let ``a=L/N`` and sample an ``L``-periodic field on left endpoints.  Assume
    the four supplied uniform derivative bounds hold on the space-time
    cylinder.  The left-Riemann error and the forward-gradient Taylor error
    give the absolute action difference bound

    ``T*L*(a*Mt*Mtx/2 + a*m**2*Mx/2 + a*Mx*Mxx
            + a**2*Mxx**2/8)``.

    Hence the sampled discrete action converges to the continuum action at
    least linearly as ``a -> 0`` under these fixed bounds.  This is an action
    statement, not a nonlinear solution-convergence theorem.
    """

    interval = _positive(length, "length")
    duration = _positive(time_duration, "time duration")
    a = _positive(spacing, "spacing")
    mx = _nonnegative(max_abs_phi_x, "max_abs_phi_x")
    mxx = _nonnegative(max_abs_phi_xx, "max_abs_phi_xx")
    mt = _nonnegative(max_abs_phi_t, "max_abs_phi_t")
    mtx = _nonnegative(max_abs_phi_tx, "max_abs_phi_tx")
    m = _nonnegative(mass, "mass")
    instantaneous = (
        a * mt * mtx / 2
        + a * m**2 * mx / 2
        + a * mx * mxx
        + a**2 * mxx**2 / 8
    )
    return sp.simplify(duration * interval * instantaneous)


@dataclass(frozen=True)
class PhysicalPhaseChainCoefficients:
    """Exact positive coefficients of a dimensionless-phase chain.

    The per-site Lagrangian is

    ``I*dot(u_j)**2/2-K*(u_(j+1)-u_j)**2/2-V0*(1-cos(u_j))``.

    Thus ``I`` is a phase inertia with dimensions energy times time squared;
    it is not a bare mass. ``K`` and ``V0`` are energies and ``a`` is the
    physical site spacing.
    """

    inertia: sp.Expr
    coupling: sp.Expr
    onsite: sp.Expr
    spacing: sp.Expr


@dataclass(frozen=True)
class PhysicalPhaseChainScales:
    """Exact linear scales of a declared physical phase chain."""

    gap_frequency: sp.Expr
    band_edge_frequency: sp.Expr
    long_wave_speed: sp.Expr


def physical_phase_chain_coefficients(
    inertia: Any,
    coupling: Any,
    onsite: Any,
    spacing: Any,
) -> PhysicalPhaseChainCoefficients:
    """Validate and return physical phase-chain coefficients."""

    return PhysicalPhaseChainCoefficients(
        inertia=_positive_exact(inertia, "inertia"),
        coupling=_positive_exact(coupling, "coupling"),
        onsite=_positive_exact(onsite, "onsite"),
        spacing=_positive_exact(spacing, "spacing"),
    )


def _validated_physical_phase_chain(
    coefficients: PhysicalPhaseChainCoefficients,
) -> PhysicalPhaseChainCoefficients:
    if not isinstance(coefficients, PhysicalPhaseChainCoefficients):
        raise TypeError(
            "coefficients must be a PhysicalPhaseChainCoefficients record"
        )
    return physical_phase_chain_coefficients(
        coefficients.inertia,
        coefficients.coupling,
        coefficients.onsite,
        coefficients.spacing,
    )


def phase_inertia_from_mass_scale(mass: Any, displacement_scale: Any) -> sp.Expr:
    """Return ``I=m*b**2`` after declaring the displacement ``q=b*u``.

    This is the load-bearing coordinate conversion omitted by a formula such
    as ``sqrt(V0/m)`` when ``u`` is dimensionless and ``V0`` is an energy.
    """

    m = _positive_exact(mass, "mass")
    b = _positive_exact(displacement_scale, "displacement_scale")
    return m * b**2


def phase_coupling_from_stiffness_scale(
    displacement_stiffness: Any,
    displacement_scale: Any,
) -> sp.Expr:
    """Return the phase coupling ``K=kappa*b**2`` for ``q=b*u``."""

    stiffness = _positive_exact(
        displacement_stiffness, "displacement_stiffness"
    )
    b = _positive_exact(displacement_scale, "displacement_scale")
    return stiffness * b**2


def physical_phase_chain_from_displacement(
    mass: Any,
    displacement_scale: Any,
    displacement_stiffness: Any,
    onsite: Any,
    spacing: Any,
) -> PhysicalPhaseChainCoefficients:
    """Lift a physical displacement chain ``q=b*u`` into phase variables."""

    return physical_phase_chain_coefficients(
        phase_inertia_from_mass_scale(mass, displacement_scale),
        phase_coupling_from_stiffness_scale(
            displacement_stiffness, displacement_scale
        ),
        onsite,
        spacing,
    )


def periodic_physical_phase_chain_lagrangian(
    field: Sequence[Any],
    velocity: Sequence[Any],
    coefficients: PhysicalPhaseChainCoefficients,
) -> sp.Expr:
    """Return the exact periodic per-site-energy phase-chain Lagrangian."""

    phi = _periodic_state(field, "field")
    speed = _periodic_state(velocity, "velocity")
    if len(phi) != len(speed):
        raise ValueError("field and velocity must have the same number of sites")
    model = _validated_physical_phase_chain(coefficients)
    return sp.simplify(
        sp.Add(
            *(
                model.inertia * speed[index] ** 2 / 2
                - model.coupling
                * (phi[(index + 1) % len(phi)] - value) ** 2
                / 2
                - model.onsite * (1 - sp.cos(value))
                for index, value in enumerate(phi)
            )
        )
    )


def periodic_physical_phase_chain_eom_residual(
    field: Sequence[Any],
    acceleration: Sequence[Any],
    coefficients: PhysicalPhaseChainCoefficients,
) -> tuple[sp.Expr, ...]:
    """Return ``I*u_j,tt-K*(u_(j+1)-2u_j+u_(j-1))+V0*sin(u_j)``."""

    phi = _periodic_state(field, "field")
    accel = _periodic_state(acceleration, "acceleration")
    if len(phi) != len(accel):
        raise ValueError("field and acceleration must have the same number of sites")
    model = _validated_physical_phase_chain(coefficients)
    return tuple(
        sp.simplify(
            model.inertia * accel[index]
            - model.coupling
            * (
                phi[(index + 1) % len(phi)]
                - 2 * value
                + phi[(index - 1) % len(phi)]
            )
            + model.onsite * sp.sin(value)
        )
        for index, value in enumerate(phi)
    )


def physical_phase_chain_dispersion_squared(
    wavenumber: Any,
    coefficients: PhysicalPhaseChainCoefficients,
) -> sp.Expr:
    r"""Return ``(V0+4*K*sin(k*a/2)**2)/I`` for the vacuum linearization."""

    k = sp.sympify(wavenumber)
    if k.has(sp.Float):
        raise ValueError("wavenumber must be exact rather than floating")
    if k.is_real is not True:
        raise ValueError("wavenumber must be provably real")
    model = _validated_physical_phase_chain(coefficients)
    return sp.simplify(
        (
            model.onsite
            + 4 * model.coupling * sp.sin(k * model.spacing / 2) ** 2
        )
        / model.inertia
    )


def physical_phase_chain_scales(
    coefficients: PhysicalPhaseChainCoefficients,
) -> PhysicalPhaseChainScales:
    """Return the gap, first-zone band edge, and long-wave signal speed."""

    model = _validated_physical_phase_chain(coefficients)
    return PhysicalPhaseChainScales(
        gap_frequency=sp.sqrt(model.onsite / model.inertia),
        band_edge_frequency=sp.sqrt(
            (model.onsite + 4 * model.coupling) / model.inertia
        ),
        long_wave_speed=model.spacing * sp.sqrt(
            model.coupling / model.inertia
        ),
    )


def physical_phase_chain_gap_ratio(
    numerator: PhysicalPhaseChainCoefficients,
    denominator: PhysicalPhaseChainCoefficients,
) -> sp.Expr:
    """Return ``omega_num/omega_den`` with both curvatures and inertias."""

    first = _validated_physical_phase_chain(numerator)
    second = _validated_physical_phase_chain(denominator)
    return sp.sqrt(
        first.onsite * second.inertia / (first.inertia * second.onsite)
    )


def mass_scale_phase_gap_ratio(
    numerator_mass: Any,
    numerator_scale: Any,
    numerator_onsite: Any,
    denominator_mass: Any,
    denominator_scale: Any,
    denominator_onsite: Any,
) -> sp.Expr:
    r"""Return the exact gap ratio after both lifts ``I=m*b**2``.

    The result is
    ``sqrt(V_num*m_den*b_den**2/(V_den*m_num*b_num**2))``. It becomes
    ``sqrt(m_den/m_num)`` only when the on-site energies and displacement
    scales are equal.
    """

    first_inertia = phase_inertia_from_mass_scale(
        numerator_mass, numerator_scale
    )
    second_inertia = phase_inertia_from_mass_scale(
        denominator_mass, denominator_scale
    )
    first_onsite = _positive_exact(numerator_onsite, "numerator_onsite")
    second_onsite = _positive_exact(denominator_onsite, "denominator_onsite")
    return sp.sqrt(
        first_onsite * second_inertia / (first_inertia * second_onsite)
    )


def physical_phase_chain_dimension_matrix() -> sp.ImmutableMatrix:
    r"""Return dimensions over rows ``(mass,length,time)``.

    Columns are ``(I,K,V0,a,m,b)``. In particular, ``I`` has dimensions
    ``M*L**2`` while a bare mass has dimensions ``M``.
    """

    return sp.ImmutableMatrix(
        [
            [1, 1, 1, 0, 1, 0],
            [2, 2, 2, 1, 0, 1],
            [0, -2, -2, 0, 0, 0],
        ]
    )
