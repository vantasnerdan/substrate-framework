"""Exact radial energies and conditional capillary constitutive ledgers.

The Frank/core and quadratic-loading helpers expose declared constitutive
premises.  They do not derive a material, coupling, amplitude convention,
dispersion relation, or nucleation rate.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


def _positive_quantity(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    return expression


def _real_quantity(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_real is not True:
        raise ValueError(f"{name} must be real")
    return expression


def _nonzero_real_quantity(value: Any, name: str) -> sp.Expr:
    expression = _real_quantity(value, name)
    if expression.is_number and expression.is_zero is not False:
        raise ValueError(f"{name} must be nonzero")
    return expression


@dataclass(frozen=True)
class CapillaryConstitutiveResult:
    """Exact result of composing declared tension and area-drive premises."""

    line_tension: sp.Expr
    area_drive: sp.Expr
    critical_radius: sp.Expr
    barrier_height: sp.Expr


@dataclass(frozen=True)
class CapillaryDimensionLedger:
    """Energy/length exponents for a monomial loading-law convention."""

    base_dimensions: tuple[str, str]
    quantity_names: tuple[str, ...]
    dimension_matrix: sp.ImmutableMatrix
    amplitude_length_exponent: sp.Expr
    amplitude_power: sp.Expr
    wavenumber_power: sp.Expr
    coupling_length_exponent: sp.Expr


@dataclass(frozen=True)
class CapillaryIdentifiabilityLedger:
    """Log-exponent diagnostics for the effective quadratic capillary map."""

    parameter_names: tuple[str, ...]
    observable_names: tuple[str, ...]
    exponent_matrix: sp.ImmutableMatrix
    rank: int
    nullspace: tuple[sp.ImmutableMatrix, ...]
    coordinate_identifiable: tuple[bool, ...]
    barrier_only_matrix: sp.ImmutableMatrix
    barrier_only_rank: int
    barrier_only_nullspace: tuple[sp.ImmutableMatrix, ...]
    barrier_only_coordinate_identifiable: tuple[bool, ...]


@dataclass(frozen=True)
class FrankCoreSensitivityLedger:
    """State-dependent logarithmic sensitivities of a Frank/core tension."""

    parameter_names: tuple[str, ...]
    line_tension: sp.Expr
    elastic_line_tension: sp.Expr
    line_tension_log_elasticities: tuple[sp.Expr, ...]
    critical_radius_log_elasticities: tuple[sp.Expr, ...]
    barrier_height_log_elasticities: tuple[sp.Expr, ...]


def line_energy(radius: Any, line_density: Any) -> sp.Expr:
    """Return circumference-weighted line energy ``2*pi*R*lambda``."""

    radial_coordinate = _positive_quantity(radius, "radius")
    density = _positive_quantity(line_density, "line_density")
    return 2 * sp.pi * radial_coordinate * density


def spherical_shell_energy(radius: Any, surface_density: Any) -> sp.Expr:
    """Return spherical-area energy ``4*pi*R**2*sigma``."""

    radial_coordinate = _positive_quantity(radius, "radius")
    density = _positive_quantity(surface_density, "surface_density")
    return 4 * sp.pi * radial_coordinate**2 * density


def capillary_energy(
    radius: Any,
    line_tension: Any,
    pressure: Any,
    core_energy: Any = 0,
) -> sp.Expr:
    """Return ``2*pi*R*T - pi*R**2*P + E_core`` for positive ``R,T,P``."""

    radial_coordinate = _positive_quantity(radius, "radius")
    tension = _positive_quantity(line_tension, "line_tension")
    drive = _positive_quantity(pressure, "pressure")
    core = _real_quantity(core_energy, "core_energy")
    return 2 * sp.pi * radial_coordinate * tension - sp.pi * radial_coordinate**2 * drive + core


def capillary_critical_radius(line_tension: Any, pressure: Any) -> sp.Expr:
    """Return the unique stationary radius ``T/P`` of the capillary energy."""

    tension = _positive_quantity(line_tension, "line_tension")
    drive = _positive_quantity(pressure, "pressure")
    return tension / drive


def capillary_barrier_height(line_tension: Any, pressure: Any) -> sp.Expr:
    """Return the barrier height relative to ``R=0``, ``pi*T**2/P``.

    A radius-independent core-energy offset cancels from this relative height;
    it remains present in the absolute capillary energy.
    """

    tension = _positive_quantity(line_tension, "line_tension")
    drive = _positive_quantity(pressure, "pressure")
    return sp.pi * tension**2 / drive


def frank_core_line_tension(
    frank_constant: Any,
    defect_strength: Any,
    outer_cutoff: Any,
    core_cutoff: Any,
    core_line_energy: Any,
) -> sp.Expr:
    r"""Return the declared Frank/core line tension.

    The premise is
    ``T=pi*K_F*s**2*log(R_o/r_c)+epsilon_core`` with positive
    energy-per-length inputs, real dimensionless ``s``, and ``R_o>r_c>0``.
    Numeric inputs that violate the cutoff ordering are rejected.  The helper
    composes the premise; it does not establish a nematic material model.
    """

    stiffness = _positive_quantity(frank_constant, "frank_constant")
    strength = _real_quantity(defect_strength, "defect_strength")
    outer = _positive_quantity(outer_cutoff, "outer_cutoff")
    core = _positive_quantity(core_cutoff, "core_cutoff")
    core_energy = _positive_quantity(core_line_energy, "core_line_energy")
    ratio = sp.simplify(outer / core)
    if ratio.is_number and bool(ratio <= 1):
        raise ValueError("outer_cutoff must exceed core_cutoff")
    return sp.pi * stiffness * strength**2 * sp.log(ratio) + core_energy


def quadratic_loading_area_drive(
    coupling: Any,
    amplitude: Any,
    wavenumber: Any,
    thickness: Any,
) -> sp.Expr:
    r"""Return the declared quadratic area drive ``g*A**2*k**2*l_m/2``.

    Positive ``g`` and ``l_m`` plus nonzero real ``A`` and ``k`` give a
    positive area drive.  The quadratic loading law is an input premise, not a
    consequence of dimensional analysis.
    """

    coupling_value = _positive_quantity(coupling, "coupling")
    amplitude_value = _nonzero_real_quantity(amplitude, "amplitude")
    wavenumber_value = _nonzero_real_quantity(wavenumber, "wavenumber")
    thickness_value = _positive_quantity(thickness, "thickness")
    return sp.simplify(
        coupling_value
        * amplitude_value**2
        * wavenumber_value**2
        * thickness_value
        / 2
    )


def frank_quadratic_capillary_map(
    frank_constant: Any,
    defect_strength: Any,
    outer_cutoff: Any,
    core_cutoff: Any,
    core_line_energy: Any,
    coupling: Any,
    amplitude: Any,
    wavenumber: Any,
    thickness: Any,
) -> CapillaryConstitutiveResult:
    """Compose declared Frank/core tension and quadratic loading exactly."""

    tension = frank_core_line_tension(
        frank_constant,
        defect_strength,
        outer_cutoff,
        core_cutoff,
        core_line_energy,
    )
    drive = quadratic_loading_area_drive(
        coupling,
        amplitude,
        wavenumber,
        thickness,
    )
    return CapillaryConstitutiveResult(
        line_tension=tension,
        area_drive=drive,
        critical_radius=capillary_critical_radius(tension, drive),
        barrier_height=capillary_barrier_height(tension, drive),
    )


def monomial_loading_dimension_ledger(
    amplitude_length_exponent: Any,
    *,
    amplitude_power: Any = 2,
    wavenumber_power: Any = 2,
) -> CapillaryDimensionLedger:
    r"""Return dimensions for ``p~g*A**m*k**n*l_m``.

    With ``[A]=L**alpha`` and ``[k]=L**-1``, area-drive closure requires
    ``[g]=E*L**(n-3-m*alpha)``.  The source specialization ``m=n=2``
    therefore has ``[g]=E*L**(-1-2*alpha)``.  Other exponent choices close
    after a corresponding coupling convention, so dimensions do not select
    the quadratic law or ``alpha``.
    """

    alpha = _real_quantity(
        amplitude_length_exponent,
        "amplitude_length_exponent",
    )
    amplitude_exponent = _real_quantity(amplitude_power, "amplitude_power")
    wavenumber_exponent = _real_quantity(
        wavenumber_power,
        "wavenumber_power",
    )
    coupling_length = sp.simplify(
        wavenumber_exponent - 3 - amplitude_exponent * alpha
    )
    names = (
        "line_tension",
        "area_drive",
        "critical_radius",
        "barrier_height",
        "frank_constant",
        "defect_strength",
        "outer_cutoff",
        "core_cutoff",
        "core_line_energy",
        "coupling",
        "amplitude",
        "wavenumber",
        "thickness",
        "bulk_bias",
    )
    matrix = sp.Matrix(
        [
            [1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1],
            [-1, -2, 1, 0, -1, 0, 1, 1, -1, coupling_length, alpha, -1, 1, -3],
        ]
    )
    return CapillaryDimensionLedger(
        base_dimensions=("E", "L"),
        quantity_names=names,
        dimension_matrix=sp.ImmutableMatrix(matrix),
        amplitude_length_exponent=alpha,
        amplitude_power=amplitude_exponent,
        wavenumber_power=wavenumber_exponent,
        coupling_length_exponent=coupling_length,
    )


def _coordinate_identifiability(
    nullspace: tuple[sp.ImmutableMatrix, ...],
    column_count: int,
) -> tuple[bool, ...]:
    return tuple(
        all(sp.simplify(vector[index]) == 0 for vector in nullspace)
        for index in range(column_count)
    )


def quadratic_capillary_identifiability_ledger() -> CapillaryIdentifiabilityLedger:
    """Return exact log-identifiability for effective quadratic inputs.

    Columns are ``(T,g,A,k,l_m)``.  Rows are critical radius and relative
    barrier height.  A barrier alone identifies no constituent.  Radius plus
    barrier identifies the effective line tension, but only the combined area
    drive—not ``g``, ``A``, ``k``, or ``l_m`` separately.  Decomposing ``T``
    into Frank/core inputs introduces further non-identifiability.
    """

    matrix = sp.Matrix(
        [
            [1, -1, -2, -2, -1],
            [2, -1, -2, -2, -1],
        ]
    )
    nullspace = tuple(sp.ImmutableMatrix(vector) for vector in matrix.nullspace())
    barrier = sp.Matrix([list(matrix.row(1))])
    barrier_nullspace = tuple(
        sp.ImmutableMatrix(vector) for vector in barrier.nullspace()
    )
    return CapillaryIdentifiabilityLedger(
        parameter_names=(
            "line_tension",
            "coupling",
            "amplitude",
            "wavenumber",
            "thickness",
        ),
        observable_names=("critical_radius", "barrier_height"),
        exponent_matrix=sp.ImmutableMatrix(matrix),
        rank=int(matrix.rank()),
        nullspace=nullspace,
        coordinate_identifiable=_coordinate_identifiability(nullspace, matrix.cols),
        barrier_only_matrix=sp.ImmutableMatrix(barrier),
        barrier_only_rank=int(barrier.rank()),
        barrier_only_nullspace=barrier_nullspace,
        barrier_only_coordinate_identifiable=_coordinate_identifiability(
            barrier_nullspace,
            barrier.cols,
        ),
    )


def equivalent_quadratic_loading_parameters(
    coupling: Any,
    amplitude: Any,
    wavenumber: Any,
    thickness: Any,
    *,
    amplitude_factor: Any = 1,
    wavenumber_factor: Any = 1,
    thickness_factor: Any = 1,
) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    """Return a positive rescaling family with unchanged quadratic drive."""

    coupling_value = _positive_quantity(coupling, "coupling")
    amplitude_value = _nonzero_real_quantity(amplitude, "amplitude")
    wavenumber_value = _nonzero_real_quantity(wavenumber, "wavenumber")
    thickness_value = _positive_quantity(thickness, "thickness")
    amplitude_scale = _positive_quantity(amplitude_factor, "amplitude_factor")
    wavenumber_scale = _positive_quantity(
        wavenumber_factor,
        "wavenumber_factor",
    )
    thickness_scale = _positive_quantity(thickness_factor, "thickness_factor")
    changed_coupling = sp.simplify(
        coupling_value
        / (amplitude_scale**2 * wavenumber_scale**2 * thickness_scale)
    )
    return (
        changed_coupling,
        sp.simplify(amplitude_scale * amplitude_value),
        sp.simplify(wavenumber_scale * wavenumber_value),
        sp.simplify(thickness_scale * thickness_value),
    )


def frank_core_log_sensitivities(
    frank_constant: Any,
    defect_strength: Any,
    outer_cutoff: Any,
    core_cutoff: Any,
    core_line_energy: Any,
) -> FrankCoreSensitivityLedger:
    """Return exact local log elasticities for the Frank/core components."""

    stiffness = _positive_quantity(frank_constant, "frank_constant")
    strength = _real_quantity(defect_strength, "defect_strength")
    outer = _positive_quantity(outer_cutoff, "outer_cutoff")
    core = _positive_quantity(core_cutoff, "core_cutoff")
    core_energy = _positive_quantity(core_line_energy, "core_line_energy")
    tension = frank_core_line_tension(
        stiffness,
        strength,
        outer,
        core,
        core_energy,
    )
    elastic = sp.pi * stiffness * strength**2 * sp.log(outer / core)
    elasticities = tuple(
        sp.simplify(value)
        for value in (
            elastic / tension,
            2 * elastic / tension,
            sp.pi * stiffness * strength**2 / tension,
            -sp.pi * stiffness * strength**2 / tension,
            core_energy / tension,
        )
    )
    return FrankCoreSensitivityLedger(
        parameter_names=(
            "frank_constant",
            "defect_strength",
            "outer_cutoff",
            "core_cutoff",
            "core_line_energy",
        ),
        line_tension=tension,
        elastic_line_tension=elastic,
        line_tension_log_elasticities=elasticities,
        critical_radius_log_elasticities=elasticities,
        barrier_height_log_elasticities=tuple(
            sp.simplify(2 * value) for value in elasticities
        ),
    )
