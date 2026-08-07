"""Exact six-family sine-Gordon wall frame and its collective scale ledger.

This module declares a local 3+1 construction rather than identifying a
scalar spherical harmonic with a graviton.  Six antipodal sine-Gordon kink
wall families use the unoriented axes of an icosahedron.  One wall separation
per accepted sine-Gordon length defines each dimensionless directional strain
``q_alpha``.  The six strains reconstruct a symmetric spatial tensor exactly.

The construction fixes its quadratic normalization from the accepted
dimensional sine-Gordon coefficients and the exact kink integral.  It does
not by itself prove the collective force-balance constraints or their
gauge-invariant relativistic completion; those are audited separately in
``emergent_spin2``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import sympy as sp

from .dimensional_sine_gordon import (
    DimensionalSineGordonCoefficients,
    dimensional_sine_gordon_coefficients,
    dimensional_sine_gordon_scales,
)
from .exact_symbolic import exact_real as _exact_real


_SYMMETRIC_COMPONENT_ORDER = ("xx", "yy", "zz", "xy", "xz", "yz")
_SYMMETRIC_INDEX_PAIRS = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))


@dataclass(frozen=True)
class IcosahedralWallFrameLedger:
    """Exact tight-frame and reconstruction data for six wall families."""

    component_order: tuple[str, ...]
    directions: tuple[sp.ImmutableMatrix, ...]
    directional_projection: sp.ImmutableMatrix
    reconstruction: sp.ImmutableMatrix
    projection_rank: int
    direction_second_moment: sp.ImmutableMatrix
    trace_quadratic_weight: sp.Expr
    traceless_quadratic_weight: sp.Expr
    tight_frame_identity_residual: sp.Expr


@dataclass(frozen=True)
class SineGordonWallNetworkScales:
    """Microscopic-to-collective normalization for the declared wall cell.

    ``wall_inertia_density`` and ``wall_modulus`` multiply one directional
    separation mode in the conventional quadratic form
    ``(rho*q_dot**2-C*|grad(q)|**2)/2``.  ``spin2_spacetime_normalization`` is
    the coefficient ``B`` after the exact traceless-frame contraction, so the
    TT density is ``-B*partial_a(h_ij)*partial^a(h_ij)/2``.
    """

    coefficients: DimensionalSineGordonCoefficients
    signal_speed: sp.Expr
    profile_length: sp.Expr
    energy_scale: sp.Expr
    kink_gradient_integral: sp.Expr
    wall_tension: sp.Expr
    walls_per_length: sp.Expr
    wall_inertia_density: sp.Expr
    wall_modulus: sp.Expr
    traceless_frame_weight: sp.Expr
    spin2_time_inertia: sp.Expr
    spin2_spacetime_normalization: sp.Expr


@dataclass(frozen=True)
class LocalSineGordonWallNetwork:
    """Local 3+1 six-field density with explicit transverse derivatives."""

    fields: tuple[sp.Expr, ...]
    spatial_coordinates: tuple[sp.Symbol, ...]
    time: sp.Symbol
    lagrangian_density: sp.Expr
    field_equation_residuals: sp.ImmutableMatrix
    spatial_gradient_order: tuple[tuple[int, int], ...]
    spatial_gradient_hessian: sp.ImmutableMatrix


def icosahedral_wall_directions() -> tuple[sp.ImmutableMatrix, ...]:
    """Return one unit representative of each antipodal icosahedral axis."""

    golden_ratio = (1 + sp.sqrt(5)) / 2
    normalization = sp.sqrt(1 + golden_ratio**2)
    raw = (
        (0, 1, golden_ratio),
        (0, -1, golden_ratio),
        (1, golden_ratio, 0),
        (-1, golden_ratio, 0),
        (golden_ratio, 0, 1),
        (golden_ratio, 0, -1),
    )
    return tuple(
        sp.ImmutableMatrix([sp.simplify(value / normalization) for value in row])
        for row in raw
    )


def _projection_row(direction: sp.ImmutableMatrix) -> list[sp.Expr]:
    return [
        (2 if row != column else 1) * direction[row] * direction[column]
        for row, column in _SYMMETRIC_INDEX_PAIRS
    ]


def _validated_wall_network_inputs(
    fields: Iterable[Any],
    spatial_coordinates: Iterable[sp.Symbol],
    time: sp.Symbol,
    coefficients: DimensionalSineGordonCoefficients,
) -> tuple[
    tuple[sp.Expr, ...],
    tuple[sp.Symbol, ...],
    DimensionalSineGordonCoefficients,
]:
    field_tuple = tuple(sp.sympify(field) for field in fields)
    if len(field_tuple) != 6:
        raise ValueError(
            f"fields must contain exactly six entries; got {len(field_tuple)}"
        )
    coordinates = tuple(spatial_coordinates)
    if len(coordinates) != 3 or not all(
        isinstance(coordinate, sp.Symbol) for coordinate in coordinates
    ):
        raise ValueError(
            "spatial_coordinates must contain exactly three SymPy symbols; "
            f"got {coordinates!r}"
        )
    if len(set(coordinates)) != 3:
        raise ValueError(
            f"spatial_coordinates must be distinct; got {coordinates!r}"
        )
    if not isinstance(time, sp.Symbol) or time in coordinates:
        raise ValueError(
            "time must be a SymPy symbol distinct from spatial_coordinates; "
            f"got {time!r}"
        )
    if not isinstance(coefficients, DimensionalSineGordonCoefficients):
        raise TypeError("coefficients must be a DimensionalSineGordonCoefficients record")
    validated = dimensional_sine_gordon_coefficients(
        coefficients.inertia,
        coefficients.gradient,
        coefficients.onsite,
    )
    return field_tuple, coordinates, validated


def _wall_field_density(
    field: sp.Expr,
    spatial_coordinates: tuple[sp.Symbol, ...],
    time: sp.Symbol,
    coefficients: DimensionalSineGordonCoefficients,
    cell_area: sp.Expr,
) -> sp.Expr:
    return (
        coefficients.inertia * sp.diff(field, time) ** 2 / 2
        - coefficients.gradient
        * sum(sp.diff(field, coordinate) ** 2 for coordinate in spatial_coordinates)
        / 2
        - coefficients.onsite * (1 - sp.cos(field))
    ) / cell_area


def _wall_field_residual(
    field: sp.Expr,
    spatial_coordinates: tuple[sp.Symbol, ...],
    time: sp.Symbol,
    coefficients: DimensionalSineGordonCoefficients,
    cell_area: sp.Expr,
) -> sp.Expr:
    return sp.simplify(
        (
            coefficients.inertia * sp.diff(field, time, 2)
            - coefficients.gradient
            * sum(sp.diff(field, coordinate, 2) for coordinate in spatial_coordinates)
            + coefficients.onsite * sp.sin(field)
        )
        / cell_area
    )


def _spatial_gradient_hessian(
    lagrangian: sp.Expr,
    gradient_variables: list[sp.Expr],
    expected_coefficient: sp.Expr,
) -> sp.ImmutableMatrix:
    hessian = sp.ImmutableMatrix(sp.hessian(sp.expand(lagrangian), gradient_variables))
    if sp.simplify(hessian - expected_coefficient * sp.eye(18)) != sp.zeros(18):
        raise AssertionError("local wall density lost an isotropic spatial derivative")
    return hessian


def local_sine_gordon_wall_network(
    fields: Iterable[Any],
    spatial_coordinates: Iterable[sp.Symbol],
    time: sp.Symbol,
    coefficients: DimensionalSineGordonCoefficients,
) -> LocalSineGordonWallNetwork:
    r"""Return the declared local six-channel 3+1 sine-Gordon density.

    Each wall family is a full 3+1 field rather than an uncoupled 1+1 line:

    ``sum_alpha [lambda*u_t^2-T*sum_i u_i^2]/(2*ell^2)
    -mu*(1-cos(u))/ell^2``.

    The accepted length ``ell`` supplies the cell cross-section.  Registry
    admissibility constraints act on the collective wall separations and are
    intentionally separate from these microscopic Euler--Lagrange equations.
    """

    field_tuple, coordinates, validated = _validated_wall_network_inputs(
        fields,
        spatial_coordinates,
        time,
        coefficients,
    )
    ell = dimensional_sine_gordon_scales(validated).length
    cell_area = ell**2
    lagrangian = sp.Integer(0)
    residuals: list[sp.Expr] = []
    gradient_variables: list[sp.Expr] = []
    gradient_order: list[tuple[int, int]] = []
    for field_index, field in enumerate(field_tuple):
        spatial_derivatives = tuple(
            sp.diff(field, coordinate) for coordinate in coordinates
        )
        lagrangian += _wall_field_density(
            field, coordinates, time, validated, cell_area
        )
        residuals.append(
            _wall_field_residual(field, coordinates, time, validated, cell_area)
        )
        for coordinate_index, derivative in enumerate(spatial_derivatives):
            gradient_variables.append(derivative)
            gradient_order.append((field_index, coordinate_index))
    gradient_hessian = _spatial_gradient_hessian(
        lagrangian,
        gradient_variables,
        -validated.gradient / cell_area,
    )
    return LocalSineGordonWallNetwork(
        fields=field_tuple,
        spatial_coordinates=coordinates,
        time=time,
        lagrangian_density=sp.simplify(lagrangian),
        field_equation_residuals=sp.ImmutableMatrix(residuals),
        spatial_gradient_order=tuple(gradient_order),
        spatial_gradient_hessian=gradient_hessian,
    )


def icosahedral_wall_frame_ledger() -> IcosahedralWallFrameLedger:
    r"""Return the exact six-channel tensor reconstruction and tight-frame law.

    For a symmetric tensor ``h`` and ``q_alpha=n_alpha^i n_alpha^j h_ij``,

    ``sum(q_alpha**2) = 2*tr(h)**2/5 + 4*tr(h**2)/5``.

    The trace-free sector therefore has weight ``4/5`` with respect to the
    Frobenius norm.  Both the rank and identity are derived, not supplied.
    """

    directions = icosahedral_wall_directions()
    projection = sp.ImmutableMatrix([_projection_row(n) for n in directions])
    reconstruction = sp.ImmutableMatrix(projection.inv())
    second_moment = sp.ImmutableMatrix(
        sp.simplify(sum((n * n.T for n in directions), sp.zeros(3)))
    )

    hxx, hyy, hzz, hxy, hxz, hyz = sp.symbols(
        "h_xx h_yy h_zz h_xy h_xz h_yz", real=True
    )
    components = sp.ImmutableMatrix([hxx, hyy, hzz, hxy, hxz, hyz])
    tensor = sp.ImmutableMatrix(
        [[hxx, hxy, hxz], [hxy, hyy, hyz], [hxz, hyz, hzz]]
    )
    channel_norm = sp.expand((projection * components).dot(projection * components))
    expected = sp.Rational(2, 5) * sp.trace(tensor) ** 2 + sp.Rational(
        4, 5
    ) * sp.trace(tensor * tensor)
    residual = sp.simplify(channel_norm - expected)
    if projection.rank() != 6 or residual != 0:
        raise AssertionError("icosahedral wall axes do not form the required frame")

    return IcosahedralWallFrameLedger(
        component_order=_SYMMETRIC_COMPONENT_ORDER,
        directions=directions,
        directional_projection=projection,
        reconstruction=reconstruction,
        projection_rank=projection.rank(),
        direction_second_moment=second_moment,
        trace_quadratic_weight=sp.Rational(2, 5),
        traceless_quadratic_weight=sp.Rational(4, 5),
        tight_frame_identity_residual=residual,
    )


def _symmetric_tensor(value: Any) -> sp.ImmutableMatrix:
    try:
        tensor = sp.ImmutableMatrix(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"spatial_tensor must be matrix-like; got {value!r}") from error
    if tensor.shape != (3, 3):
        raise ValueError(
            f"spatial_tensor must have shape (3, 3); got {tensor.shape}"
        )
    if tensor != tensor.T:
        raise ValueError(f"spatial_tensor must be symmetric; got {tensor!r}")
    entries: list[sp.Expr] = []
    for row in range(3):
        for column in range(3):
            entry = tensor[row, column]
            try:
                entries.append(
                    _exact_real(entry, f"spatial_tensor[{row},{column}]")
                )
            except ValueError as error:
                raise ValueError(f"{error}; got {entry!r}") from error
    return sp.ImmutableMatrix(
        3,
        3,
        entries,
    )


def _symmetric_components(tensor: sp.ImmutableMatrix) -> sp.ImmutableMatrix:
    return sp.ImmutableMatrix([tensor[row, column] for row, column in _SYMMETRIC_INDEX_PAIRS])


def _tensor_from_symmetric_components(components: Iterable[sp.Expr]) -> sp.ImmutableMatrix:
    tensor = sp.zeros(3)
    for value, (row, column) in zip(components, _SYMMETRIC_INDEX_PAIRS, strict=True):
        tensor[row, column] = value
        tensor[column, row] = value
    return sp.ImmutableMatrix(tensor)


def wall_channels_from_spatial_tensor(spatial_tensor: Any) -> sp.ImmutableMatrix:
    """Project an exact symmetric spatial tensor onto six wall separations."""

    tensor = _symmetric_tensor(spatial_tensor)
    return sp.ImmutableMatrix(
        sp.simplify(
            icosahedral_wall_frame_ledger().directional_projection
            * _symmetric_components(tensor)
        )
    )


def spatial_tensor_from_wall_channels(channels: Iterable[Any]) -> sp.ImmutableMatrix:
    """Reconstruct the unique exact symmetric tensor from six wall channels."""

    raw_values = tuple(channels)
    if len(raw_values) != 6:
        raise ValueError(
            f"channels must contain exactly six entries; got {len(raw_values)}"
        )
    values: list[sp.Expr] = []
    for index, value in enumerate(raw_values):
        try:
            values.append(_exact_real(value, f"channels[{index}]"))
        except ValueError as error:
            raise ValueError(f"{error}; got {value!r}") from error
    components = sp.simplify(
        icosahedral_wall_frame_ledger().reconstruction * sp.ImmutableMatrix(values)
    )
    return _tensor_from_symmetric_components(components)


def sine_gordon_wall_network_scales(
    coefficients: DimensionalSineGordonCoefficients,
) -> SineGordonWallNetworkScales:
    r"""Derive the wall and collective tensor coefficients.

    The normalized kink ``u=4*atan(exp(s/ell))`` obeys
    ``integral u_s**2 ds=8/ell``.  Lifting the accepted 1+1 density to one
    wall-cell cross-section ``ell**2`` and placing one wall separation per
    length ``ell`` gives ``rho_wall=8*lambda/ell**2`` and
    ``C_wall=8*T/ell**2=8*mu``.  The exact frame then fixes
    ``B=(4/5)*C_wall=32*mu/5``.  The separate gauge-completion ledger owns any
    gravitational interpretation of ``B``.
    """

    if not isinstance(coefficients, DimensionalSineGordonCoefficients):
        raise TypeError("coefficients must be a DimensionalSineGordonCoefficients record")
    validated = dimensional_sine_gordon_coefficients(
        coefficients.inertia,
        coefficients.gradient,
        coefficients.onsite,
    )
    scales = dimensional_sine_gordon_scales(validated)
    ell = scales.length
    kink_integral = sp.simplify(8 / ell)
    wall_tension = sp.simplify(8 * scales.energy / ell**2)
    wall_density = sp.simplify(1 / ell)
    wall_inertia = sp.simplify(8 * validated.inertia / ell**2)
    wall_modulus = sp.simplify(8 * validated.gradient / ell**2)
    frame_weight = icosahedral_wall_frame_ledger().traceless_quadratic_weight
    time_inertia = sp.simplify(frame_weight * wall_inertia)
    spacetime_normalization = sp.simplify(frame_weight * wall_modulus)
    if sp.simplify(time_inertia * scales.signal_speed**2 - spacetime_normalization) != 0:
        raise AssertionError("wall kinetic and gradient scales do not share signal speed")
    if sp.simplify(wall_tension * ell - wall_modulus) != 0:
        raise AssertionError(
            "wall tension and fractional displacement do not reproduce the modulus"
        )

    return SineGordonWallNetworkScales(
        coefficients=validated,
        signal_speed=scales.signal_speed,
        profile_length=ell,
        energy_scale=scales.energy,
        kink_gradient_integral=kink_integral,
        wall_tension=wall_tension,
        walls_per_length=wall_density,
        wall_inertia_density=wall_inertia,
        wall_modulus=wall_modulus,
        traceless_frame_weight=frame_weight,
        spin2_time_inertia=time_inertia,
        spin2_spacetime_normalization=spacetime_normalization,
    )
