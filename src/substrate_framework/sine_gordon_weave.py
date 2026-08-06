"""Collective topological sine-Gordon weave and its Einstein coupling.

The microscopic links use the accepted dimensional sine-Gordon scales and
the exact parity-related winding sectors ``Q=+1`` and ``Q=-1``.  Four shared
links per cell form a tetrahedral coordination frame.  Coarse cells have
uniformly distributed orientations, making the line ensemble statistically
isotropic.  An equiprobable topological doublet on every horizon crossing
then fixes the area entropy density passed to local horizon equilibrium.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .dimensional_sine_gordon import (
    DimensionalSineGordonCoefficients,
    DimensionalSineGordonScales,
    dimensional_sine_gordon_scales,
)
from .local_horizon_gravity import (
    LocalHorizonEinsteinLedger,
    local_horizon_einstein_ledger,
)


@dataclass(frozen=True)
class TetrahedralWeaveGeometry:
    """Exact cell and isotropic crossing geometry of the link weave."""

    cell_length: sp.Expr
    unit_directions: tuple[sp.ImmutableMatrix, ...]
    direction_sum: sp.ImmutableMatrix
    direction_second_moment: sp.ImmutableMatrix
    coordination: int
    shared_line_length_per_cell: sp.Expr
    line_length_density: sp.Expr
    rotational_mean_absolute_projection: sp.Expr
    crossing_density: sp.Expr


@dataclass(frozen=True)
class CollectiveWeaveMetric:
    """Lorentzian metric and measure reconstructed from one link coframe."""

    coframe: sp.ImmutableMatrix
    clock_covector: sp.ImmutableMatrix
    link_covectors: tuple[sp.ImmutableMatrix, ...]
    link_second_moment: sp.ImmutableMatrix
    covariant_metric: sp.ImmutableMatrix
    contravariant_metric: sp.ImmutableMatrix
    volume_density: sp.Expr


@dataclass(frozen=True)
class SineGordonWeaveGravity:
    """Microscopic weave data and the resulting nonlinear GR coefficient."""

    coefficients: DimensionalSineGordonCoefficients
    scales: DimensionalSineGordonScales
    geometry: TetrahedralWeaveGeometry
    topological_sector_count: int
    entropy_per_crossing: sp.Expr
    entropy_area_density: sp.Expr
    equilibrium_metric: CollectiveWeaveMetric
    horizon_equilibrium: LocalHorizonEinsteinLedger


def tetrahedral_weave_geometry(cell_length: sp.Expr) -> TetrahedralWeaveGeometry:
    """Return the exact coordination and plane-crossing density of one cell.

    Four length-``a`` segments are shared by neighboring cells, hence the
    line length density is ``(4*a/2)/a**3=2/a**2``.  Haar averaging a line
    direction against a fixed plane normal gives ``<|cos(theta)|>=1/2``, so
    the statistically isotropic crossing density is ``1/a**2``.
    """

    length = sp.sympify(cell_length)
    if length.has(sp.Float):
        raise ValueError("cell_length must be exact rather than floating")
    if length.is_real is not True or length.is_positive is not True:
        raise ValueError("cell_length must be provably positive and real")
    directions = tuple(
        sp.ImmutableMatrix(vector) / sp.sqrt(3)
        for vector in (
            (1, 1, 1),
            (1, -1, -1),
            (-1, 1, -1),
            (-1, -1, 1),
        )
    )
    direction_sum = sp.ImmutableMatrix(sum(directions, sp.zeros(3, 1)))
    second_moment = sp.ImmutableMatrix(
        sum((direction * direction.T for direction in directions), sp.zeros(3))
        / len(directions)
    ).applyfunc(sp.simplify)
    if direction_sum != sp.zeros(3, 1):
        raise AssertionError("tetrahedral directions do not close")
    if second_moment != sp.eye(3) / 3:
        raise AssertionError("tetrahedral direction covariance is not isotropic")
    shared_length = 2 * length
    line_density = sp.simplify(shared_length / length**3)
    mean_projection = sp.Rational(1, 2)
    return TetrahedralWeaveGeometry(
        cell_length=length,
        unit_directions=directions,
        direction_sum=direction_sum,
        direction_second_moment=second_moment,
        coordination=4,
        shared_line_length_per_cell=shared_length,
        line_length_density=line_density,
        rotational_mean_absolute_projection=mean_projection,
        crossing_density=sp.simplify(line_density * mean_projection),
    )


def collective_weave_metric(coframe: sp.MatrixBase) -> CollectiveWeaveMetric:
    r"""Reconstruct the collective metric from clock and link deformations.

    Rows of the invertible coframe ``e**A_mu`` are one phase-clock covector and
    three spatial link-frame covectors.  Contracting the latter with the four
    tetrahedral directions gives observable link covectors ``l_alpha``.  Their
    exact second moment obeys

    ``<l_alpha l_alpha> = e_spatial.T*e_spatial/3``.

    Therefore the link ensemble itself reconstructs
    ``g=clock.T*clock-3*<l l>=e.T*diag(1,-1,-1,-1)*e`` and its invariant volume
    density.  This makes the 3+1 metric and measure collective weave degrees of
    freedom rather than inputs to the sine-Gordon channel stress.
    """

    frame = sp.Matrix(coframe)
    if frame.shape != (4, 4):
        raise ValueError("coframe must be 4 by 4")
    if any(entry.has(sp.Float) for entry in frame):
        raise ValueError("coframe must be exact rather than floating")
    determinant = sp.simplify(frame.det())
    if determinant.is_zero is not False:
        raise ValueError("coframe determinant must be provably nonzero")
    immutable_frame = sp.ImmutableMatrix(frame)
    clock = sp.ImmutableMatrix(frame[0, :])
    spatial_frame = frame[1:, :]
    geometry = tetrahedral_weave_geometry(sp.Integer(1))
    links = tuple(
        sp.ImmutableMatrix(direction.T * spatial_frame)
        for direction in geometry.unit_directions
    )
    link_second_moment = sp.ImmutableMatrix(
        sum((link.T * link for link in links), sp.zeros(4)) / len(links)
    ).applyfunc(sp.simplify)
    internal_metric = sp.diag(1, -1, -1, -1)
    metric = sp.ImmutableMatrix(frame.T * internal_metric * frame).applyfunc(sp.simplify)
    reconstructed = sp.ImmutableMatrix(
        clock.T * clock - 3 * link_second_moment
    ).applyfunc(sp.simplify)
    if reconstructed != metric:
        raise AssertionError("tetrahedral links did not reconstruct the metric")
    inverse = sp.ImmutableMatrix(metric.inv()).applyfunc(sp.simplify)
    volume_density = sp.simplify(sp.sqrt(-metric.det()))
    return CollectiveWeaveMetric(
        coframe=immutable_frame,
        clock_covector=clock,
        link_covectors=links,
        link_second_moment=link_second_moment,
        covariant_metric=metric,
        contravariant_metric=inverse,
        volume_density=volume_density,
    )


def sine_gordon_weave_gravity(
    coefficients: DimensionalSineGordonCoefficients,
) -> SineGordonWeaveGravity:
    r"""Construct the weave and derive its total Newton coupling.

    The cell length is the unique sine-Gordon length ``ell`` and the Unruh
    action is the same model action ``J=sqrt(lambda*T)``.  No independent
    cutoff, action quantum, Newton constant, or additive inverse-coupling
    baseline enters the API.  The two exact parity-related winding sectors
    give ``log(2)`` entropy per crossing, so

    ``eta_A=log(2)/ell**2``,
    ``8*pi*G/c**4=2*pi/(mu*log(2))``, and
    ``G=T**2/(4*lambda**2*mu*log(2))``.

    The flat equilibrium branch fixes the theorem's cosmological integration
    constant to zero.
    """

    scales = dimensional_sine_gordon_scales(coefficients)
    geometry = tetrahedral_weave_geometry(scales.length)
    sector_count = 2
    entropy_per_crossing = sp.log(sector_count)
    entropy_density = sp.simplify(
        geometry.crossing_density * entropy_per_crossing
    )
    horizon = local_horizon_einstein_ledger(
        entropy_density,
        scales.action,
        scales.signal_speed,
        cosmological_constant=0,
    )
    expected_coupling = sp.simplify(2 * sp.pi / (coefficients.onsite * sp.log(2)))
    if sp.simplify(horizon.einstein_stress_coupling - expected_coupling) != 0:
        raise AssertionError("sine-Gordon weave coupling did not close upstream")
    equilibrium_metric = collective_weave_metric(sp.eye(4))
    return SineGordonWeaveGravity(
        coefficients=coefficients,
        scales=scales,
        geometry=geometry,
        topological_sector_count=sector_count,
        entropy_per_crossing=entropy_per_crossing,
        entropy_area_density=entropy_density,
        equilibrium_metric=equilibrium_metric,
        horizon_equilibrium=horizon,
    )
