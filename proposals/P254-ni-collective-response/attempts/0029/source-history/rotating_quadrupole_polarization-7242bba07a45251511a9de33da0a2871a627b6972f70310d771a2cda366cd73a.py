"""Generic-observer TT ellipse algebra for a prescribed perpendicular rotation.

The phase and angular speed are separate inputs: ``phase`` locates the
prescribed tensor on its rotation orbit, while ``angular_speed**2`` supplies
the second-time-derivative scale.  This separation makes a fixed-phase scale
cancellation impossible to misread as independence of a physical-time signal
from angular speed.  All waveform quantities remain conditional on C-GW-001's
declared far-zone field premise; this module supplies no source dynamics,
gravity theory, propagation model, or detector observable.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .conditional_triaxial_radiation import (
    ConditionalScaledSTFWaveform,
)


@dataclass(frozen=True)
class PerpendicularRotationTTPolarization:
    """Conditional phase-parametrized TT readout and ellipse invariants.

    ``coefficient_matrix`` has conventional ``(plus, cross)`` coordinates as
    rows and ``(cos(2*(phase-azimuth)), sin(2*(phase-azimuth)))`` as columns.
    Its determinant decides temporal rank when the common waveform scale is
    nonzero. ``phase_gram_matrix = M.T*M`` is invariant under an orthogonal
    rotation of the transverse polarization frame.
    """

    transverse_eigenvalue: sp.Expr
    quadrupole_scale: sp.Expr
    angular_speed: sp.Expr
    phase: sp.Expr
    inclination: sp.Expr
    observer_azimuth: sp.Expr
    transverse_frame_angle: sp.Expr
    gravitational_coupling: sp.Expr
    distance: sp.Expr
    harmonic_phase: sp.Expr
    line_of_sight: sp.Matrix
    meridian_transverse: sp.Matrix
    azimuthal_transverse: sp.Matrix
    first_transverse: sp.Matrix
    second_transverse: sp.Matrix
    source_second_derivative: sp.Matrix
    waveform: ConditionalScaledSTFWaveform
    common_waveform_scale: sp.Expr
    coefficient_matrix: sp.Matrix
    normalized_coefficient_matrix: sp.Matrix
    phase_gram_matrix: sp.Matrix
    coefficient_determinant: sp.Expr


def perpendicular_axisymmetric_stf_second_derivative(
    transverse_eigenvalue: Any,
    angular_speed: Any,
    phase: Any,
    quadrupole_scale: Any = 1,
) -> sp.Matrix:
    """Return the exact second-time derivative at a declared orbit phase.

    The body moment has symmetry axis ``z`` and is prescribed to rotate about
    ``x``.  Differentiation is performed with respect to a dimensionless phase
    and multiplied by ``angular_speed**2``.  Thus angular speed zero gives the
    static limit without dividing by angular speed.
    """

    value = sp.sympify(transverse_eigenvalue)
    speed = sp.sympify(angular_speed)
    angle = sp.sympify(phase)
    scale = sp.sympify(quadrupole_scale)
    if sp.simplify(scale) == 0:
        raise ValueError("quadrupole_scale must be nonzero")
    cosine = sp.cos(2 * angle)
    sine = sp.sin(2 * angle)
    return 6 * scale * value * speed**2 * sp.Matrix(
        [
            [0, 0, 0],
            [0, -cosine, -sine],
            [0, -sine, cosine],
        ]
    )


def conditional_perpendicular_rotation_polarization(
    transverse_eigenvalue: Any,
    angular_speed: Any,
    phase: Any,
    inclination: Any,
    gravitational_coupling: Any,
    distance: Any,
    *,
    observer_azimuth: Any = 0,
    transverse_frame_angle: Any = 0,
    quadrupole_scale: Any = 1,
) -> PerpendicularRotationTTPolarization:
    """Return the exact generic-observer conditional TT polarization ellipse.

    The rotation axis is ``x``. ``inclination`` is measured from that axis and
    ``observer_azimuth`` locates the sightline in the ``y-z`` plane.  The
    natural first transverse vector is the projected positive rotation axis;
    ``transverse_frame_angle`` then rotates the transverse pair while
    preserving its orientation.

    The returned coefficient matrix is the exact contraction in that oriented
    frame. Its two singular directions determine the polarization ellipse.
    Full projector composition is intentionally kept in the independent
    verifier rather than recomputed on every API call.
    """

    value = sp.sympify(transverse_eigenvalue)
    speed = sp.sympify(angular_speed)
    orbit_phase = sp.sympify(phase)
    angle = sp.sympify(inclination)
    azimuth = sp.sympify(observer_azimuth)
    frame_angle = sp.sympify(transverse_frame_angle)
    coupling = sp.sympify(gravitational_coupling)
    radius = sp.sympify(distance)
    scale = sp.sympify(quadrupole_scale)
    if sp.simplify(radius) == 0:
        raise ValueError("distance must be nonzero")
    if sp.simplify(scale) == 0:
        raise ValueError("quadrupole_scale must be nonzero")

    cosine = sp.cos(angle)
    sine = sp.sin(angle)
    line_of_sight = sp.Matrix(
        [cosine, sine * sp.cos(azimuth), sine * sp.sin(azimuth)]
    )
    meridian = sp.Matrix(
        [sine, -cosine * sp.cos(azimuth), -cosine * sp.sin(azimuth)]
    )
    azimuthal = sp.Matrix([0, sp.sin(azimuth), -sp.cos(azimuth)])
    first = (
        sp.cos(frame_angle) * meridian
        + sp.sin(frame_angle) * azimuthal
    )
    second = (
        -sp.sin(frame_angle) * meridian
        + sp.cos(frame_angle) * azimuthal
    )

    coordinate_rotation = sp.Matrix(
        [
            [sp.cos(2 * frame_angle), sp.sin(2 * frame_angle)],
            [-sp.sin(2 * frame_angle), sp.cos(2 * frame_angle)],
        ]
    )
    inclination_cosine = sp.cos(angle)
    semimajor_factor = (1 + inclination_cosine**2) / 2
    common_scale = 12 * coupling * value * speed**2 / radius
    source_scale = 6 * scale * value * speed**2
    coefficient_matrix = coordinate_rotation * sp.diag(
        -common_scale * semimajor_factor,
        -common_scale * inclination_cosine,
    )
    normalized_coefficients = sp.sqrt(2) * coefficient_matrix
    source_coefficient_matrix = coordinate_rotation * sp.diag(
        -source_scale * semimajor_factor,
        -source_scale * inclination_cosine,
    )
    harmonic_phase = 2 * (orbit_phase - azimuth)
    phase_vector = sp.Matrix(
        [sp.cos(harmonic_phase), sp.sin(harmonic_phase)]
    )
    actual_conventional = coefficient_matrix * phase_vector
    actual_normalized = sp.sqrt(2) * actual_conventional
    actual_source_coordinates = source_coefficient_matrix * phase_vector
    conventional_plus_tensor = first * first.T - second * second.T
    conventional_cross_tensor = first * second.T + second * first.T
    actual_derivative = perpendicular_axisymmetric_stf_second_derivative(
        value,
        speed,
        orbit_phase,
        scale,
    )
    actual_projected = (
        actual_source_coordinates[0] * conventional_plus_tensor
        + actual_source_coordinates[1] * conventional_cross_tensor
    )
    actual_tensor = (
        actual_conventional[0] * conventional_plus_tensor
        + actual_conventional[1] * conventional_cross_tensor
    )
    waveform_prefactor = 2 * coupling / (scale * radius)
    waveform = ConditionalScaledSTFWaveform(
        quadrupole_scale=scale,
        gravitational_coupling=coupling,
        distance=radius,
        waveform_prefactor=waveform_prefactor,
        quadrupole_second_derivative=actual_derivative,
        projected_quadrupole_second_derivative=actual_projected,
        waveform_tensor=actual_tensor,
        normalized_plus=actual_normalized[0],
        normalized_cross=actual_normalized[1],
        conventional_plus=actual_conventional[0],
        conventional_cross=actual_conventional[1],
    )
    phase_gram = sp.diag(
        common_scale**2 * semimajor_factor**2,
        common_scale**2 * inclination_cosine**2,
    )
    determinant = common_scale**2 * semimajor_factor * inclination_cosine
    return PerpendicularRotationTTPolarization(
        transverse_eigenvalue=value,
        quadrupole_scale=scale,
        angular_speed=speed,
        phase=orbit_phase,
        inclination=angle,
        observer_azimuth=azimuth,
        transverse_frame_angle=frame_angle,
        gravitational_coupling=coupling,
        distance=radius,
        harmonic_phase=harmonic_phase,
        line_of_sight=line_of_sight,
        meridian_transverse=meridian,
        azimuthal_transverse=azimuthal,
        first_transverse=first,
        second_transverse=second,
        source_second_derivative=waveform.quadrupole_second_derivative,
        waveform=waveform,
        common_waveform_scale=common_scale,
        coefficient_matrix=coefficient_matrix,
        normalized_coefficient_matrix=normalized_coefficients,
        phase_gram_matrix=phase_gram,
        coefficient_determinant=determinant,
    )
