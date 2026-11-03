"""Exact angular reduction for a declared transverse-traceless wave model.

This module proves projector and normalization algebra.  It deliberately does
not supply a gravitational field equation, a source-to-waveform map, or an
energy-flux law; callers must pass those conditional prefactors explicitly.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .conserved_moments import symmetric_trace_free


def _symmetric_three(matrix: Any, name: str) -> sp.Matrix:
    value = sp.Matrix(matrix)
    if value.shape != (3, 3):
        raise ValueError(f"{name} must be 3 by 3")
    if sp.simplify(value - value.T) != sp.zeros(3):
        raise ValueError(f"{name} must be symmetric")
    return value


def _column_three(vector: Any, name: str) -> sp.Matrix:
    value = sp.Matrix(vector)
    if value.shape not in ((3, 1), (1, 3)):
        raise ValueError(f"{name} must have three components")
    return value if value.shape == (3, 1) else value.T


def transverse_projector(direction: Any) -> sp.Matrix:
    """Return ``I-n*n.T/(n.n)`` for a nonzero three-direction."""

    vector = _column_three(direction, "direction")
    norm_squared = sp.simplify(vector.dot(vector))
    if norm_squared == 0:
        raise ValueError("direction must be nonzero")
    return sp.simplify(sp.eye(3) - vector * vector.T / norm_squared)


def tt_project_symmetric(tensor: Any, direction: Any) -> sp.Matrix:
    """Project a symmetric spatial tensor transverse and trace-free.

    With ``P`` the rank-two transverse projector, this returns
    ``P*S*P - P*trace(P*S*P)/2``.  The factor one-half is the trace removal in
    the two-dimensional plane orthogonal to the line of sight.
    """

    source = _symmetric_three(tensor, "tensor")
    projector = transverse_projector(direction)
    transverse = sp.simplify(projector * source * projector)
    return sp.simplify(transverse - projector * sp.trace(transverse) / 2)


def frobenius_norm_squared(tensor: Any) -> sp.Expr:
    """Return the exact Cartesian contraction ``S_ij*S_ij``."""

    value = sp.Matrix(tensor)
    return sp.simplify(
        sum(
            (value[i, j] ** 2 for i in range(value.rows) for j in range(value.cols)),
            sp.Integer(0),
        )
    )


def frobenius_inner_product(left: Any, right: Any) -> sp.Expr:
    """Return the exact Cartesian contraction ``A_ij*B_ij``."""

    left_matrix = sp.Matrix(left)
    right_matrix = sp.Matrix(right)
    if left_matrix.shape != right_matrix.shape:
        raise ValueError("tensor shapes must match")
    return sp.simplify(
        sum(
            (
                left_matrix[i, j] * right_matrix[i, j]
                for i in range(left_matrix.rows)
                for j in range(left_matrix.cols)
            ),
            sp.Integer(0),
        )
    )


def _orthonormal_symmetric_basis() -> tuple[sp.Matrix, ...]:
    diagonal = tuple(
        sp.Matrix(3, 3, lambda i, j, axis=axis: int(i == axis and j == axis))
        for axis in range(3)
    )
    off_diagonal = []
    for first, second in ((0, 1), (0, 2), (1, 2)):
        tensor = sp.zeros(3)
        tensor[first, second] = 1 / sp.sqrt(2)
        tensor[second, first] = 1 / sp.sqrt(2)
        off_diagonal.append(tensor)
    return diagonal + tuple(off_diagonal)


def tt_operator_matrix(direction: Any) -> sp.Matrix:
    """Return the TT operator on an orthonormal six-tensor basis."""

    basis = _orthonormal_symmetric_basis()
    return sp.simplify(
        sp.Matrix(
            6,
            6,
            lambda row, column: frobenius_inner_product(
                basis[row], tt_project_symmetric(basis[column], direction)
            ),
        )
    )


@dataclass(frozen=True)
class TTPolarizationBasis:
    """An oriented transverse frame and normalized real TT tensor basis."""

    direction: sp.Matrix
    first_transverse: sp.Matrix
    second_transverse: sp.Matrix
    plus: sp.Matrix
    cross: sp.Matrix


@dataclass(frozen=True)
class AxisymmetricSTFReadout:
    """Natural-frame TT readout of an arbitrary-axis axisymmetric STF tensor.

    ``axial_amplitude`` is the coefficient ``alpha`` in
    ``S=scale*alpha*(e*e.T-I/3)``; it is not the eigenvalue along ``e``.
    The axial eigenvalue is ``2*scale*alpha/3``. Basis coordinates use the
    normalized tensors of :class:`TTPolarizationBasis`; conventional matrix
    readouts are smaller by ``sqrt(2)``.
    """

    axial_amplitude: sp.Expr
    quadrupole_scale: sp.Expr
    symmetry_axis: sp.Matrix
    direction: sp.Matrix
    inclination_cosine: sp.Expr
    inclination_sine_squared: sp.Expr
    first_transverse: sp.Matrix
    second_transverse: sp.Matrix
    source_tensor: sp.Matrix
    projected_tensor: sp.Matrix
    normalized_plus_coordinate: sp.Expr
    normalized_cross_coordinate: sp.Expr
    conventional_plus_readout: sp.Expr
    conventional_cross_readout: sp.Expr


def tt_polarization_basis(
    direction: Any,
    reference: Any | None = None,
) -> TTPolarizationBasis:
    """Construct a normalized plus/cross basis for a nonzero direction.

    The reference must not be parallel to the direction.  If it is omitted, a
    usable Cartesian axis is selected piecewise.  This is a deterministic
    all-direction construction, not a claim of a globally continuous frame on
    the sphere.
    """

    vector = _column_three(direction, "direction")
    norm_squared = sp.simplify(vector.dot(vector))
    if norm_squared == 0:
        raise ValueError("direction must be nonzero")
    unit_direction = sp.simplify(vector / sp.sqrt(norm_squared))

    if reference is None:
        selected_reference = None
        for axis in (sp.Matrix([1, 0, 0]), sp.Matrix([0, 1, 0]), sp.Matrix([0, 0, 1])):
            transverse = sp.simplify(
                axis - unit_direction * unit_direction.dot(axis)
            )
            if sp.simplify(transverse.dot(transverse)) != 0:
                selected_reference = axis
                break
        if selected_reference is None:
            raise ValueError("could not select a transverse reference")
    else:
        selected_reference = _column_three(reference, "reference")

    first_raw = sp.simplify(
        selected_reference
        - unit_direction * unit_direction.dot(selected_reference)
    )
    first_norm_squared = sp.simplify(first_raw.dot(first_raw))
    if first_norm_squared == 0:
        raise ValueError("reference must not be parallel to direction")
    first = sp.simplify(first_raw / sp.sqrt(first_norm_squared))
    second_raw = sp.simplify(unit_direction.cross(first))
    second_norm_squared = sp.simplify(second_raw.dot(second_raw))
    second = sp.simplify(second_raw / sp.sqrt(second_norm_squared))
    plus = sp.simplify((first * first.T - second * second.T) / sp.sqrt(2))
    cross = sp.simplify((first * second.T + second * first.T) / sp.sqrt(2))
    return TTPolarizationBasis(
        direction=unit_direction,
        first_transverse=first,
        second_transverse=second,
        plus=plus,
        cross=cross,
    )


def axisymmetric_stf_tensor(
    axial_amplitude: Any,
    symmetry_axis: Any,
    quadrupole_scale: Any = 1,
) -> sp.Matrix:
    """Return ``scale*alpha*(e*e.T-I/3)`` for a nonzero symmetry axis.

    Scale one is the normalized STF convention. Scale three is the
    triple-normalized convention ``Q=3*I_STF``. Keeping the scale explicit
    prevents a tensor constructed in one convention from silently receiving
    the waveform or power coefficient of another.
    """

    amplitude = sp.sympify(axial_amplitude)
    scale = sp.sympify(quadrupole_scale)
    if sp.simplify(scale) == 0:
        raise ValueError("quadrupole_scale must be nonzero")
    axis = _column_three(symmetry_axis, "symmetry_axis")
    axis_norm_squared = sp.simplify(axis.dot(axis))
    if axis_norm_squared == 0:
        raise ValueError("symmetry_axis must be nonzero")
    unit_axis = sp.simplify(axis / sp.sqrt(axis_norm_squared))
    return sp.simplify(
        scale * amplitude * (unit_axis * unit_axis.T - sp.eye(3) / 3)
    )


def axisymmetric_stf_readout(
    axial_amplitude: Any,
    symmetry_axis: Any,
    direction: Any,
    quadrupole_scale: Any = 1,
) -> AxisymmetricSTFReadout:
    """Project an arbitrary-axis axisymmetric STF tensor in its natural frame.

    The first transverse vector is the projected symmetry axis, so the natural
    meridian frame has zero cross coordinate. When the line of sight is
    parallel to the symmetry axis, the TT tensor is zero and a deterministic
    fallback transverse frame is used only to represent that null tensor.
    """

    amplitude = sp.sympify(axial_amplitude)
    scale = sp.sympify(quadrupole_scale)
    axis = _column_three(symmetry_axis, "symmetry_axis")
    line_of_sight = _column_three(direction, "direction")
    axis_norm_squared = sp.simplify(axis.dot(axis))
    direction_norm_squared = sp.simplify(line_of_sight.dot(line_of_sight))
    if axis_norm_squared == 0:
        raise ValueError("symmetry_axis must be nonzero")
    if direction_norm_squared == 0:
        raise ValueError("direction must be nonzero")
    unit_axis = sp.simplify(axis / sp.sqrt(axis_norm_squared))
    unit_direction = sp.simplify(
        line_of_sight / sp.sqrt(direction_norm_squared)
    )
    cosine = sp.simplify(unit_axis.dot(unit_direction))
    sine_squared = sp.simplify(1 - cosine**2)
    first_raw = sp.simplify(unit_axis - cosine * unit_direction)
    if sp.simplify(first_raw.dot(first_raw)) == 0:
        basis = tt_polarization_basis(unit_direction)
    else:
        basis = tt_polarization_basis(unit_direction, first_raw)
    tensor = axisymmetric_stf_tensor(amplitude, unit_axis, scale)
    projected = tt_project_symmetric(tensor, unit_direction)
    plus = sp.simplify(frobenius_inner_product(projected, basis.plus))
    cross = sp.simplify(frobenius_inner_product(projected, basis.cross))
    return AxisymmetricSTFReadout(
        axial_amplitude=amplitude,
        quadrupole_scale=scale,
        symmetry_axis=unit_axis,
        direction=unit_direction,
        inclination_cosine=cosine,
        inclination_sine_squared=sine_squared,
        first_transverse=basis.first_transverse,
        second_transverse=basis.second_transverse,
        source_tensor=tensor,
        projected_tensor=projected,
        normalized_plus_coordinate=plus,
        normalized_cross_coordinate=cross,
        conventional_plus_readout=sp.simplify(plus / sp.sqrt(2)),
        conventional_cross_readout=sp.simplify(cross / sp.sqrt(2)),
    )


def conditional_axisymmetric_stf_power(
    axial_third_derivative: Any,
    gravitational_coupling: Any,
    quadrupole_scale: Any = 1,
) -> sp.Expr:
    """Return the conditional power for an axisymmetric STF derivative.

    This imports the *declared* normalized waveform coefficient ``2G`` and
    flux coefficient ``1/(32*pi*G)`` used by :func:`conditional_tt_power`.
    ``axial_third_derivative`` is the derivative of the normalized underlying
    ``alpha`` in ``I_STF=alpha*(e*e.T-I/3)``. A scale-three tensor therefore
    carries three times the tensor components but receives waveform coefficient
    ``2G/3``; the result is convention invariant and equals
    ``2*G*alpha'''**2/15``. This algebra does not supply a gravity theory.
    """

    derivative = sp.sympify(axial_third_derivative)
    coupling = sp.sympify(gravitational_coupling)
    scale = sp.sympify(quadrupole_scale)
    if sp.simplify(coupling) == 0:
        raise ValueError("gravitational_coupling must be nonzero")
    tensor = axisymmetric_stf_tensor(derivative, [1, 0, 0], scale)
    waveform_prefactor = waveform_prefactor_for_quadrupole_convention(
        2 * coupling, scale
    )
    return conditional_tt_power(
        tensor,
        waveform_prefactor,
        1 / (32 * sp.pi * coupling),
    )


def tt_basis_reconstruct(tensor: Any, basis: TTPolarizationBasis) -> sp.Matrix:
    """Project a symmetric tensor using normalized plus/cross coefficients."""

    source = _symmetric_three(tensor, "tensor")
    return sp.simplify(
        frobenius_inner_product(source, basis.plus) * basis.plus
        + frobenius_inner_product(source, basis.cross) * basis.cross
    )


def rotated_tt_polarizations(
    basis: TTPolarizationBasis,
    angle: Any,
) -> tuple[sp.Matrix, sp.Matrix]:
    """Return plus/cross after rotating the transverse frame by ``angle``.

    With ``u'=cos(angle)u+sin(angle)v`` and
    ``v'=-sin(angle)u+cos(angle)v``, the tensor basis rotates through twice the
    frame angle.
    """

    value = sp.sympify(angle)
    cosine = sp.cos(2 * value)
    sine = sp.sin(2 * value)
    plus = sp.simplify(cosine * basis.plus + sine * basis.cross)
    cross = sp.simplify(-sine * basis.plus + cosine * basis.cross)
    return plus, cross


def circular_tt_polarizations(
    basis: TTPolarizationBasis,
) -> tuple[sp.Matrix, sp.Matrix]:
    """Return ``(plus+i*cross)/sqrt(2)`` and its conjugate convention."""

    right = sp.simplify((basis.plus + sp.I * basis.cross) / sp.sqrt(2))
    left = sp.simplify((basis.plus - sp.I * basis.cross) / sp.sqrt(2))
    return right, left


def integrated_tt_norm_squared(tensor: Any) -> sp.Expr:
    """Return the exact full-sphere integral of ``|TT_n(tensor)|^2``.

    For a symmetric spatial tensor ``S`` in three dimensions,

    ``integral_S2 |TT_n(S)|^2 dOmega = (8*pi/5)*|STF(S)|^2``.

    The result is invariant under adding a pure trace.
    """

    source = _symmetric_three(tensor, "tensor")
    trace_free = symmetric_trace_free(source)
    return sp.simplify(sp.Rational(8, 5) * sp.pi * frobenius_norm_squared(trace_free))


def conditional_tt_power(
    quadrupole_third_derivative: Any,
    waveform_prefactor: Any,
    flux_prefactor: Any,
) -> sp.Expr:
    """Reduce declared waveform and flux inputs to total angular power.

    The declared inputs are
    ``h_TT=(waveform_prefactor/r)*TT(Q_ddot)`` and
    ``dP/dOmega=flux_prefactor*r**2*<hdot_TT:hdot_TT>``.  This function
    performs only the exact angular contraction for the supplied instantaneous
    or already-averaged third derivative.
    """

    wave = sp.sympify(waveform_prefactor)
    flux = sp.sympify(flux_prefactor)
    return sp.simplify(
        flux * wave**2 * integrated_tt_norm_squared(quadrupole_third_derivative)
    )


def waveform_prefactor_for_quadrupole_convention(
    normalized_stf_waveform_prefactor: Any,
    quadrupole_scale: Any,
) -> sp.Expr:
    """Convert a waveform coefficient when ``Q=scale*I_STF``.

    Preserving the same field requires the coefficient multiplying ``Q`` to be
    the coefficient multiplying ``I_STF`` divided by ``scale``.
    """

    prefactor = sp.sympify(normalized_stf_waveform_prefactor)
    scale = sp.sympify(quadrupole_scale)
    if scale == 0:
        raise ValueError("quadrupole_scale must be nonzero")
    return sp.simplify(prefactor / scale)


def harmonic_stf_third_derivative_average(
    cosine_amplitude: Any,
    sine_amplitude: Any,
    angular_frequency: Any,
) -> sp.Expr:
    """Average ``|d^3/dt^3(C cos wt + S sin wt)|^2`` over one cycle.

    Only the symmetric trace-free parts enter a TT contraction.  Exact cycle
    orthogonality gives ``w**6*(|STF(C)|**2+|STF(S)|**2)/2``.
    """

    cosine = symmetric_trace_free(_symmetric_three(cosine_amplitude, "cosine_amplitude"))
    sine = symmetric_trace_free(_symmetric_three(sine_amplitude, "sine_amplitude"))
    frequency = sp.sympify(angular_frequency)
    return sp.simplify(
        frequency**6
        * (frobenius_norm_squared(cosine) + frobenius_norm_squared(sine))
        / 2
    )
