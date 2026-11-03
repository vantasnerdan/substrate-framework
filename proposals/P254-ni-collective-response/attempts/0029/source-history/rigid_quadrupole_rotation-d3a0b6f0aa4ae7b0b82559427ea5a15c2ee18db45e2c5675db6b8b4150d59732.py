"""Exact kinematics of rigidly rotated axisymmetric STF moments.

The APIs in this module prescribe a tensor path under a declared spatial
rotation.  They do not establish that the path solves a field equation, that
the angular speed is dynamically selected, or that any gravitational theory
couples to the moment.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


def _column_three(vector: Any, name: str) -> sp.Matrix:
    value = sp.Matrix(vector)
    if value.shape not in ((3, 1), (1, 3)):
        raise ValueError(f"{name} must have three components")
    return value if value.shape == (3, 1) else value.T


def _unit_vector(vector: Any, name: str) -> sp.Matrix:
    value = _column_three(vector, name)
    norm_squared = sp.simplify(value.dot(value))
    if norm_squared == 0:
        raise ValueError(f"{name} must be nonzero")
    return sp.simplify(value / sp.sqrt(norm_squared))


def rodrigues_rotation_matrix(axis: Any, angle: Any) -> sp.Matrix:
    """Return the exact right-handed Rodrigues rotation about ``axis``.

    For the positive Cartesian x axis this convention gives
    ``[[1,0,0],[0,cos,-sin],[0,sin,cos]]``.
    """

    unit_axis = _unit_vector(axis, "axis")
    x, y, z = unit_axis
    cross = sp.Matrix([[0, -z, y], [z, 0, -x], [-y, x, 0]])
    cosine = sp.cos(sp.sympify(angle))
    sine = sp.sin(sp.sympify(angle))
    return sp.simplify(
        cosine * sp.eye(3)
        + (1 - cosine) * unit_axis * unit_axis.T
        + sine * cross
    )


def axisymmetric_stf_from_transverse_eigenvalue(
    transverse_eigenvalue: Any,
    symmetry_axis: Any,
    quadrupole_scale: Any = 1,
) -> sp.Matrix:
    """Return an axisymmetric STF tensor with repeated transverse eigenvalue.

    If ``q`` is ``transverse_eigenvalue`` and ``e`` is the unit symmetry axis,
    the normalized tensor is ``q*(I-3*e*e.T)`` and has eigenvalues
    ``(q,q,-2q)``.  ``quadrupole_scale=3`` returns the corresponding
    ``Q=3*I_STF`` convention without relabeling the normalized tensor.
    """

    value = sp.sympify(transverse_eigenvalue)
    scale = sp.sympify(quadrupole_scale)
    if sp.simplify(scale) == 0:
        raise ValueError("quadrupole_scale must be nonzero")
    axis = _unit_vector(symmetry_axis, "symmetry_axis")
    return sp.simplify(scale * value * (sp.eye(3) - 3 * axis * axis.T))


@dataclass(frozen=True)
class RigidAxisymmetricSTFRotation:
    """A prescribed rigid rotation and its first three exact derivatives."""

    transverse_eigenvalue: sp.Expr
    quadrupole_scale: sp.Expr
    angular_speed: sp.Expr
    time: sp.Symbol
    rotation_axis: sp.Matrix
    body_symmetry_axis: sp.Matrix
    instantaneous_symmetry_axis: sp.Matrix
    rotation_matrix: sp.Matrix
    body_tensor: sp.Matrix
    tensor: sp.Matrix
    first_derivative: sp.Matrix
    second_derivative: sp.Matrix
    third_derivative: sp.Matrix


def rigid_axisymmetric_stf_rotation(
    transverse_eigenvalue: Any,
    angular_speed: Any,
    time: sp.Symbol,
    *,
    rotation_axis: Any,
    body_symmetry_axis: Any,
    quadrupole_scale: Any = 1,
) -> RigidAxisymmetricSTFRotation:
    """Return ``R(t) S_body R(t).T`` and three time derivatives exactly.

    The rotation is kinematic input.  Nonzero derivatives require both a
    nonzero anisotropy and a nonzero angular speed, and do not by themselves
    establish a dynamically allowed rotating source or physical radiation.
    """

    if not isinstance(time, sp.Symbol):
        raise ValueError("time must be a SymPy Symbol")
    value = sp.sympify(transverse_eigenvalue)
    speed = sp.sympify(angular_speed)
    scale = sp.sympify(quadrupole_scale)
    axis = _unit_vector(rotation_axis, "rotation_axis")
    body_axis = _unit_vector(body_symmetry_axis, "body_symmetry_axis")
    rotation = rodrigues_rotation_matrix(axis, speed * time)
    instantaneous_axis = sp.simplify(rotation * body_axis)
    body = axisymmetric_stf_from_transverse_eigenvalue(
        value, body_axis, scale
    )
    tensor = sp.simplify(rotation * body * rotation.T)
    derivatives = tuple(
        sp.simplify(sp.diff(tensor, time, order)) for order in (1, 2, 3)
    )
    return RigidAxisymmetricSTFRotation(
        transverse_eigenvalue=value,
        quadrupole_scale=scale,
        angular_speed=speed,
        time=time,
        rotation_axis=axis,
        body_symmetry_axis=body_axis,
        instantaneous_symmetry_axis=instantaneous_axis,
        rotation_matrix=rotation,
        body_tensor=body,
        tensor=tensor,
        first_derivative=derivatives[0],
        second_derivative=derivatives[1],
        third_derivative=derivatives[2],
    )


def tilted_axisymmetric_stf_rotation_about_z(
    transverse_eigenvalue: Any,
    angular_speed: Any,
    time: sp.Symbol,
    tilt: Any,
    quadrupole_scale: Any = 1,
) -> RigidAxisymmetricSTFRotation:
    """Return the canonical constant-tilt rotation about the z axis.

    The body symmetry axis starts in the x-z plane at angle ``tilt`` from the
    rotation axis.  Generic tilt contains both angular-speed and twice-angular-
    speed tensor harmonics; the aligned path is constant and the exactly
    perpendicular path contains only the twice-angular-speed harmonic plus DC.
    """

    angle = sp.sympify(tilt)
    return rigid_axisymmetric_stf_rotation(
        transverse_eigenvalue,
        angular_speed,
        time,
        rotation_axis=[0, 0, 1],
        body_symmetry_axis=[sp.sin(angle), 0, sp.cos(angle)],
        quadrupole_scale=quadrupole_scale,
    )


def symmetric_tensor_characteristic_polynomial(
    tensor: Any,
    spectral_parameter: Any,
) -> sp.Expr:
    """Return ``det(lambda*I-S)`` for a symmetric three-tensor exactly."""

    value = sp.Matrix(tensor)
    if value.shape != (3, 3):
        raise ValueError("tensor must be 3 by 3")
    if sp.simplify(value - value.T) != sp.zeros(3):
        raise ValueError("tensor must be symmetric")
    parameter = sp.sympify(spectral_parameter)
    return sp.factor((parameter * sp.eye(3) - value).det())
