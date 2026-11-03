"""Conditional radial Abelian-Higgs vortex equations and BVP evidence.

The convention is ``phi=f(r) exp(i*n*theta)`` and
``A_theta=a(r)/(g*r)``.  This module supplies no substrate, chromoelectric,
QCD, dual-superconductor, or confinement interpretation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray
import sympy as sp

from .numerics import BVPEvidence, solve_bvp_evidence, trapezoid_integral

FloatArray = NDArray[np.float64]


def _positive_symbolic(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    return expression


def _positive_float(value: float, name: str) -> float:
    result = float(value)
    if not np.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return result


def _positive_winding(value: int) -> int:
    if isinstance(value, bool) or int(value) != value or int(value) <= 0:
        raise ValueError("winding must be a positive integer")
    return int(value)


@dataclass(frozen=True)
class VortexParameters:
    """Positive parameters for the declared radial vortex convention."""

    vacuum_scale: float = 1.0
    winding: int = 1
    self_coupling: float = 2.0
    gauge_coupling: float = 1.0

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "vacuum_scale", _positive_float(self.vacuum_scale, "vacuum_scale")
        )
        object.__setattr__(self, "winding", _positive_winding(self.winding))
        object.__setattr__(
            self,
            "self_coupling",
            _positive_float(self.self_coupling, "self_coupling"),
        )
        object.__setattr__(
            self,
            "gauge_coupling",
            _positive_float(self.gauge_coupling, "gauge_coupling"),
        )


@dataclass(frozen=True)
class VortexSolution:
    """A converged truncated-domain solution and its declared setup."""

    parameters: VortexParameters
    inner_radius: float
    outer_radius: float
    tolerance: float
    initial_points: int
    guess_family: str
    evidence: BVPEvidence

    def state_at(self, coordinate: ArrayLike) -> FloatArray:
        """Linearly interpolate ``[f,f',a,a']`` on the converged mesh."""

        points = np.asarray(coordinate, dtype=np.float64)
        if points.ndim != 1 or not np.all(np.isfinite(points)):
            raise ValueError("coordinate must be a finite one-dimensional array")
        if np.any(points < self.inner_radius) or np.any(points > self.outer_radius):
            raise ValueError("coordinate lies outside the solved radial interval")
        return np.vstack(
            [
                np.interp(points, self.evidence.coordinate, component)
                for component in self.evidence.state
            ]
        )


def radial_energy_lagrangian(
    radius: sp.Symbol,
    scalar_profile: sp.Expr,
    gauge_profile: sp.Expr,
    winding: Any,
    self_coupling: Any,
    vacuum_scale: Any,
    gauge_coupling: Any,
) -> sp.Expr:
    """Return the radial energy integrand including the cylindrical measure."""

    n = sp.sympify(winding)
    lam = _positive_symbolic(self_coupling, "self_coupling")
    vacuum = _positive_symbolic(vacuum_scale, "vacuum_scale")
    gauge = _positive_symbolic(gauge_coupling, "gauge_coupling")
    f, a = sp.sympify(scalar_profile), sp.sympify(gauge_profile)
    return sp.simplify(
        radius
        * (
            sp.diff(f, radius) ** 2 / 2
            + f**2 * (n - a) ** 2 / (2 * radius**2)
            + sp.diff(a, radius) ** 2 / (2 * gauge**2 * radius**2)
            + lam * (f**2 - vacuum**2) ** 2 / 4
        )
    )


def euler_lagrange_residuals(
    radius: sp.Symbol,
    scalar_profile: sp.Expr,
    gauge_profile: sp.Expr,
    winding: Any,
    self_coupling: Any,
    vacuum_scale: Any,
    gauge_coupling: Any,
) -> tuple[sp.Expr, sp.Expr]:
    """Return canonical scalar and gauge Euler-Lagrange residuals."""

    n = sp.sympify(winding)
    lam = _positive_symbolic(self_coupling, "self_coupling")
    vacuum = _positive_symbolic(vacuum_scale, "vacuum_scale")
    gauge = _positive_symbolic(gauge_coupling, "gauge_coupling")
    f, a = sp.sympify(scalar_profile), sp.sympify(gauge_profile)
    scalar = (
        sp.diff(f, radius, 2)
        + sp.diff(f, radius) / radius
        - f * (n - a) ** 2 / radius**2
        - lam * f * (f**2 - vacuum**2)
    )
    gauge_residual = (
        sp.diff(a, radius, 2)
        - sp.diff(a, radius) / radius
        + gauge**2 * (n - a) * f**2
    )
    return sp.simplify(scalar), sp.simplify(gauge_residual)


def angular_log_coefficient(
    vacuum_scale: Any, winding: Any, asymptotic_gauge_profile: Any
) -> sp.Expr:
    """Return the large-radius logarithmic angular-energy coefficient."""

    vacuum = _positive_symbolic(vacuum_scale, "vacuum_scale")
    return sp.simplify(
        vacuum**2
        * (sp.sympify(winding) - sp.sympify(asymptotic_gauge_profile)) ** 2
    )


def quantized_flux(winding: Any, gauge_coupling: Any) -> sp.Expr:
    """Return ``2*pi*n/g`` for ``A_theta=a/(g*r)``."""

    n = sp.sympify(winding)
    if n.is_number and n.is_integer is not True:
        raise ValueError("winding must be an integer")
    gauge = _positive_symbolic(gauge_coupling, "gauge_coupling")
    return sp.simplify(2 * sp.pi * n / gauge)


def asymptotic_masses(
    vacuum_scale: Any, self_coupling: Any, gauge_coupling: Any
) -> tuple[sp.Expr, sp.Expr]:
    """Return exact linearized vector and scalar inverse lengths."""

    vacuum = _positive_symbolic(vacuum_scale, "vacuum_scale")
    lam = _positive_symbolic(self_coupling, "self_coupling")
    gauge = _positive_symbolic(gauge_coupling, "gauge_coupling")
    return sp.simplify(gauge * vacuum), sp.simplify(
        vacuum * sp.sqrt(2 * lam)
    )


def vortex_bvp_equations(
    coordinate: FloatArray, state: FloatArray, parameters: VortexParameters
) -> FloatArray:
    """Return the first-order radial system for ``[f,f',a,a']``."""

    f, fp, a, ap = state
    fpp = (
        -fp / coordinate
        + f * (parameters.winding - a) ** 2 / coordinate**2
        + parameters.self_coupling
        * f
        * (f**2 - parameters.vacuum_scale**2)
    )
    app = (
        ap / coordinate
        - parameters.gauge_coupling**2
        * (parameters.winding - a)
        * f**2
    )
    return np.vstack((fp, fpp, ap, app))


def vortex_boundary_residual(
    left_state: FloatArray, right_state: FloatArray, parameters: VortexParameters
) -> FloatArray:
    """Return truncated-domain Dirichlet boundary residuals."""

    return np.asarray(
        [
            left_state[0],
            left_state[2],
            right_state[0] - parameters.vacuum_scale,
            right_state[2] - parameters.winding,
        ],
        dtype=np.float64,
    )


def _initial_guess(
    coordinate: FloatArray,
    parameters: VortexParameters,
    family: Literal["exponential", "rational"],
) -> FloatArray:
    vacuum, winding = parameters.vacuum_scale, float(parameters.winding)
    scalar_width = 1.0 / (vacuum * np.sqrt(parameters.self_coupling))
    gauge_width = 1.0 / (vacuum * parameters.gauge_coupling)
    if family == "exponential":
        z = coordinate / scalar_width
        f = vacuum * np.tanh(z)
        fp = vacuum / scalar_width / np.cosh(z) ** 2
        tail = np.exp(-(coordinate / gauge_width) ** 2)
        a = winding * (1.0 - tail)
        ap = 2.0 * winding * coordinate * tail / gauge_width**2
    elif family == "rational":
        scalar_denominator = np.sqrt(coordinate**2 + scalar_width**2)
        f = vacuum * coordinate / scalar_denominator
        fp = vacuum * scalar_width**2 / scalar_denominator**3
        gauge_denominator = coordinate**2 + gauge_width**2
        a = winding * coordinate**2 / gauge_denominator
        ap = 2.0 * winding * coordinate * gauge_width**2 / gauge_denominator**2
    else:
        raise ValueError("guess_family must be 'exponential' or 'rational'")
    return np.vstack((f, fp, a, ap))


def solve_vortex_bvp(
    parameters: VortexParameters = VortexParameters(),
    *,
    inner_radius: float = 1.0e-4,
    outer_radius: float = 25.0,
    initial_points: int = 200,
    tolerance: float = 1.0e-7,
    max_nodes: int = 100_000,
    guess_family: Literal["exponential", "rational"] = "exponential",
) -> VortexSolution:
    """Solve the declared truncated radial BVP with retained diagnostics."""

    inner = _positive_float(inner_radius, "inner_radius")
    outer = _positive_float(outer_radius, "outer_radius")
    if outer <= inner:
        raise ValueError("outer_radius must exceed inner_radius")
    if initial_points < 8:
        raise ValueError("initial_points must be at least eight")
    coordinate = np.linspace(inner, outer, int(initial_points))
    guess = _initial_guess(coordinate, parameters, guess_family)
    evidence = solve_bvp_evidence(
        lambda r, y: vortex_bvp_equations(r, y, parameters),
        lambda left, right: vortex_boundary_residual(left, right, parameters),
        coordinate,
        guess,
        tolerance=tolerance,
        max_nodes=max_nodes,
    )
    return VortexSolution(
        parameters,
        inner,
        outer,
        float(tolerance),
        int(initial_points),
        guess_family,
        evidence,
    )


def vortex_energy_density(
    coordinate: ArrayLike, state: ArrayLike, parameters: VortexParameters
) -> FloatArray:
    """Return the nonnegative energy density inside the radial measure."""

    radius = np.asarray(coordinate, dtype=np.float64)
    values = np.asarray(state, dtype=np.float64)
    if radius.ndim != 1 or values.shape != (4, radius.size):
        raise ValueError("state must have shape (4, coordinate size)")
    f, fp, a, ap = values
    return np.asarray(
        fp**2 / 2
        + f**2 * (parameters.winding - a) ** 2 / (2 * radius**2)
        + (ap / radius) ** 2 / (2 * parameters.gauge_coupling**2)
        + parameters.self_coupling
        * (f**2 - parameters.vacuum_scale**2) ** 2
        / 4,
        dtype=np.float64,
    )


def vortex_tension(solution: VortexSolution, *, quadrature_points: int = 20_001) -> float:
    """Integrate energy per length on a declared uniform quadrature mesh."""

    if quadrature_points < 101:
        raise ValueError("quadrature_points must be at least 101")
    radius = np.linspace(
        solution.inner_radius, solution.outer_radius, int(quadrature_points)
    )
    state = solution.state_at(radius)
    density = vortex_energy_density(radius, state, solution.parameters)
    tension = 2.0 * np.pi * trapezoid_integral(radius * density, radius)
    if not np.isfinite(tension) or tension <= 0.0:
        raise ValueError("computed vortex tension must be positive and finite")
    return tension
