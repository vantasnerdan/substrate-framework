"""Nonlinear spherical compactness for typed confined-clock sources.

Units and signs (``c=1``):

- ``r`` is the areal radius (length).
- ``rho`` is the inertial-frame mass-energy density (mass / length**3).
- ``m(r)`` is the enclosed matter mass, not the geometric mass ``G*m``.
- ``G`` has units length / mass and is positive.
- ``m'(r)=4*pi*r**2*rho`` and ``f(r)=1-2*G*m(r)/r``.

Thus positive density increases the Misner--Sharp compactness
``C(r)=2*G*m(r)/r`` and a sign change of ``f=1-C`` is a trapped-surface
obstruction to a globally static horizonless areal-gauge metric.  The numeric
finite-volume API integrates density itself; it never obtains a small mass by
subtracting large energies.

The homothetic API keeps two separately scaling nonnegative source ledgers:
a quartic-curvature contribution with integrated energy proportional to
``R**-1`` and a potential contribution proportional to ``R**3``.  It does not
claim arbitrary profile relaxation or localize a global fixed-J constraint.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import sympy as sp
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import minimize_scalar

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class ExactSphericalMassConstraint:
    """Exact areal-gauge Hamiltonian constraint and metric function."""

    radius: sp.Expr
    enclosed_mass: sp.Expr
    energy_density: sp.Expr
    newton_constant: sp.Expr
    mass_constraint_residual: sp.Expr
    compactness: sp.Expr
    radial_metric_function: sp.Expr


@dataclass(frozen=True)
class SphericalCompactnessProfile:
    """Finite-volume Misner--Sharp evidence on radial cell edges."""

    radius_edges: FloatArray
    density_cells: FloatArray
    enclosed_mass: FloatArray
    compactness: FloatArray
    radial_metric_function: FloatArray
    newton_constant: float

    @property
    def total_mass(self) -> float:
        return float(self.enclosed_mass[-1])

    @property
    def maximum_compactness(self) -> float:
        return float(np.max(self.compactness))

    @property
    def minimum_radial_metric_function(self) -> float:
        return float(np.min(self.radial_metric_function))

    @property
    def exterior_horizon_radius(self) -> float:
        """Return ``2*G*M`` for the matched vacuum Schwarzschild exterior."""

        return float(2.0 * self.newton_constant * self.total_mass)

    @property
    def critical_newton_constant(self) -> float:
        """Largest ``G`` with ``2*G*m(r)/r < 1`` on sampled positive edges."""

        radii = self.radius_edges[1:]
        masses = self.enclosed_mass[1:]
        positive = masses > 0.0
        if not np.any(positive):
            return float("inf")
        return float(np.min(radii[positive] / (2.0 * masses[positive])))

    def first_trapped_surface_radius(self) -> float | None:
        """Linearly bracket the first inward-static ``f=0`` crossing."""

        f = self.radial_metric_function
        exact = np.flatnonzero(f == 0.0)
        if exact.size:
            return float(self.radius_edges[int(exact[0])])
        crossings = np.flatnonzero((f[:-1] > 0.0) & (f[1:] < 0.0))
        if not crossings.size:
            return None
        index = int(crossings[0])
        left_r, right_r = self.radius_edges[index:index + 2]
        left_f, right_f = f[index:index + 2]
        fraction = left_f / (left_f - right_f)
        return float(left_r + fraction * (right_r - left_r))


@dataclass(frozen=True)
class HomotheticCompactnessMinimum:
    """Global-in-bracket minimum of maximum compactness versus log scale."""

    scale_radius: float
    maximum_compactness: float
    critical_newton_constant: float
    lower_scale_compactness: float
    upper_scale_compactness: float
    optimizer_success: bool
    function_evaluations: int


def _positive_exact(value: Any, name: str) -> sp.Expr:
    result = sp.sympify(value)
    if result.has(sp.Float):
        raise ValueError(f"{name} must be exact")
    if result.is_positive is not True:
        raise ValueError(f"{name} must be provably positive")
    return sp.simplify(result)


def static_spherical_mass_constraint(
    radius: Any,
    enclosed_mass: Any,
    energy_density: Any,
    newton_constant: Any,
) -> ExactSphericalMassConstraint:
    r"""Return ``m'-4*pi*r**2*rho`` and ``f=1-2*G*m/r`` exactly."""

    r = sp.sympify(radius)
    if not isinstance(r, sp.Symbol) or r.is_positive is not True:
        raise ValueError("radius must be a positive exact symbol")
    mass = sp.sympify(enclosed_mass)
    density = sp.sympify(energy_density)
    gravity = _positive_exact(newton_constant, "newton_constant")
    compactness = sp.simplify(2 * gravity * mass / r)
    return ExactSphericalMassConstraint(
        radius=r,
        enclosed_mass=mass,
        energy_density=density,
        newton_constant=gravity,
        mass_constraint_residual=sp.simplify(
            sp.diff(mass, r) - 4 * sp.pi * r**2 * density
        ),
        compactness=compactness,
        radial_metric_function=sp.simplify(1 - compactness),
    )


def integrate_spherical_density_cells(
    radius_edges: ArrayLike,
    density_cells: ArrayLike,
    newton_constant: float,
) -> SphericalCompactnessProfile:
    """Integrate piecewise-constant cell densities with exact shell volumes."""

    edges = np.asarray(radius_edges, dtype=np.float64)
    density = np.asarray(density_cells, dtype=np.float64)
    gravity = float(newton_constant)
    if edges.ndim != 1 or density.ndim != 1 or edges.size != density.size + 1:
        raise ValueError(
            "radius_edges must have exactly one more entry than density_cells"
        )
    if edges.size < 2 or edges[0] != 0.0 or np.any(np.diff(edges) <= 0.0):
        raise ValueError(
            "radius_edges must start at zero and increase strictly"
        )
    if not np.all(np.isfinite(edges)) or not np.all(np.isfinite(density)):
        raise ValueError("radius and density must be finite")
    if np.any(density < 0.0):
        raise ValueError("energy density must be nonnegative")
    if not np.isfinite(gravity) or gravity <= 0.0:
        raise ValueError("newton_constant must be positive and finite")

    shell_volume = (4.0 * np.pi / 3.0) * (edges[1:] ** 3 - edges[:-1] ** 3)
    shell_mass = density * shell_volume
    enclosed = np.concatenate(([0.0], np.cumsum(shell_mass, dtype=np.float64)))
    compactness = np.zeros_like(edges)
    compactness[1:] = 2.0 * gravity * enclosed[1:] / edges[1:]
    radial_metric = 1.0 - compactness
    return SphericalCompactnessProfile(
        radius_edges=edges,
        density_cells=density,
        enclosed_mass=enclosed,
        compactness=compactness,
        radial_metric_function=radial_metric,
        newton_constant=gravity,
    )


def homothetic_compactness_profile(
    dimensionless_radius: ArrayLike,
    cumulative_curvature_mass: ArrayLike,
    cumulative_potential_mass: ArrayLike,
    *,
    reference_radius: float,
    scale_radius: float,
    newton_constant: float,
) -> FloatArray:
    """Return compactness for the exact curvature/potential scale ledger."""

    x = np.asarray(dimensionless_radius, dtype=np.float64)
    curvature = np.asarray(cumulative_curvature_mass, dtype=np.float64)
    potential = np.asarray(cumulative_potential_mass, dtype=np.float64)
    if x.ndim != 1 or curvature.shape != x.shape or potential.shape != x.shape:
        raise ValueError(
            "dimensionless radius and cumulative masses must be "
            "equal 1-D arrays"
        )
    if x.size < 2 or x[0] != 0.0 or np.any(np.diff(x) <= 0.0):
        raise ValueError(
            "dimensionless_radius must start at zero and increase"
        )
    if np.any(curvature < 0.0) or np.any(potential < 0.0):
        raise ValueError("component masses must be nonnegative")
    if np.any(np.diff(curvature) < 0.0) or np.any(np.diff(potential) < 0.0):
        raise ValueError(
            "component masses must be cumulative and nondecreasing"
        )
    reference = float(reference_radius)
    scale = float(scale_radius)
    gravity = float(newton_constant)
    if not all(
        np.isfinite(value) and value > 0.0
        for value in (reference, scale, gravity)
    ):
        raise ValueError(
            "radii and newton_constant must be positive and finite"
        )

    scaled_mass = (
        curvature * (reference / scale)
        + potential * (scale / reference) ** 3
    )
    compactness = np.zeros_like(x)
    compactness[1:] = 2.0 * gravity * scaled_mass[1:] / (scale * x[1:])
    return compactness


def minimize_homothetic_max_compactness(
    dimensionless_radius: ArrayLike,
    cumulative_curvature_mass: ArrayLike,
    cumulative_potential_mass: ArrayLike,
    *,
    reference_radius: float,
    newton_constant: float,
    scale_bounds: tuple[float, float],
    log_tolerance: float = 1.0e-10,
) -> HomotheticCompactnessMinimum:
    """Minimize convex maximum compactness over a declared scale bracket."""

    lower, upper = map(float, scale_bounds)
    if not 0.0 < lower < upper:
        raise ValueError("scale_bounds must be positive and ordered")

    def objective(
        log_scale: float,
        gravity: float = float(newton_constant),
    ) -> float:
        profile = homothetic_compactness_profile(
            dimensionless_radius,
            cumulative_curvature_mass,
            cumulative_potential_mass,
            reference_radius=reference_radius,
            scale_radius=float(np.exp(log_scale)),
            newton_constant=gravity,
        )
        return float(np.max(profile))

    result = minimize_scalar(
        objective,
        bounds=(float(np.log(lower)), float(np.log(upper))),
        method="bounded",
        options={"xatol": float(log_tolerance), "maxiter": 1000},
    )
    scale = float(np.exp(result.x))
    maximum = float(result.fun)
    gravity = float(newton_constant)
    return HomotheticCompactnessMinimum(
        scale_radius=scale,
        maximum_compactness=maximum,
        critical_newton_constant=float(gravity / maximum),
        lower_scale_compactness=objective(float(np.log(lower))),
        upper_scale_compactness=objective(float(np.log(upper))),
        optimizer_success=bool(result.success),
        function_evaluations=int(result.nfev),
    )
