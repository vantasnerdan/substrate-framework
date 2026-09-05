"""Exact conditional Bernoulli lift of stationary planar Euler (P251/0136).

This additive infrastructure is unpromoted. Real smooth planar fields and
constant positive density are assumed. The caller selects a spatial domain
where C > B = p/rho + |v|²/2. A global lift needs that inequality globally;
the symbolic construction does not prove boundedness on an unspecified domain.
The result is generalized Beltrami, not necessarily constant-factor Beltrami,
finite total energy, knotted, or an autonomous Cosserat continuum.
"""

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class PlanarBernoulliLift:
    """Derived field and explicit domain condition, with physical pressure."""

    velocity: sp.ImmutableMatrix
    pressure: sp.Expr
    planar_bernoulli: sp.Expr
    axial_speed: sp.Expr
    curl_factor: sp.Expr
    domain_condition: sp.Expr


def planar_bernoulli_lift(velocity, pressure, density, coordinates, level):
    """Lift verified steady planar Euler to ``curl u = f*u`` in three dimensions.

    Return u=(v_x,v_y,sqrt(2*(C-B))), f=(d_x v_y-d_y v_x)/u_z.
    Pressure is physical p; C and B have velocity-squared units. Density
    and C are spatial constants. Incompressibility and stationary Euler
    are checked symbolically; unresolved residuals are rejected rather than
    treated as proved. Positivity of 2*(C-B) on the caller's domain is an
    explicit hypothesis (a provably nonpositive expression is rejected).

    No globally single-valued relation zeta=F(psi) is needed. The actual
    global Bernoulli function supplies the first integral, including when
    different streamline components carry different vorticity laws.
    """
    coords = tuple(coordinates)
    if len(coords) != 2 or len(set(coords)) != 2 or not all(
        isinstance(c, sp.Symbol) for c in coords
    ):
        raise ValueError("two distinct Cartesian coordinate symbols are required")
    flow = sp.ImmutableMatrix(velocity)
    if flow.shape != (2, 1):
        raise ValueError("velocity must be a planar column vector")
    p, rho, constant = map(sp.sympify, (pressure, density, level))
    if rho.is_positive is False or rho.is_finite is False or (
        rho.is_number and rho.is_positive is not True
    ):
        raise ValueError("density must be finite and positive")
    if constant.has(sp.nan, sp.zoo, sp.oo, -sp.oo) or (
        constant.is_real is False or constant.is_finite is False
    ):
        raise ValueError("Bernoulli level must be finite and real")
    if any(sp.diff(value, c) != 0 for value in (rho, constant) for c in coords):
        raise ValueError("density and Bernoulli level must be spatial constants")
    if any(value.has(sp.nan, sp.zoo, sp.oo, -sp.oo) or value.is_real is False
           for value in (*flow, p)):
        raise ValueError("velocity and pressure must be finite real symbolic fields")
    gradient_p = sp.ImmutableMatrix([sp.diff(p, c) for c in coords])
    divergence = sum(sp.diff(flow[i], coords[i]) for i in range(2))
    residual = rho*flow.jacobian(coords)*flow+gradient_p
    if sp.simplify(divergence) != 0 or any(sp.simplify(v) != 0 for v in residual):
        raise ValueError("stationary incompressible planar Euler residual is not verified")
    bernoulli = sp.simplify(p/rho+flow.dot(flow)/2)
    radicand = sp.simplify(2*(constant-bernoulli))
    if radicand.is_nonpositive is True:
        raise ValueError("Bernoulli level must be strictly above B on the chosen domain")
    axial = sp.sqrt(radicand)
    vorticity = sp.diff(flow[1], coords[0])-sp.diff(flow[0], coords[1])
    return PlanarBernoulliLift(
        sp.ImmutableMatrix([*flow, axial]), p, bernoulli, axial,
        sp.simplify(vorticity/axial), sp.Gt(radicand, 0),
    )
