"""Exact finite-packet coefficient constructions from P251/0147.

Coefficient infrastructure supporting C-CST-011 at its registered parameters;
generic parameter extensions remain unpromoted infrastructure. These functions
derive finite polynomial, Gaussian and angular moment identities. They do
not solve Euler, prove packet localization, supply a physical tag, or bound
the PDE/geometry remainder. Those licenses belong to the declared analytic
construction. Full phase actions use :mod:`euler_phase`, not a new inertia.
"""

from dataclasses import dataclass

import sympy as sp


def _order(value, name):
    item = sp.sympify(value)
    if not isinstance(item, sp.Integer) or item < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(item)


def _finite(value, name, *, positive=False):
    item = sp.sympify(value)
    if item.has(sp.nan, sp.zoo, sp.oo, -sp.oo) or item.is_finite is False:
        raise ValueError(f"{name} must be finite")
    if positive and (item.is_nonpositive is True or item.is_real is False):
        raise ValueError(f"{name} must be positive and real")
    return item


def laguerre_packet_angle(excitation, laplace_variable):
    """Normalized integral of x**2 exp(-w*x/2) L_n^1(x) on x>=0.

The integral representation assumes Re(w)>0. The returned rational
function also gives its meromorphic continuation; w=0 is rejected.
Normalization is its actual value at w=1. Integer n>=0 is required,
including the n=0 case where a factored closed form has cancellations.
No physical clock or frequency is passed in as a target value.
"""
    n = _order(excitation, "excitation")
    w = _finite(laplace_variable, "Laplace variable")
    if w.is_zero is True:
        raise ValueError("Laplace variable must be nonzero")
    x = sp.Dummy("x")
    z = sp.Dummy("w")
    polynomial = sp.Poly(sp.assoc_laguerre(n, 1, x), x)
    integral = sum(coefficient * sp.factorial(power[0] + 2)
                   * (2 / z) ** (power[0] + 3)
                   for power, coefficient in polynomial.terms())
    normalized = sp.cancel(integral / integral.subs(z, 1))
    return sp.factor(normalized.subs(z, w))


@dataclass(frozen=True)
class GaussianCarrierFilter:
    """Uncut packet/marker overlap and whole-packet Plancherel measure.

    The overlap is L exp(-L²(q-p)²/2) times the normalized material
    marker's Fourier factor exp(-ell_z²(q-p_ref)²/2). Band and spatial
    cutoffs are not silently included. ``plancherel_weight`` is the
    integral of the squared spectral envelope; a varying fiber KKS must
    still be integrated against that envelope, not multiplied by this
    constant-profile result.
    """

    center: sp.Expr
    variance: sp.Expr
    carrier_factor: sp.Expr
    amplitude: sp.Expr
    envelope_cost: sp.Expr
    plancherel_weight: sp.Expr


def gaussian_carrier_filter(packet_length, marker_length, carrier, reference):
    """Complete the actual Gaussian overlap square, keeping finite lengths.

    Lengths are fixed during carrier differentiation and strictly positive.
    Symbolic undetermined positivity remains a caller hypothesis; explicit
    invalid lengths are rejected. Carrier and reference are real wave
    numbers. The amplitude includes the carrier-dependent real envelope,
    which cancels from the phase but not from action connections.
    """
    length = _finite(packet_length, "packet length", positive=True)
    marker = _finite(marker_length, "marker length", positive=True)
    p = _finite(carrier, "carrier")
    p0 = _finite(reference, "reference")
    if p.is_real is False or p0.is_real is False:
        raise ValueError("carrier and reference must be real")
    precision = length**2 + marker**2
    center = (length**2 * p + marker**2 * p0) / precision
    cost = sp.factor(length**2 * marker**2 / precision)
    amplitude = length / sp.sqrt(precision) * sp.exp(-cost * (p-p0)**2 / 2)
    return GaussianCarrierFilter(
        sp.factor(center), 1/precision, length**2/precision,
        amplitude, cost, sp.sqrt(sp.pi)*length,
    )


def packet_material_moment_rows(excitation, time_order, coordinate, *, reference_order=3):
    """Construct reference and carrier/time rows of the radial tag moment map.

    Return (1,x,...,x**reference_order) followed by each f_j and D_r f_j,
    with P_n=L_n^1-2 d_x L_n^1, f_j=exp(-x/2)x**j P_n, and
    D_r=-1/2+(3/2)x d_x. Independence of the full family is licensed by
    the analytic n>time_order argument, not presumed by this function.
    Degenerate parameter choices remain available as counterexamples.
    """
    n = _order(excitation, "excitation")
    order = _order(time_order, "time order")
    ref = _order(reference_order, "reference order")
    if not isinstance(coordinate, sp.Symbol):
        raise ValueError("coordinate must be a symbol")
    x = coordinate
    laguerre = sp.assoc_laguerre(n, 1, x)
    pressure = laguerre - 2*sp.diff(laguerre, x)
    rows = [x**j for j in range(ref+1)]
    for j in range(order+1):
        value = sp.exp(-x/2)*x**j*pressure
        rows.extend((value, -value/2 + sp.Rational(3, 2)*x*sp.diff(value, x)))
    return tuple(sp.factor(row) for row in rows)


@dataclass(frozen=True)
class AngularMomentRule:
    """Normalized squared radial frequencies and their exact moment weights.

    Nodes stand for kappa_j**2/lambda**2. The weights match the arcsine
    angular measure; mapping nodes to shared-circle Bessel zeros and
    controlling the resulting field is a separate analytic construction.
    """

    nodes: tuple[sp.Expr, ...]
    weights: tuple[sp.Expr, ...]


def common_circle_angular_rule(count):
    """Positive Gauss-Chebyshev rule on [0,1] through degree 2*count-1."""
    size = _order(count, "count")
    if not size:
        raise ValueError("count must be positive")
    nodes = tuple((1 + sp.cos(sp.pi*(2*j+1)/(2*size)))/2 for j in range(size))
    return AngularMomentRule(nodes, (sp.Rational(1, size),)*size)


def common_circle_moment_weights(nodes):
    """Solve the exact first len(nodes) arcsine moments on specified nodes.

    Distinct nodes lie strictly between zero and one. The returned weights
    need not be positive: positivity after Bessel-node perturbation is a
    separate continuity condition and is not enforced by clipping or
    replacing the actual moment solution. Symbolic domain/determinant
    conditions that cannot be decided remain explicit caller hypotheses.
    """
    points = tuple(_finite(node, "node") for node in nodes)
    if not points:
        raise ValueError("at least one node is required")
    if any(point.is_real is False or point.is_nonpositive is True
           or (point-1).is_nonnegative is True for point in points):
        raise ValueError("nodes must lie strictly between zero and one")
    matrix = sp.Matrix([[point**j for point in points] for j in range(len(points))])
    if sp.simplify(matrix.det()).is_zero is True:
        raise ValueError("nodes must be distinct")
    moments = sp.Matrix([sp.binomial(2*j, j)/4**j for j in range(len(points))])
    weights = matrix.inv()*moments
    return AngularMomentRule(points, tuple(sp.factor(weight) for weight in weights))
