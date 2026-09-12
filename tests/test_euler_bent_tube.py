"""Exact pushforward, divergence and relative moment checks."""

import sympy as sp

from substrate_framework.euler_bent_tube import (
    bent_column, bent_column_relative_impulse_density,
)


def test_volume_preserving_pushforward_is_solenoidal_in_eulerian_coordinates():
    x, y, z, R = sp.symbols("x y z R", real=True)
    a, b = sp.Function("a")(z), sp.Function("b")(z)
    w = sp.Function("w")(R)
    data = bent_column((x, y, z), (a, b), w, R)
    assert data.deformation_gradient.det() == 1
    direct = data.deformation_gradient*sp.Matrix([0, 0, w.subs(R, x*x+y*y)])
    direct = direct.subs({x: x-a, y: y-b}, simultaneous=True)
    assert direct == data.vorticity
    divergence = sum(sp.diff(data.vorticity[j], q) for j, q in enumerate((x, y, z)))
    assert sp.simplify(divergence) == 0
    # Omitting transverse pushforward terms leaves a nonzero divergence.
    assert sp.diff(data.vorticity[2], z) != 0


def test_relative_impulse_from_actual_gaussian_cross_section():
    x, y, z, R = sp.symbols("x y z R", real=True)
    a, b, da, db = sp.symbols("a b da db", real=True)
    g = sp.Symbol("Gamma", positive=True)
    A, B = sp.Function("a")(z), sp.Function("b")(z)
    data = bent_column((x, y, z), (A, B), g*sp.exp(-R)/sp.pi, R)
    jets = {sp.diff(A, z): da, sp.diff(B, z): db, A: a, B: b}
    # Translate the integration variables for the bent part, separately from
    # the unbent part, so no cancellation of divergent axial integrals is used.
    omega = data.vorticity.xreplace(jets).subs({x: x+a, y: y+b}, simultaneous=True)
    bent = sp.Matrix([x+a, y+b, z]).cross(omega)/2
    base = sp.Matrix([x, y, z]).cross(sp.Matrix([0, 0, g*sp.exp(-x*x-y*y)/sp.pi]))/2
    actual = sp.Matrix([sp.integrate(sp.expand(entry), (x, -sp.oo, sp.oo),
                                     (y, -sp.oo, sp.oo)) for entry in bent-base])
    expected = bent_column_relative_impulse_density((A, B), z, g).xreplace(jets)
    assert (actual-expected).applyfunc(sp.simplify) == sp.zeros(3, 1)
    # A localized bend has nonzero vector impulse even though int delta_omega=0.
    center = (sp.exp(-z*z), z*sp.exp(-z*z))
    density = bent_column_relative_impulse_density(center, z, g)
    total = density.applyfunc(lambda f: sp.integrate(f, (z, -sp.oo, sp.oo)))
    assert total == sp.Matrix([0, -g*sp.sqrt(sp.pi), g*sp.sqrt(sp.pi/2)/2])
