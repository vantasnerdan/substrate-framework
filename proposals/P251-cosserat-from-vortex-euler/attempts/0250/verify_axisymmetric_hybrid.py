"""Expose the leading n=0 thin-ring hybrid coefficients used in 0250."""

import sympy as sp


theta, phi = sp.symbols("theta phi", real=True)
R, s = sp.symbols("R s", positive=True, real=True)
r = R + s * sp.cos(theta)


def fourier(expr, m):
    return sp.simplify(sp.integrate(expr * sp.exp(sp.I * m * theta),
                                    (theta, 0, 2 * sp.pi)))


for m in (1, -1):
    g_u = fourier(r, m)
    g_p = fourier(r**2, m)
    b_parallel = fourier(
        sp.sin(theta) * (1 - 3 * sp.cos(theta)**2), m)
    # Differentiate the Cartesian material expression (8) before either
    # angular integral.  x=r*cos(phi), xi_x=-C*V*sin(theta)*cos(phi)+...
    # and the omitted azimuthal part is odd in phi.
    b_transverse_cartesian = (
        (sp.sin(theta) / s) * (r * sp.cos(phi))**2
        + 2 * (-sp.cos(theta)) * (r * sp.cos(phi))
        * (-sp.sin(theta) * sp.cos(phi))
    )
    b_transverse = sp.simplify(sp.integrate(
        sp.integrate(b_transverse_cartesian * sp.exp(sp.I * m * theta),
                     (phi, 0, 2 * sp.pi)),
        (theta, 0, 2 * sp.pi),
    ))
    velocity_only_mutation = fourier(sp.sin(theta)**3, m)

    assert sp.simplify(g_u - sp.pi * s) == 0
    assert sp.simplify(g_p - 2 * sp.pi * R * s) == 0
    assert sp.simplify(b_parallel - sp.I * sp.pi * m / 4) == 0
    assert sp.simplify(
        b_transverse
        - sp.I * sp.pi**2 * m * (R**2 / s + 3 * s / 4)
    ) == 0
    assert sp.simplify(velocity_only_mutation - b_parallel) != 0

ratio = sp.simplify(
    (sp.pi**2 * (R**2 / s + 3 * s / 4))
    / (sp.pi**2 * s / 2)
)
assert sp.simplify(ratio - (2 * R**2 / s**2 + sp.Rational(3, 2))) == 0

print("PASS axial G coefficients and exact G-cancellation ratio")
print("PASS complete moved-material B_zzz coefficient; velocity-only mutation exposed")
print("PASS transverse/parallel hybrid ratio =", ratio)
print("Scope: leading circular thin-ring coefficients; exact lift and pressure ordering are analytic inputs")
