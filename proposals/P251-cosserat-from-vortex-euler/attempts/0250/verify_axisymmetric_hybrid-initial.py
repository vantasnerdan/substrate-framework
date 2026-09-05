"""Initial check preserved after its transverse phi-normalization bug was found."""

import sympy as sp


theta = sp.symbols("theta", real=True)
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
    b_transverse = fourier(
        r**2 * sp.sin(theta) / (2 * s)
        + r * sp.cos(theta) * sp.sin(theta), m)
    velocity_only_mutation = fourier(sp.sin(theta)**3, m)

    assert sp.simplify(g_u - sp.pi * s) == 0
    assert sp.simplify(g_p - 2 * sp.pi * R * s) == 0
    assert sp.simplify(b_parallel - sp.I * sp.pi * m / 4) == 0
    assert sp.simplify(
        b_transverse
        - sp.I * sp.pi * m * (R**2 / (2 * s) + 3 * s / 8)
    ) == 0
    assert sp.simplify(velocity_only_mutation - b_parallel) != 0

ratio = sp.simplify(
    (sp.pi * (R**2 / (2 * s) + 3 * s / 8))
    / (2 * sp.pi * s / 4)
)
assert sp.simplify(ratio - (R**2 / s**2 + sp.Rational(3, 4))) == 0

print("PASS axial G coefficients and exact G-cancellation ratio")
print("PASS complete moved-material B_zzz coefficient; velocity-only mutation exposed")
print("PASS transverse/parallel hybrid ratio =", ratio)
print("Scope: leading circular thin-ring coefficients; exact lift and pressure ordering are analytic inputs")
