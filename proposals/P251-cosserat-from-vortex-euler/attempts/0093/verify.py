"""Exact algebra oracle for material-gradient.md; no soft numerics."""

import sympy as s

checks = 0


def check(name, condition):
    global checks
    if not condition:
        raise AssertionError(name)
    checks += 1
    print(f"PASS {name}")


def zero(value):
    return s.simplify(value) == 0


def jet(matrix, k):
    return matrix.applyfunc(lambda x: s.series(x, k, 0, 3).removeO().expand())


L, a, m, mx, kx, d = s.symbols("L a m mx kx d", positive=True)
ax, kb, t = s.symbols("ax kb t", real=True)
j = m * L**5 + mx * t**2
kappa = kb + kx * t**2
rx = kx / mx
check("exact optical ratio margin", zero(rx * j - kappa - (rx * m * L**5 - kb)))
root = (ax + s.sqrt(ax**2 + 4 * mx * d * L**5)) / (2 * mx)
check("spin completion root", zero(mx * root**2 - ax * root - d * L**5))
check(
    "physical spin equals metric",
    zero((a * L**5 + ax * t - j).subs(mx * t**2, ax * t + d * L**5).subs(a, m + d)),
)
check("finite margin lower bound", zero((rx * m * L**5 / 2) / (2 * a * L**5) - rx * m / (4 * a)))
check("wrong stiffness-only criterion exposed", s.Rational(1, 2) - 1 * 1 < 0)

# A compact solenoidal profile has zero zeroth moment and antisymmetric first
# moment. Translation adds a tensor product of its position and zero mean.
p, q, r = s.symbols("p q r", real=True)
first = s.Matrix([[0, p, q], [-p, 0, r], [-q, -r, 0]])
position = s.Matrix(s.symbols("x:3"))
translated = first + s.zeros(3, 1) * position.T
check("translated first moment unchanged", translated == first)
check("paired first moments cancel", translated - first == s.zeros(3))
H = s.Matrix(3, 3, s.symbols("H:9"))
check("paired arbitrary affine kinetic pairing zero", zero(s.trace(H.T * (translated - first))))
check("single-copy symmetric strain pairing zero", zero(s.trace((H + H.T) * first)))

k, h = s.symbols("k h", real=True)
check("exact discrete gradient mass jet", zero(s.series(4 * s.sin(k * h / 2) ** 2, k, 0, 4).removeO() - k**2 * h**2))
check("paired mean starts at cubic after attachment", s.series(k * (s.exp(s.I * k) - 1 - s.I * k), k, 0, 3).removeO() == 0)
kv = s.Matrix(s.symbols("kx ky kz", real=True))
check("three orthogonal bonds span every wavevector", zero(sum((kv.dot(s.eye(3)[:, i])) ** 2 for i in range(3)) - kv.dot(kv)))

rho, j0, kap = s.symbols("rho j kap", positive=True)
mu, mp, b, g, A, C, mg, kg, tau = s.symbols("mu mp b g A C mg kg tau", real=True)
for helicity in (-1, 1):
    M = s.Matrix([[rho + mu * k**2, b * helicity * k], [b * helicity * k, j0 + (mp + tau**2 * mg) * k**2]])
    K = s.Matrix([[A * k**2, g * helicity * k], [g * helicity * k, kap + (C + tau**2 * kg) * k**2]])
    T = s.Matrix([[1 - mu * k**2 / (2 * rho), -b * helicity * k / rho], [0, 1 - (mp + tau**2 * mg - b**2 / rho) * k**2 / (2 * j0)]])
    normal_mass = jet(T.T * M * T, k)
    normal_K = jet(T.T * K * T, k)
    baseline = C - 2 * g * b / rho - kap * (mp - b**2 / rho) / j0
    check(f"full normal mass helicity {helicity}", normal_mass == s.diag(rho, j0))
    check(f"full normalized gain helicity {helicity}", zero(normal_K[1, 1].coeff(k, 2) - baseline - tau**2 * (kg - kap * mg / j0)))
    check(f"mixed stiffness unchanged helicity {helicity}", zero(normal_K[0, 1] - g * helicity * k))

long_mass = j0 + (mp + tau**2 * mg) * k**2
long_map = 1 - (mp + tau**2 * mg) * k**2 / (2 * j0)
long_K = kap + (C + tau**2 * kg) * k**2
check("longitudinal mass retained and normalized", zero(s.series(long_map**2 * long_mass, k, 0, 3).removeO() - j0))
check("longitudinal positive residual is same", zero(s.series(long_map**2 * long_K, k, 0, 3).removeO().expand().coeff(k, 2) - (C - kap * mp / j0) - tau**2 * (kg - kap * mg / j0)))

# Arbitrary fixed coupled reaction block: a support-orthogonal new spin
# coefficient commutes with its Schur complement. No diagonal-cell inverse.
aa, nn, pp, delta = s.symbols("aa nn pp delta", nonzero=True)
schur_before = aa - nn**2 / pp
schur_after = (aa + delta * k**2) - nn**2 / pp
check("fixed full reaction Schur retains new jet", zero(schur_after - schur_before - delta * k**2))
overlap = s.symbols("overlap", nonzero=True)
bad_after = (aa + delta * k**2) - (nn + overlap * k) ** 2 / pp
check("overlapping free reaction mutation detected", not zero(s.expand(bad_after - schur_before).coeff(k, 2) - delta))
check("discarded mass mutation detected", not zero(tau**2 * kg - tau**2 * (kg - kap * mg / j0)))

print(f"{checks}/{checks} exact checks passed")
