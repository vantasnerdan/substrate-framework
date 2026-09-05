"""Exact reaction and moment algebra for compact-action.md."""

import sympy as s

checks = 0


def check(name, condition):
    global checks
    if not condition:
        raise AssertionError(name)
    checks += 1
    print(f"PASS {name}")


def zero(expression):
    return s.simplify(expression) == 0


P, B, H, N = s.symbols("P B H N", nonzero=True, real=True)
z, a, adot, e, C, F, A = s.symbols("z a adot e C F A", real=True)
lagrangian = B * z * adot - (P * z**2 + 2 * z * (N * a + C * e) + H * a**2 + 2 * F * a * e + A * e**2) / 2
reaction = (B * adot - N * a - C * e) / P
reduced = s.expand(lagrangian.subs(z, reaction))
expected = B**2 * adot**2 / (2 * P) - B * adot * (N * a + C * e) / P - (H * a**2 + 2 * F * a * e + A * e**2 - (N * a + C * e) ** 2 / P) / 2
check("full reaction stationary equation", zero(s.diff(lagrangian, z).subs(z, reaction)))
check("full Routh square including macro forcing", zero(reduced - expected))
paired = s.expand((reduced + reduced.subs(B, -B)) / 2)
check("time reversal cancels gyro after reaction", zero(paired.coeff(adot, 1)))
check("time reversal retains positive inertia", zero(2 * paired.coeff(adot, 2) - B**2 / P))
check("macro forcing Schur term retained", zero(2 * paired.coeff(e, 2) - (C**2 / P - A)))
check("reduced cage ratio is determinant over B squared", zero((H - N**2 / P) / (B**2 / P) - (H * P - N**2) / B**2))

t, edot = s.symbols("t edot", real=True)
attached = s.expand(paired.subs({a: t * e, adot: t * edot}))
check("STF leading full Schur stiffness", zero(-2 * attached.coeff(e, 2).coeff(t, 2) - (H - N**2 / P)))
check("STF linear macro correction", zero(-attached.coeff(e, 2).coeff(t, 1) - (F - N * C / P)))
check("STF actual gradient inertia", zero(2 * attached.coeff(edot, 2).coeff(t, 2) - B**2 / P))

# A genuinely coupled reaction block: never replace its inverse by diagonal
# entries or invert averaged cell coefficients.
p, q, r = s.symbols("p q r", real=True)
matrix = s.Matrix([[p, r], [r, q]])
f = s.Matrix(s.symbols("f:2"))
sol = matrix.inv() * f
check("complete two-reaction inverse", s.simplify(matrix * sol - f) == s.zeros(2, 1))
check("offdiagonal reaction contributes", not zero(matrix.inv()[0, 0] - 1 / p))
check("averaging isolated scalar coefficients mutation", s.Rational(1, 2) * (s.Rational(1, 1) + s.Rational(1, 4)) != 1 / s.Rational(5, 2))

# Response supports have zero pairwise KKS, so arbitrary spin projections
# leave the raw nonzero canonical pairing intact.
raw_B = s.symbols("raw_B", nonzero=True)
omega = s.diag(s.Matrix([[0, raw_B], [-raw_B, 0]]), s.zeros(3))
lq = s.Matrix(s.symbols("lq:3"))
ls = s.Matrix(s.symbols("ls:3"))
Q = s.Matrix([1, 0, *(-lq)])
S = s.Matrix([0, 1, *(-ls)])
spin = s.Matrix.hstack(lq, ls, s.eye(3))
check("both exact mechanical spin projections", spin * Q == spin * S == s.zeros(3, 1))
check("spin projection preserves exact KKS", zero((Q.T * omega * S)[0] - raw_B))

h0, k, b0 = s.symbols("h0 k b0", positive=True)
check("finite high-carrier gap lower bound", zero((h0 / 2) ** 2 / (2 * b0 / k) ** 2 - h0**2 * k**2 / (16 * b0**2)))
kk, hh = s.symbols("kk hh", real=True)
check("exact gradient difference square", zero(s.series(4 * s.sin(kk * hh / 2) ** 2, kk, 0, 3).removeO() - kk**2 * hh**2))
check("zero velocity moments imply cubic attached mean", s.series(kk * (s.exp(s.I * kk) - 1 - s.I * kk), kk, 0, 3).removeO() == 0)

kap, j, kg, jg, rho, b, g, cp, mp = s.symbols("kap j kg jg rho b g cp mp", nonzero=True)
normal = cp - 2 * g * b / rho - kap * (mp - b**2 / rho) / j
new = normal.subs({cp: cp + t**2 * kg, mp: mp + t**2 * jg})
check("full transverse normalized gain", zero(new - normal - t**2 * (kg - kap * jg / j)))
longitudinal = cp - kap * mp / j
check("full longitudinal normalized gain", zero(longitudinal.subs({cp: cp + t**2 * kg, mp: mp + t**2 * jg}) - longitudinal - t**2 * (kg - kap * jg / j)))
check("omitted reaction cross mutation detected", not zero(H - (H - N**2 / P)))
check("omitted gradient mass mutation detected", not zero(kg - (kg - kap * jg / j)))

phi = s.symbols("phi", real=True)
base_spin = B * reaction.subs({adot: phi, C: 0})
reversed_spin = -B * reaction.subs({B: -B, adot: phi, C: 0})
check("actual time-paired base spin equals reduced momentum", zero((base_spin + reversed_spin) / 2 - B**2 * phi / P))

print(f"{checks}/{checks} exact checks passed")
