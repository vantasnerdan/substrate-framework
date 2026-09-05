"""Exact leading compact-sector energy/KKS check, including quadrature returns."""

import sympy as s

from substrate_framework.verification import CheckLedger

ledger = CheckLedger("P251/0085 compact WKB")
z, frequency, delta = s.symbols("z k delta", real=True)
cosine, sine = s.cos(z), s.sin(z)
weight = 1 + delta * cosine
# omega=(cos z,sin z,0), curl omega=-omega.
omega = s.Matrix([cosine, sine, 0])
e1 = omega
e2 = s.Matrix([-sine, cosine, 0])
d1 = s.I * frequency * cosine
d2 = -s.I * frequency * sine
# Exact differential syzygy applied to weight(z) exp(i k x).
bb = d1**3 * d2**2 * weight
cc = (3 * d1**2 * d2**2 - 2 * d1**4) * weight + d1**3 * d2 * s.diff(weight, z)
aa = -(d1**2 * d2**3 + 6 * d2**3 - 14 * d1**2 * d2) * weight
aa -= (6 * d1 * d2**2 - 3 * d1**3) * s.diff(weight, z)
aa -= d1**2 * d2 * s.diff(weight, z, 2)
xi = s.expand(aa * e1 + bb * e2 + cc * s.Matrix([0, 0, 1]))
velocity = xi.cross(omega)
ledger.check("the explicit variable-Beltrami syzygy is divergence free exactly",
             s.trigsimp(s.I * frequency * xi[0] + s.diff(xi[2], z)) == 0)
ledger.check("its complete induced velocity is divergence free exactly",
             s.trigsimp(s.I * frequency * velocity[0] + s.diff(velocity[2], z)) == 0)
curl_velocity = s.Matrix([-s.diff(velocity[1], z),
                          s.diff(velocity[0], z) - s.I * frequency * velocity[2],
                          s.I * frequency * velocity[1]])
energy = s.expand((s.conjugate(velocity).dot(velocity)
                   + s.re(s.conjugate(velocity).dot(curl_velocity))) / 2)
energy10 = s.expand(energy).coeff(frequency, 10)
coefficient = -cosine**3 * sine * weight
predicted = coefficient**2 * (1 + sine**2) / 2


def circle_mean(expression):
    """Exact Laurent constant term; avoids numerical quadrature."""
    t = s.symbols("t", nonzero=True)
    expression = s.expand(expression).subs({s.cos(z): (t + 1 / t) / 2,
                                            s.sin(z): (t - 1 / t) / (2 * s.I)})
    return s.expand(expression).coeff(t, 0)


ledger.check("complete energy includes the derivative return and matches its integrated symbol",
             s.simplify(circle_mean(energy10 - predicted)) == 0)
ledger.check("the explicit leading energy is strictly positive at the symmetric envelope",
             circle_mean(energy10).subs(delta, 0) == s.Rational(13, 512))

# Rotate the input phase by -i so the leading fifth-order displacement is real.
shifted = s.expand(-s.I * xi)
principal = s.Matrix([s.expand(s.re(v)).coeff(frequency, 5) for v in shifted])
quadrature = s.Matrix([-s.expand(s.im(v)).coeff(frequency, 4) for v in shifted])
pairing9 = -omega.dot(principal.cross(quadrature))
predicted_pairing = coefficient**2 * cosine / 2  # -lambda/2 with lambda=-1.
ledger.check("the complete compact cosine/sine KKS coefficient has the derived sign",
             s.simplify(circle_mean(pairing9 - predicted_pairing)) == 0)
ledger.check("a non-symmetric envelope gives an exact nonzero KKS pair",
             circle_mean(pairing9) == 7 * delta / 256)
ledger.check("positive energy survives the same nonzero-pair envelope",
             circle_mean(energy10).subs(delta, s.Rational(1, 2)) > 0)

# Universal triad sum from only curl omega=lambda omega.
w = s.Matrix(s.symbols("w0:3", real=True))
lam = s.symbols("lambda", nonzero=True, real=True)
g = s.Matrix(3, 3, s.symbols("g0:9", real=True))
curl_omega = s.Matrix([g[2, 1] - g[1, 2],
                       g[0, 2] - g[2, 0], g[1, 0] - g[0, 1]])
symbols = []
cross_sum = 0
for axis in range(3):
    normal = s.eye(3)[:, axis]
    transverse = normal.cross(w)
    cross_sum += transverse.dot(g * normal)
    symbols.append(2 * transverse.dot(transverse) + normal.dot(w)**2
                   + 2 * transverse.dot(g * normal) / lam)
ledger.check("the triad derivative contraction has the exact curl sign",
             s.expand(cross_sum + w.dot(curl_omega)) == 0)
ledger.check("full corrected energy symbols sum to three times vorticity squared",
             s.expand(sum(symbols) - 3 * w.dot(w)
                      - 2 * w.dot(lam * w - curl_omega) / lam) == 0)
print(f"exact plane-witness energy coefficient={s.factor(circle_mean(energy10))}")
print(f"exact plane-witness KKS coefficient={s.factor(circle_mean(pairing9))}")
raise SystemExit(ledger.finish())
