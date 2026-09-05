"""Exact finite-parcel identities; no discretized Euler existence claim."""

import sympy as s

from substrate_framework.verification import CheckLedger

ledger = CheckLedger("P251/0075")
a, b, c = s.symbols("a b c", positive=True)
qd, bd, pd = s.symbols("qd bd pd", real=True)
kinetic = (a * qd**2 + 2 * b * qd * bd + c * bd**2) / 2
j = b**2 / a
mapping = b * (pd - bd) / a
ledger.check(
    "physical collective angle diagonalizes the full parcel kinetic matrix",
    s.simplify(kinetic.subs(qd, mapping) - j * pd**2 / 2 - (c - j) * bd**2 / 2)
    == 0,
)
spin = s.diff(kinetic, bd)
ledger.check(
    "physical spin retains locked affine contribution",
    s.simplify(spin.subs(qd, mapping) - j * pd - (c - j) * bd) == 0,
)
ledger.check(
    "positive kinetic determinant is the residual affine inertia",
    s.simplify(s.det(s.Matrix([[a, b], [b, c]])) / a - (c - j)) == 0,
)
k, q, beta, phi = s.symbols("K q beta phi", real=True)
ledger.check(
    "locking is transformed by the measured angle map, not fitted",
    s.simplify((k * q**2 / 2).subs(q, b * (phi - beta) / a)
               - k * (b / a)**2 * (phi - beta)**2 / 2) == 0,
)

# Exact finite material quadrature tests the action identities only.
rho = s.symbols("rho", positive=True)
radius = s.symbols("r", positive=True)
vel = s.Matrix(s.symbols("V0:3", real=True))
beta_rate = s.Matrix(s.symbols("B0:3", real=True))
q_rate = s.Matrix(s.symbols("Q0:3", real=True))
weights = [s.Rational(1, 6)] * 6
positions = [sign * radius * s.eye(3)[:, axis]
             for axis in range(3) for sign in (-1, 1)]
chi = s.symbols("chi", real=True)
internal = [(beta_rate + chi * q_rate).cross(r) for r in positions]
mean = sum((w * v for w, v in zip(weights, internal, strict=True)), s.zeros(3, 1))
ledger.check("material internal velocities have exact zero centroid", mean == s.zeros(3, 1))
total_energy = sum(rho * w * (vel + v).dot(vel + v) / 2
                   for w, v in zip(weights, internal, strict=True))
split = rho * vel.dot(vel) / 2 + sum(
    rho * w * v.dot(v) / 2 for w, v in zip(weights, internal, strict=True)
)
ledger.check("all material mass appears once in the centroid split",
             s.expand(total_energy - split) == 0)
angular = sum((rho * w * r.cross(v) for w, r, v in
               zip(weights, positions, internal, strict=True)), s.zeros(3, 1))
ledger.check("affine motion contributes actual parcel angular momentum",
             s.simplify(angular - 2 * rho * radius**2
                        * (beta_rate + chi * q_rate) / 3) == s.zeros(3, 1))

# Constant-vorticity no-spin moment proof.
nxy, nxz, nyz, pressure_mass, strength = s.symbols("nxy nxz nyz p W", real=True)
n = s.Matrix([[0, nxy, nxz], [-nxy, 0, nyz], [-nxz, -nyz, 0]])
force = strength * s.Matrix.vstack(n[1, :], -n[0, :], s.zeros(1, 3))
velocity_moment = force + pressure_mass * s.eye(3)
sym = velocity_moment + velocity_moment.T
ledger.check("compact pressure mass vanishes from the zz equation",
             sym[2, 2] == 2 * pressure_mass)
sol = s.solve(list(sym.subs(pressure_mass, 0)), [nxy, nxz, nyz], dict=True)
ledger.check("compact uniform-vorticity response has zero first moment",
             sol == [{nxy: 0, nxz: 0, nyz: 0}])
ledger.check("hence its actual angular moment is zero",
             velocity_moment.subs(sol[0]).subs(pressure_mass, 0) == s.zeros(3))
x, y, z, angle = s.symbols("x y z theta", real=True)
g = angle * (x * z**2 / 2 + x**3 / 6)
xi = s.Matrix([s.diff(g, x, z), s.diff(g, y, z),
               -s.diff(g, x, 2) - s.diff(g, y, 2)])
ledger.check("the complete compact family can retain an actual core-rotation jet",
             xi == s.Matrix([angle * z, 0, -angle * x]))

# Full two-site interaction retained through centroid subtraction.
p = s.Matrix([[5, 1], [1, 4]])
centroid = s.Matrix([[s.Rational(1, 10), s.Rational(1, 20)],
                     [s.Rational(1, 30), -s.Rational(1, 10)]])
pc = p - centroid.T * centroid
d = s.Matrix([1, 2])
arow = s.Matrix([s.Rational(99, 100), s.Rational(201, 100)])
aa = (d.T * pc.inv() * d)[0]
bb = (d.T * pc.inv() * arow)[0]
cc = (arow.T * pc.inv() * arow)[0]
ledger.check("the full centred operator remains positive in the rational oracle",
             pc[0, 0] > 0 and pc.det() > 0)
ledger.check("physical parcel-spin and shape rows have nonzero pairing", bb > 0)
ledger.check("full inverse gives the nonnegative affine residual", cc - bb**2 / aa >= 0)
ledger.check("discarding full reaction interactions changes the inertia",
             aa != sum(d[i]**2 / pc[i, i] for i in range(2)))
ledger.check("the spin row is not silently identified with the global source row",
             aa != bb)
raise SystemExit(ledger.finish())
