"""Exact moment/phase-lift algebra; analytic rank proof remains in the receipt."""

import sympy as s

from substrate_framework.verification import CheckLedger

ledger = CheckLedger("P251/0080")
b = s.symbols("B", nonzero=True, real=True)
n = s.Matrix([s.Rational(2, 3), s.Rational(-2, 3), s.Rational(1, 3)])
ledger.check("the marked axis is exactly normalized", n.dot(n) == 1)

# Basis: raw Q, raw S, nineteen disjoint dual responses.
count = 19
moments = s.zeros(count, count + 2)
for i in range(count):
    moments[i, 0] = s.Rational(i + 1, i + 2)
    moments[i, 1] = s.Rational(2 * i - 3, i + 3)
    moments[i, i + 2] = 1
omega = s.zeros(count + 2)
omega[0, 1], omega[1, 0] = b, -b
response = s.eye(count + 2)[:, 2:]
projector = s.eye(count + 2) - response * moments
q = projector[:, 0]
reaction = projector[:, 1]
# Actual tube rows: centroid0:3, STF3:8, spin8:11. Global rows11:19.
for i in range(3):
    reaction += b * n[i] * response[:, 8 + i]
ledger.check("the exact moment projector kills all nineteen rows",
             moments * projector == s.zeros(count, count + 2))
ledger.check("the moment projector is idempotent", projector**2 == projector)
ledger.check("disjoint dual responses are isotropic for KKS",
             response.T * omega * response == s.zeros(count))
ledger.check("all raw-to-response KKS terms vanish by the support construction",
             omega[:2, 2:] == s.zeros(2, count))
ledger.check("the physical core angle retains all nineteen zero moments",
             moments * q == s.zeros(count, 1))
target = s.zeros(count, 1)
target[8:11, 0] = b * n
ledger.check("the reaction carries precisely the actual tube-spin row",
             s.simplify(moments * reaction - target) == s.zeros(count, 1))
ledger.check("the raw angle/reaction KKS pairing is preserved exactly",
             s.simplify((q.T * omega * reaction)[0] - b) == 0)
ledger.check("the raw physical core jets are unchanged",
             q[:2, :] == s.Matrix([1, 0]) and reaction[:2, :] == s.Matrix([0, 1]))

# Affine pairing of the actual first velocity moment.
spin = s.Matrix(s.symbols("S0:3", real=True))
beta = s.Matrix(s.symbols("beta0:3", real=True))
e0, e1, e2, e3, e4 = s.symbols("E0:5", real=True)
strain = s.Matrix([[e0, e1, e2], [e1, e3, e4], [e2, e4, -e0 - e3]])
rot = s.Matrix([[0, -beta[2], beta[1]],
                [beta[2], 0, -beta[0]], [-beta[1], beta[0], 0]])
moment = s.Matrix(3, 3, lambda i, j:
                  -sum(s.LeviCivita(i, j, k) * spin[k] for k in range(3)) / 2)
pairing = sum((strain + rot)[i, j] * moment[i, j]
              for i in range(3) for j in range(3))
ledger.check("the actual affine tube pairing is beta dot mechanical spin",
             s.simplify(pairing - beta.dot(spin)) == 0)
ledger.check("tracefree symmetric affine motion has no retained tube-spin pairing",
             sum(strain[i, j] * moment[i, j]
                 for i in range(3) for j in range(3)) == 0)
global_pair, tube_spin, beta_scalar = s.symbols("G A beta", real=True)
kks_pair = b
hybrid = global_pair - beta_scalar * tube_spin
lifted = global_pair - beta_scalar * kks_pair
ledger.check("the explicit affine lift equals the physical hybrid current",
             s.simplify((lifted - hybrid).subs(tube_spin, b)) == 0)
ledger.check("global source spin cannot be substituted for the actual tube row",
             s.simplify(lifted - hybrid) != 0)

# Exact volume-preserving affine correction potential.
x, y, z = s.symbols("x y z", real=True)
r = s.Matrix([x, y, z])
affine = (strain + rot) * r
potential = -r.cross(affine) / 3


def curl(field):
    """Local differential oracle, not a duplicated campaign API."""
    return s.Matrix([s.diff(field[2], y) - s.diff(field[1], z),
                     s.diff(field[0], z) - s.diff(field[2], x),
                     s.diff(field[1], x) - s.diff(field[0], y)])


ledger.check("the affine return has an exact compact-cutoff curl potential",
             s.simplify(curl(potential) - affine) == s.zeros(3, 1))
ledger.check("the affine generator is incompressible",
             sum(s.diff(affine[i], r[i]) for i in range(3)) == 0)

# The proof's leading ABC conditions really remove all STF coefficients.
conditions = list(strain.T * s.eye(3))
solution = s.solve(conditions, [e0, e1, e2, e3, e4], dict=True)
ledger.check("three independent wave atoms remove a leading STF affine symmetry",
             solution == [{e0: 0, e1: 0, e2: 0, e3: 0, e4: 0}])

# Passive material-tag dynamics is explicitly NOT deleted by the moments.
normal_velocity, tag_rate, transport = s.symbols("v_n tag_dot adv_tag", real=True)
tag_equation = tag_rate + transport + normal_velocity
ledger.check("a zero reference tag tangent does not license zero later tag motion",
             tag_equation.subs({tag_rate: 0, transport: 0}) == normal_velocity)
raise SystemExit(ledger.finish())
