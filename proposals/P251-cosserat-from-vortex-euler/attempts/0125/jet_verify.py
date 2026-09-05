"""Exact algebraic joins for the actual Euler/Lin response in 0125.

The generic matrices below verify operator identities, not an invented
oscillator model; their actual Euler coefficients are derived in the note.
"""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0125-physical-jet-joins")
    t, eps, perturbation = s.symbols("t eps perturbation", real=True)
    observation = s.Matrix([[1, t, 0], [0, 1, t]])
    lifting = s.Matrix([[1, -t], [0, 1], [0, 0]])
    generator = s.Matrix(3, 3, s.symbols("l:9"))
    complement = s.eye(3)-lifting*observation
    measured_generator = observation.diff(t)+observation*generator
    retained = measured_generator*lifting
    source = generator*lifting-lifting.diff(t)-lifting*retained
    hidden = generator-lifting*measured_generator
    ledger.check("physical observation lifting is a genuine right inverse",
                 observation*lifting == s.eye(2))
    ledger.check("complement source remains in the moving physical null space",
                 s.simplify(observation*source) == s.zeros(2))
    ledger.check("null-space propagator respects the moving physical rows",
                 s.simplify((observation.diff(t)+observation*hidden)*complement)
                 == s.zeros(2, 3))
    x = s.Matrix(s.symbols("x:2"))
    y = complement*s.Matrix(s.symbols("y:3"))
    reconstructed = (lifting.diff(t)*x+lifting*(retained*x+measured_generator*y)
                     +source*x+hidden*y)
    ledger.check("physical retained plus hidden equations reconstruct the full operator",
                 s.simplify(reconstructed-generator*(lifting*x+y)) == s.zeros(3, 1))

    # A moving particle proves the current identity before integration. Keeping
    # every first variation tests centroid phases, material displacement and
    # material velocity simultaneously; linear integration gives the tag law.
    center, radius, velocity, dc, dr, dv = s.symbols("center radius velocity dc dr dv")
    moved_center = center+perturbation*dc
    moved_radius = radius+perturbation*dr
    moved_velocity = velocity+perturbation*dv
    difference = moved_velocity*(s.exp(-s.I*eps*(moved_center+moved_radius))
                                 -s.exp(-s.I*eps*moved_center))
    exact_response = s.series(s.diff(difference, perturbation).subs(perturbation, 0),
                             eps, 0, 3).removeO()
    moment_expression = s.exp(-s.I*eps*moved_center)*(
        -s.I*eps*moved_velocity*moved_radius
        -eps**2*moved_velocity*moved_radius**2/2)
    moment_response = s.series(s.diff(moment_expression, perturbation).subs(perturbation, 0),
                              eps, 0, 3).removeO()
    ledger.check("hybrid current quadrupole retains all material and centroid variations",
                 s.expand(exact_response-moment_response) == 0)
    ledger.check("quadrupole response is generally nonzero at second spatial order",
                 s.expand(exact_response).coeff(eps, 2) != 0)

    # Ordered prepared-transfer closure, with no commutativity assumption.
    transfer0, transfer1, transfer2, dot0, dot1, dot2 = s.symbols(
        "transfer0 transfer1 transfer2 dot0 dot1 dot2", commutative=False)
    g0 = dot0*transfer0**-1
    g1 = (dot1-g0*transfer1)*transfer0**-1
    g2 = (dot2-g0*transfer2-g1*transfer1)*transfer0**-1
    ledger.check("physical prepared-transfer first generator coefficient",
                 s.expand(g1*transfer0+g0*transfer1-dot1) == 0)
    ledger.check("physical prepared-transfer second generator coefficient",
                 s.expand(g2*transfer0+g1*transfer1+g0*transfer2-dot2) == 0)

    # Exact leading isotropic material-tag spin tensor identity.
    omega = s.Symbol("omega")
    d = s.Matrix(3, 3, s.symbols("d:9"))
    d[2, 2] = -d[0, 0]-d[1, 1]
    velocity_gradient = s.Matrix(3, 3, s.symbols("v:9"))
    axis = s.Matrix([0, 0, 1])
    rotation = s.Matrix([[0, -omega, 0], [omega, 0, 0], [0, 0, 0]])
    spin = s.zeros(3, 1)
    for j in range(3):
        unit = s.eye(3)[:, j]
        spin += d[:, j].cross(rotation[:, j])
        spin += unit.cross((rotation*d+velocity_gradient)[:, j])
    curl = s.Matrix([velocity_gradient[2, 1]-velocity_gradient[1, 2],
                     velocity_gradient[0, 2]-velocity_gradient[2, 0],
                     velocity_gradient[1, 0]-velocity_gradient[0, 1]])
    ledger.check("complete material spin is curl plus the retained shape-gradient row",
                 s.expand(spin-curl+omega*(d+d.T)*axis) == s.zeros(3, 1))
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
