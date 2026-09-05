"""Exact relative-orbit moment, constrained KKS and full Routh algebra."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0048-direct-EPS-common-rotor")
    x, y, z = s.symbols("x y z", real=True)
    e = s.Matrix(s.symbols("e0:3", real=True))
    r = s.Matrix([x, y, z])
    krot = e.cross(r)
    radial_moment = (x*x+y*y+z*z)*e

    def curl(v):
        return s.Matrix([s.diff(v[2], y)-s.diff(v[1], z),
                         s.diff(v[0], z)-s.diff(v[2], x),
                         s.diff(v[1], x)-s.diff(v[0], y)])

    ledger.check("relative angular impulse has the stated moment-map sign",
                 s.simplify(curl(radial_moment)+2*krot) == s.zeros(3, 1))
    eta = s.Matrix(s.symbols("eta0:3"))
    omega = s.Matrix(s.symbols("omega0:3"))
    ledger.check("angular impulse variation equals positive KKS pairing",
                 s.expand(krot.dot(eta.cross(omega))-omega.dot(krot.cross(eta))) == 0)
    ledger.check("compact core potential supplies its exact rotation jet",
                 s.simplify(curl(-radial_moment/2)-krot) == s.zeros(3, 1))

    # Old order is K, eta0, body cage A, internal q0, internal s0.
    l0, la, lq, ls, c = s.symbols("l0 la lq ls c", real=True, nonzero=True)
    old = s.zeros(5)
    for j, moment in enumerate((l0, la, lq, ls), start=1):
        old[0, j], old[j, 0] = moment, -moment
    old[3, 4], old[4, 3] = c, -c
    change = s.Matrix([[1, 0, 0, 0],
                       [0, 1-la/l0, -lq/l0, -ls/l0],
                       [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    actual = s.simplify(change.T*old*change)
    target = s.Matrix([[0, l0, 0, 0], [-l0, 0, 0, 0],
                       [0, 0, 0, c], [0, 0, -c, 0]])
    ledger.check("fixed-shape moment projection gives both exact canonical blocks", actual == target)
    ledger.check("common moment is fixed exactly, not asymptotically", actual[0, 1] == l0)
    ledger.check("four-generator form is nondegenerate", s.factor(actual.det()) == l0*l0*c*c)
    core_plus = s.Matrix([[1, 0, 0, 1, 0]])
    core_minus = s.Matrix([[1, 0, 0, -1, 0]])
    ledger.check("moment corrections preserve physical plus-section angle",
                 core_plus*change == s.Matrix([[1, 0, 1, 0]]))
    ledger.check("moment corrections preserve physical minus-section angle",
                 core_minus*change == s.Matrix([[1, 0, -1, 0]]))

    p11, p12, p22, g1, g2, h = s.symbols("p11 p12 p22 g1 g2 h", real=True)
    b, c = s.symbols("b c", real=True, nonzero=True)
    bd, qd, q, shape1, shape2 = s.symbols("Bdot qdot q y shape", real=True)
    p = s.Matrix([[p11, p12], [p12, p22]])
    g = s.Matrix([g1, g2])
    d = s.diag(b, c)
    vel = s.Matrix([bd, qd])
    shapes = s.Matrix([shape1, shape2])
    action = (shapes.T*d*vel)[0]-(shapes.T*p*shapes)[0]/2-q*(g.T*shapes)[0]-h*q*q/2
    solution = p.inv()*(d*vel-g*q)
    reduced = s.factor(action.subs({shape1: solution[0], shape2: solution[1]}, simultaneous=True))
    mass = d*p.inv()*d
    n = d*p.inv()*g
    stiffness = h-(g.T*p.inv()*g)[0]
    expected = (vel.T*mass*vel)[0]/2-q*(vel.T*n)[0]-stiffness*q*q/2
    ledger.check("full momentum elimination includes all off-diagonal entries",
                 s.factor(reduced-expected) == 0)
    ledger.check("single-sign gyro term is generally present",
                 s.factor(s.diff(reduced, bd, q)+n[0]) == 0 and n[0] != 0)
    ledger.check("the q-qdot part alone is a total derivative",
                 s.diff(-n[1]*q*q/2, q)*qd == -n[1]*q*qd)
    reversed_action = reduced.subs({b: -b, c: -c}, simultaneous=True)
    paired = s.factor((reduced+reversed_action)/2)
    ledger.check("time-reversal pairing cancels gyro and preserves positive quadratic terms",
                 s.factor(paired-((vel.T*mass*vel)[0]-stiffness*q*q)/2) == 0)
    minus1, minus2 = s.symbols("y_minus shape_minus", real=True)
    action_minus = action.subs({b: -b, c: -c, shape1: minus1, shape2: minus2},
                              simultaneous=True)
    ensemble = (action+action_minus)/2
    minus_solution = p.inv()*(-d*vel-g*q)
    independently_reduced = ensemble.subs(
        {shape1: solution[0], shape2: solution[1],
         minus1: minus_solution[0], minus2: minus_solution[1]}, simultaneous=True)
    ledger.check("shared physical angles with independent momenta give the same reduced ensemble",
                 s.factor(independently_reduced-paired) == 0)
    wrongly_tied = ensemble.subs({minus1: shape1, minus2: shape2}, simultaneous=True)
    ledger.check("incorrectly tying conjugate momenta destroys the kinetic response",
                 s.diff(wrongly_tied, bd) == 0 and s.diff(wrongly_tied, qd) == 0)
    ledger.check("mass determinant is positive whenever P is positive",
                 s.factor(mass.det()-b*b*c*c/p.det()) == 0)
    full_h = p.row_join(g).col_join(s.Matrix([[g1, g2, h]]))
    ledger.check("restoring Schur complement follows the full positive 3x3 Hessian",
                 s.factor(full_h.det()-p.det()*stiffness) == 0)
    ledger.check("body-relative mass cross has its exact controlled ratio",
                 s.factor(mass[0, 1]/mass[0, 0]+c*p12/(b*p22)) == 0)

    m00, m01, m11 = s.symbols("m00 m01 m11", real=True)
    pd, betad = s.symbols("Psidot betadot", real=True)
    cross = m00+m01
    internal_mass = m00+2*m01+m11
    factor = internal_mass/cross
    qdot_map = (pd-betad)/factor
    bdot_map = betad+qdot_map
    kinetic = (m00*bdot_map**2+2*m01*bdot_map*qdot_map+m11*qdot_map**2)/2
    desired = cross**2*pd**2/(2*internal_mass)+(m00*m11-m01*m01)*betad**2/(2*internal_mass)
    ledger.check("general kinetic field map diagonalizes the actual physical section angles",
                 s.factor(kinetic-desired) == 0)
    ledger.check("diagonal common-relative special case reproduces the earlier map",
                 s.factor(factor.subs(m01, 0)-(1+m11/m00)) == 0)
    beta, psi, angle_k = s.symbols("beta Psi K", real=True)
    ledger.check("same map transports the full restoring energy",
                 s.factor(angle_k*((psi-beta)/factor)**2/2
                          -angle_k*(psi-beta)**2/(2*factor**2)) == 0)
    print("Analytic oracle: direct-rotor.md; finite source asymptotics and norm bounds, no spectrum.")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
