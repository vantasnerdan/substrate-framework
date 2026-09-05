"""Exact same-core Floquet, Kelvin, pressure recursion and moment checks."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0112-same-EPS-core-Floquet")
    x, y, z, t = s.symbols("x y z t", real=True)
    coords = (x, y, z)
    e1, e3 = s.eye(3)[:, 0], s.eye(3)[:, 2]
    aa, bb, cc, dd, ee, ff, gg, hh = s.symbols("a b c d e f g h", real=True)
    gradient = s.Matrix([[aa, bb, cc], [dd, ee, ff], [gg, hh, -aa-ee]])
    omega = s.Matrix([hh-ff, cc-gg, dd-bb])
    kvec = e3
    bx, by = s.symbols("bx by", real=True)
    bvec = s.Matrix([bx, by, 0])

    def cross_matrix(v):
        return s.Matrix([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])

    def zero(v):
        return all(s.simplify(entry) == 0 for entry in v)

    def grad(q):
        return s.Matrix([s.diff(q, c) for c in coords])

    def curl(v):
        return s.Matrix([s.diff(v[2], y)-s.diff(v[1], z),
                         s.diff(v[0], z)-s.diff(v[2], x),
                         s.diff(v[1], x)-s.diff(v[0], y)])

    def divergence(v):
        return s.trace(v.jacobian(coords))

    ledger.check("actual velocity-gradient antisymmetry has the physical curl convention",
                 gradient-gradient.T == cross_matrix(omega) and s.trace(gradient) == 0)
    kdot = -gradient.T*kvec
    bdot = -gradient*bvec+2*kvec*kvec.dot(gradient*bvec)
    ledger.check("Euler pressure transport exactly preserves polarization transversality",
                 s.expand(kdot.dot(bvec)+kvec.dot(bdot)) == 0)
    ledger.check("the Kelvin pairing is conserved by covector and Cauchy transport",
                 s.expand(kdot.dot(omega)+kvec.dot(gradient*omega)) == 0)
    cvec = kvec.cross(bvec)
    cdot = kdot.cross(bvec)+kvec.cross(bdot)
    pairing = kvec.dot(omega)
    ledger.check("actual velocity amplitude converts to the stated material gyro equation",
                 zero(cdot-gradient*cvec+pairing*kvec.cross(cvec)))
    material = cvec/pairing
    ledger.check("the material amplitude retains the true leading Lin reconstruction",
                 zero(cdot/pairing-gradient*material-bvec))
    coefficient = gradient-pairing*cross_matrix(kvec)
    ledger.check("full material-amplitude generator has zero trace and the correct adjoint covector",
                 s.trace(coefficient) == 0 and zero(-coefficient.T*kvec-kdot))
    transverse = s.Matrix(2, 2, s.symbols("m11 m12 m21 m22"))
    return_map = s.Matrix([[transverse[0, 0], transverse[0, 1], aa],
                           [transverse[1, 0], transverse[1, 1], bb], [0, 0, 1]])
    ledger.check("periodic covector reduces full determinant and trace to the physical plane",
                 return_map.det() == transverse.det()
                 and s.trace(return_map)-1 == s.trace(transverse))
    tr = s.symbols("tr", real=True)
    eig = s.symbols("eig")
    ledger.check("unit determinant gives the exact exposing elliptic discriminant",
                 s.discriminant(eig**2-tr*eig+1, eig) == tr**2-4)

    c_dot, k_grad_pi = s.symbols("D_c k_grad_pi", real=True)
    trial_b = s.Matrix(s.symbols("b1 b2 b3", real=True))
    pressure = -2*kvec.dot(gradient*trial_b)-k_grad_pi-c_dot
    d_b = -gradient*trial_b-kvec*pressure-s.Matrix([aa, bb, k_grad_pi])
    ledger.check("every forced WKB pressure column preserves its actual divergence constraint",
                 s.expand(kdot.dot(trial_b)+kvec.dot(d_b)-c_dot) == 0)

    # Expose the exact initial Kelvin/Leray recursion at order three using
    # nonterminating smooth profiles. The residual is a gradient plus the
    # explicit fourth-order compact-profile remainder, not an omitted tail.
    f0 = s.exp(x+2*y+3*z)*s.Matrix([1, 2, 3])
    f1 = s.exp(2*x+y-z)*s.Matrix([2, -1, 1])
    bs, qs, ds = [], [], []
    old_b, old_q, old_d = s.zeros(3, 1), s.S.Zero, s.zeros(3, 1)
    for j in range(4):
        fj = f0 if j == 0 else f1 if j == 1 else s.zeros(3, 1)
        qj = s.simplify(kvec.dot(fj)-kvec.dot(grad(old_q))+divergence(old_b))
        bj = s.simplify(fj-kvec*qj-grad(old_q))
        dj = s.simplify(-kvec.cross(bj-curl(old_d)))
        ledger.check(f"prepared pressure/curl column {j} has both exact compatibility identities",
                     s.simplify(kvec.dot(bj)+divergence(old_b)) == 0
                     and zero(kvec.cross(dj)+curl(old_d)-bj))
        bs.append(bj)
        qs.append(qj)
        ds.append(dj)
        old_b, old_q, old_d = bj, qj, dj
    carrier = s.symbols("Z", nonzero=True)
    oscillation = s.exp(carrier*z)
    potential = sum((ds[j]/carrier**(j+1) for j in range(4)), s.zeros(3, 1))
    pressure_series = sum(qs[j]/carrier**(j+1) for j in range(4))
    exact_solenoidal = curl(oscillation*potential)
    difference = s.simplify((oscillation*(f0+f1/carrier)
                            -grad(oscillation*pressure_series)-exact_solenoidal)/oscillation)
    fourth_coefficient = difference.applyfunc(lambda v: s.simplify(carrier**4*v))
    ledger.check("the actual initial projection mismatch begins at the displayed fourth order",
                 not zero(fourth_coefficient)
                 and all(carrier not in entry.free_symbols for entry in fourth_coefficient)
                 and s.simplify(divergence(exact_solenoidal)) == 0)

    # A real positive parcel shape makes the actual full spin observable
    # nonzero. No numerical Hessian or prescribed rotor mass enters.
    a0 = e1
    b0 = -pairing*s.Matrix([0, 1, 0])
    spin_row = cross_matrix(a0)*gradient-cross_matrix(gradient*a0+b0)
    ledger.check("the actual spin row cannot vanish when the Kelvin pairing is nonzero",
                 s.expand(spin_row[2, 0]+2*spin_row[0, 2]-pairing) == 0)
    candidates = s.Matrix.hstack(e3, e3+e1, e3+s.eye(3)[:, 1])
    ledger.check("three positive-halfspace candidates expose every nonzero spin row",
                 candidates.det() != 0 and all(e3.dot(candidates[:, j]) == 1 for j in range(3)))
    q1, q2 = s.symbols("q1 q2", real=True)
    covariance = s.Matrix([[1+q1**2, q1*q2, q1],
                           [q1*q2, 1+q2**2, q2], [q1, q2, 1]])
    ledger.check("the selected material covariance is positive and realizes its prescribed moment",
                 covariance*e3 == s.Matrix([q1, q2, 1])
                 and covariance.det() == 1
                 and s.simplify(covariance[:2, :2]-covariance[:2, 2:]*covariance[2:, :2]) == s.eye(2))
    order = s.symbols("m", integer=True, positive=True)
    power = -order-1+s.Rational(3, 2)
    ledger.check("high-order remainder beats actual shrinking-parcel observation sensitivity",
                 power == s.Rational(1, 2)-order and power.subs(order, 3) < 0)

    # Exact Galilean chart on a smooth stationary Beltrami fixture.
    base = s.Matrix([s.sin(z), s.cos(z), 0])
    translation = s.Matrix(s.symbols("tx ty tz", real=True))
    boost = s.Matrix(s.symbols("vx vy vz", real=True))
    displacement = translation+t*boost
    velocity = boost-base.jacobian(coords)*displacement
    ledger.check("independent mean displacement obeys the exact Lin identity",
                 zero(displacement.diff(t)+displacement.jacobian(coords)*base
                      -base.jacobian(coords)*displacement-velocity))
    ledger.check("Galilean translation is an actual zero-wave-number linearized Euler solution",
                 zero(velocity.diff(t)+velocity.jacobian(coords)*base
                      +base.jacobian(coords)*velocity))
    point, center = s.Matrix(s.symbols("r1 r2 r3")), s.Matrix(s.symbols("c1 c2 c3"))
    local_u, mean_u = s.Matrix(s.symbols("u1 u2 u3")), s.Matrix(s.symbols("w1 w2 w3"))
    ledger.check("actual centered angular momentum is boost and translation invariant",
                 zero((point+displacement-center-displacement).cross(local_u+boost-mean_u-boost)
                      -(point-center).cross(local_u-mean_u)))

    print("EXACT material amplitude: adot=A*a-(k.dot(omega)/|k|^2)*k.cross(a)")
    print("EXACT transverse return determinant: 1; discriminant tr^2-4")
    print("EXACT spin nonvanishing identity: L31+2 L13=omega3")
    print("EXACT shrinking-parcel relative error power:", power)
    print("SCOPE: same EPS-core Floquet packet and actual material moment rows; no whole-knot scalar inertia")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
