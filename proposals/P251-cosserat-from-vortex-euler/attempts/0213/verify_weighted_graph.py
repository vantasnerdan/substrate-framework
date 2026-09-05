"""Exact flat-tail reconstruction identities, with a pressure-term mutation."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0213-weighted-graph")
    r, k, c = s.symbols("r k c", positive=True)
    w, b, om = (s.Function(n)(r) for n in ("w", "b", "Omega"))
    er, et, ez = (s.Function(n)(r) for n in ("eta_r", "eta_theta", "eta_z"))
    d = c+w
    kap = -2*w*s.diff(w, r)/r
    zeta = kap/(2*om)
    y = b/d
    sturm = (s.diff(d**2*s.diff(y, r)/r, r)
             +(kap-k**2*d**2)*y/r-zeta*et-d*(s.diff(ez, r)-s.I*k*er))
    coefficient = s.diff(d, r, 2)/d-s.diff(d, r)/(r*d)+k**2-kap/d**2
    transformed = (s.diff(b, r, 2)-s.diff(b, r)/r-coefficient*b
                   -r*zeta*et/d-r*(s.diff(ez, r)-s.I*k*er))
    checks.check("weighted displacement normal form follows from the complete forced Sturm equation",
                 s.simplify(r*sturm/d-transformed) == 0)
    checks.check("the transformed first-derivative coefficient is exactly geometric minus one over r",
                 s.simplify(s.diff(transformed, s.diff(b, r))) == -1/r)
    checks.check("deleting radial d-curvature changes the real response equation",
                 s.simplify(r*sturm/d-(transformed+s.diff(d, r, 2)*b/d)) != 0)
    phi, delta = s.symbols("phi delta", positive=True)
    p1, p2 = s.symbols("phi_prime phi_second", real=True, nonzero=True)
    leading = phi**2*s.exp(-delta/phi)
    ratio = s.simplify(leading*s.diff(leading, phi, 2)/s.diff(leading, phi)**2)
    checks.check("new primitive tail has the required phi-derivative ratio",
                 s.limit(ratio, phi, 0, dir="+") == 1)
    wr = s.diff(leading, phi)*p1
    wrr = s.diff(leading, phi, 2)*p1**2+s.diff(leading, phi)*p2
    checks.check("the true radial chain rule preserves the tail Hardy limit",
                 s.limit(s.simplify(leading*wrr/wr**2), phi, 0, dir="+") == 1)
    checks.check("the primitive and source derivative are not silently identified",
                 s.limit(s.simplify(leading/s.diff(leading, phi)/phi**2),
                         phi, 0, dir="+") == 1/delta)
    xx, yy, zz = s.symbols("x y z", real=True)
    position = s.Matrix([xx, yy, zz])
    omega = s.Matrix(s.symbols("w1 w2 w3"))
    axis = s.Matrix(s.symbols("a1 a2 a3"))
    vel = omega.cross(position)/(position.dot(position))**s.Rational(3, 2)
    covariance = (axis.cross(vel)-vel.jacobian((xx, yy, zz))*axis.cross(position)
                  -(axis.cross(omega)).cross(position)
                  /(position.dot(position))**s.Rational(3, 2))
    checks.check("full Biot-Savart kernel rigid covariance supplies the commutator subtraction",
                 s.simplify(covariance) == s.zeros(3, 1))
    generator = s.Matrix([y/r, s.I*(et-2*om*y/r)/(k*d), s.I*s.diff(y, r)/(k*r)])
    checks.check("complete response generator is exactly solenoidal",
                 s.simplify(s.diff(r*generator[0], r)/r+s.I*k*generator[2]) == 0)
    force = generator.cross(s.Matrix([0, s.diff(w, r), zeta]))
    pressure = s.I*(d*s.diff(y, r)/r-ez)/k
    velocity = s.Matrix([-s.I*k*d*y/r-er, -zeta*y/r, s.diff(d*y, r)/r-ez])
    defect = force-s.Matrix([s.diff(pressure, r), 0, s.I*k*pressure])-velocity
    checks.check("azimuthal and axial coadjoint response rows retain the whole pressure",
                 s.simplify(defect[1]) == 0 and s.simplify(defect[2]) == 0)
    checks.check("radial coadjoint response follows from the actual forced Sturm equation",
                 s.simplify(defect[0]+s.I*sturm/(k*d)) == 0)
    checks.check("multiplying the angular generator by physical Z removes its false inverse-c estimate",
                 s.simplify(force[0]-(s.I*zeta*et/(k*d)
                                      -s.I*2*om*zeta*y/(k*r*d)
                                      -s.I*s.diff(w, r)*s.diff(y, r)/(k*r))) == 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
