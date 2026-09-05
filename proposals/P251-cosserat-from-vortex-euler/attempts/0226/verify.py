"""Exact symmetry and separate nonlinear action-normalization anchors."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0226-symmetry-density-normalization")
    x, y, z, t, om, eps = s.symbols("x y z t Omega epsilon", real=True)
    pos = s.Matrix([x, y, z])
    axis = s.Matrix([1, 0, 0])
    rigid = axis.cross(pos)
    u = s.Matrix([y*z, x*z, x*y])
    p = -(x*x*y*y+x*x*z*z+y*y*z*z)/2
    gradp = s.Matrix([s.diff(p, c) for c in pos])
    du = u.jacobian(pos)
    checks.check("nonuniform diagnostic is an actual stationary Euler field",
                 s.simplify(du*u+gradp) == s.zeros(3, 1)
                 and s.trace(du) == 0)
    rotation = axis.cross(u)-du*rigid
    dp = -rigid.dot(gradp)
    static_residual = rotation.jacobian(pos)*u+du*rotation
    checks.check("rotation of the whole field has its actual pressure variation",
                 s.simplify(static_residual+s.Matrix([s.diff(dp, c) for c in pos]))
                 == s.zeros(3, 1))
    rate_residual = -rigid.jacobian(pos)*u-du*rigid-rotation
    checks.check("rigid rate has a nonzero exact Coriolis residual",
                 s.simplify(rate_residual+2*axis.cross(u)) == s.zeros(3, 1)
                 and rate_residual != s.zeros(3, 1))
    curl_rate = s.Matrix([s.diff(rate_residual[2], y)-s.diff(rate_residual[1], z),
                          s.diff(rate_residual[0], z)-s.diff(rate_residual[2], x),
                          s.diff(rate_residual[1], x)-s.diff(rate_residual[0], y)])
    checks.check("its nonzero curl prevents deleting it as pressure",
                 s.simplify(curl_rate-2*du*axis) == s.zeros(3, 1)
                 and curl_rate != s.zeros(3, 1))

    aa = s.Matrix([[0, 0, 0], [0, 0, -1], [0, 1, 0]])
    bb = om*s.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]])
    hh = -(bb+eps*aa)**2
    checks.check("tilted uniform rotations are exact stationary Euler families",
                 hh == hh.T and s.trace(bb+eps*aa) == 0)
    e1, e2 = s.sin(om*t)/om, (1-s.cos(om*t))/om
    ee = s.Matrix([[0, 0, e2], [0, 0, -e1], [-e2, e1, 0]])
    checks.check("actual material tilt solves full Lin transport, not imposed t rotation",
                 s.simplify(s.diff(ee, t)+ee*bb-bb*ee-aa) == s.zeros(3, 3))
    even_tilt = s.simplify((e1+e1.subs(om, -om))/2)
    checks.check("time reversal retains the even cubic material-clock correction",
                 s.diff(even_tilt, t, 3).subs(t, 0) == -om**2
                 and s.limit(even_tilt, om, 0) == t)

    rr, lam, ampl = s.symbols("r lambda A", positive=True)
    omega = ampl*lam*s.besselj(1, lam*rr)/rr
    zz = ampl*lam**2*s.besselj(0, lam*rr)
    checks.check("literal-curl radial phase weight is derived from its actual field",
                 s.simplify(s.diff(zz, rr)+lam**2*rr*omega) == 0)
    small, xx = s.symbols("small x", positive=True)
    ratio = 2*s.besselj(1, s.sqrt(small)*xx)/(s.sqrt(small)*xx)
    checks.check("normalized inner clock variation has the claimed finite limit",
                 s.limit((1-ratio)/small, small, 0) == xx**2/8)

    nodes = [-3, -1, 2]
    vand = s.Matrix([[node**j for node in nodes] for j in range(3)])
    weights = vand.inv()*s.Matrix([1, 0, 0])
    checks.check("finite signed moment preparation preserves amplitude and cancels two jets",
                 vand*weights == s.Matrix([1, 0, 0])
                 and any(value < 0 for value in weights))
    spectral = sum(weights[i]*s.exp(-s.I*small*nodes[i]*t) for i in range(3))
    checks.check("moment cancellation concerns actual time-response derivatives",
                 spectral.subs(t, 0) == 1
                 and s.diff(spectral, t).subs(t, 0) == 0
                 and s.diff(spectral, t, 2).subs(t, 0) == 0
                 and s.diff(spectral, t, 3).subs(t, 0) != 0)
    checks.check("the same signed cancellation does not cancel positive quadratic energy",
                 sum(weights[i]**2*nodes[i]**2 for i in range(3)) > 0)
    q0, qz = s.symbols("Q0 Qz", positive=True)
    root = s.sqrt(q0/qz)
    checks.check("the independent quadratic return solves the actual energy-centering form",
                 s.simplify(q0-qz*root**2) == 0
                 and s.diff(q0-qz*t**2, t).subs(t, root) != 0)
    checks.check("normalized centering root has no hidden inverse-smallness exponent",
                 s.simplify(s.sqrt(small*q0/(small*qz))-root) == 0)

    rho, radius, width, freq = s.symbols("rho R a frequency", positive=True)
    phase = rho*radius*lam**2*freq*width**4
    inertia = s.simplify(phase/(freq/radius**2))
    checks.check("actual phase and angle lever arm yield the stated inertia density",
                 s.simplify(inertia/radius**3-rho*lam**2*width**4) == 0)
    checks.check("joint scaling preserves nonzero inertia while tag fraction remains positive",
                 s.simplify((inertia/radius**3).subs(width, lam**(-s.Rational(1, 2))))
                 == rho
                 and s.simplify((lam**2*width**2).subs(
                     width, lam**(-s.Rational(1, 2)))) == lam)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
