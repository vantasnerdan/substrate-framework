"""Exact affine Euler, physical annular-clock and joint-current checks."""

import sympy as s

from substrate_framework.euler_phase import physical_configuration_chart
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0164-all-time-physical-joint-jet")
    x, y, zc, t = s.symbols("x y zc t", real=True)
    xyz = s.Matrix([x, y, zc])
    u = s.Matrix([s.sin(y), -s.sin(x), s.cos(x)+s.cos(y)])
    rot = s.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]])

    def curl(vector):
        return s.Matrix([s.diff(vector[2], y)-s.diff(vector[1], zc),
                         s.diff(vector[0], zc)-s.diff(vector[2], x),
                         s.diff(vector[1], x)-s.diff(vector[0], y)])

    def zero(matrix):
        return matrix.applyfunc(s.simplify) == s.zeros(*matrix.shape)

    jac = u.jacobian(xyz)
    pressure = -s.cos(x)*s.cos(y)
    checks.check("the exposing periodic 2D3C field is actual stationary Euler",
                 zero(jac*u+s.Matrix([s.diff(pressure, q) for q in xyz])))
    checks.check("the full rotational Coriolis field is a pressure gradient",
                 zero(rot*u-s.Matrix([s.diff(-s.cos(x)-s.cos(y), q) for q in xyz]))
                 and zero(curl(rot*u)))
    xi = rot*xyz
    static_velocity = rot*u-jac*xi
    checks.check("the static rotation is the actual Kelvin affine column",
                 zero(xi.cross(curl(u))-static_velocity
                      -s.Matrix([s.diff(u.dot(xi), q) for q in xyz])))
    eta = t*xi
    velocity = xi+t*static_velocity
    checks.check("the common-V affine column obeys full material Lin transport",
                 zero(eta.diff(t)+eta.jacobian(xyz)*u-jac*eta-velocity))
    acceleration = velocity.diff(t)+velocity.jacobian(xyz)*u+jac*velocity
    checks.check("the rotating affine column solves Euler with its full pressure",
                 zero(curl(acceleration)))
    checks.check("omitting the microscopic rotating-cell velocity is exposed",
                 not zero(curl(xi.jacobian(xyz)*u+jac*xi)))

    # Non-axis-independent exact Beltrami extension: lambda=-1.
    wave = s.Matrix([s.cos(zc), s.sin(zc), 0])
    derivative = wave.diff(zc)
    checks.check("the Beltrami affine correction uses the actual signed curl eigenvalue",
                 zero(curl(wave)+wave) and zero(rot*wave-derivative))
    corrected_velocity = xi-t*derivative
    corrected_displacement = t*xi-t*t*derivative
    wave_jac = wave.jacobian(xyz)
    checks.check("the nonzero-axis-derivative affine velocity solves full Euler exactly",
                 zero(corrected_velocity.diff(t)+corrected_velocity.jacobian(xyz)*wave
                      +wave_jac*corrected_velocity))
    checks.check("the local transport correction solves actual Lin displacement",
                 zero(corrected_displacement.diff(t)
                      +corrected_displacement.jacobian(xyz)*wave
                      -wave_jac*corrected_displacement-corrected_velocity))
    checks.check("dropping the Beltrami material correction is exposed",
                 not zero((t*xi).diff(t)+(t*xi).jacobian(xyz)*wave
                          -wave_jac*(t*xi)-corrected_velocity))

    kx, ky, kz = s.symbols("kx ky kz", real=True)
    kv = s.Matrix([kx, ky, kz])

    def cross_matrix(v):
        return s.Matrix([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])

    mirror = s.diag(-1, 1, 1)
    checks.check("the real curl symbol has the correct whole-field axial parity",
                 zero(cross_matrix(mirror*kv)
                      -mirror.det()*mirror*cross_matrix(kv)*mirror.T))
    checks.check("the chiral scalar second-jet term is removed by mirror pairing",
                 s.eye(3)+mirror.det()*mirror*s.eye(3)*mirror.T == s.zeros(3))

    h, node, c1, c2, c3, p = s.symbols("h node c1 c2 c3 p", positive=True)
    xx = s.Symbol("xx", real=True)
    nodal = c1*(xx-node)+c2*(xx-node)**2+c3*(xx-node)**3
    xp, xm = node+h, node-h
    weights = (s.Integer(2), s.Integer(1))
    points = (xp, xm)
    fs = [nodal.subs(xx, q) for q in points]
    fs1 = [s.Rational(3, 2)*q*s.diff(nodal, xx).subs(xx, q) for q in points]
    fs2 = [(s.Rational(9, 4)*q*q*s.diff(nodal, xx, 2)
            +s.Rational(3, 4)*q*s.diff(nodal, xx)).subs(xx, q) for q in points]
    denominator = sum(w*f for w, f in zip(weights, fs))
    numerator = sum(w*q*f for w, q, f in zip(weights, points, fs))
    den1 = sum(w*f for w, f in zip(weights, fs1))
    den2 = sum(w*f for w, f in zip(weights, fs2))
    num1 = sum(w*q*f for w, q, f in zip(weights, points, fs1))
    num2 = sum(w*q*f for w, q, f in zip(weights, points, fs2))
    mean = s.cancel(numerator/denominator)
    mean1 = s.cancel(num1/denominator-numerator*den1/denominator**2)
    mean2 = s.cancel(num2/denominator-numerator*den2/denominator**2
                     -2*num1*den1/denominator**2+2*numerator*den1**2/denominator**3)
    checks.check("fixed-radius nodal annuli derive their finite first carrier response",
                 s.limit(mean1, h, 0) == -12*node)
    checks.check("the second carrier quotient has the positive nodal susceptibility",
                 s.limit(h*mean2, h, 0) == 108*node**2)
    second = sum(w*q*q*f for w, q, f in zip(weights, points, fs))/denominator
    checks.check("positive physical annuli have the actual signed small clock variance",
                 s.limit(s.cancel((second-mean**2)/h**2), h, 0) == -8)
    laguerre = s.assoc_laguerre(8, 1, xx)
    checks.check("the selected actual Laguerre mode has simple nodal roots",
                 s.gcd(laguerre, s.diff(laguerre, xx)) == 1)
    om, cd = s.symbols("Omega cD", positive=True)
    sigma = -2*om+18*cd*p**s.Rational(-1, 2)
    ring = sigma-cd*p*node/2
    checks.check("a fixed single annulus really has negative leading curvature",
                 s.diff(s.diff(ring**2, p, 2).subs(p, 1), cd).subs(cd, 0) == -54*om)
    gamma2_leading = -4*om*(s.Rational(27, 2)-mean1-mean2/2)
    checks.check("two physical annuli retain positive actual squared-clock curvature",
                 s.limit(h*gamma2_leading, h, 0) == 216*om*node**2)

    rho, j, inertia, k, acoustic, optical, freq = s.symbols(
        "rho j I k a B freq", positive=True)
    measured_map = s.Matrix([[1-inertia*k*k/(4*rho), -j*k/(2*rho)], [k/2, 1]])
    inverse = measured_map.inv()
    mass0 = s.diag(rho, j)
    stiffness0 = s.diag(rho*acoustic*k*k, j*(om*om+optical*k*k))
    mass = inverse.T*mass0*inverse
    stiffness = inverse.T*stiffness0*inverse

    def jet(matrix):
        return matrix.applyfunc(lambda value: s.series(value, k, 0, 3).removeO().expand())

    checks.check("the true acoustic polar moment remains in translational gradient inertia",
                 zero(jet(mass)-s.Matrix([[rho+(inertia/2-j/4)*k*k, 0],
                                         [0, j-j*j*k*k/(4*rho)]])))
    expected_stiffness = s.Matrix([
        [(rho*acoustic+j*om*om/4)*k*k, -j*om*om*k/2],
        [-j*om*om*k/2, j*om*om+j*(optical-j*om*om/(2*rho))*k*k],
    ])
    checks.check("the full physical potential retains the leading locking coupling",
                 zero(jet(stiffness)-expected_stiffness))
    checks.check("both measured pullbacks preserve the entire dispersion determinant",
                 s.factor((stiffness-freq*mass).det()*measured_map.det()**2
                          -(stiffness0-freq*mass0).det()) == 0)

    jj = s.Matrix([[0, 1], [-1, 0]])
    symplectic = s.diag(rho*jj, j*jj)
    generator = s.diag(s.Matrix([[0, 1], [0, 0]]), s.Matrix([[0, 1], [-om*om, 0]]))
    eps = s.Symbol("eps", real=True)
    rows = s.Matrix([[1-inertia*k*k/(4*rho), 0, -j*k/(2*rho), -j*k*eps/(2*rho*om)],
                     [k/2, 0, 1, 0]])
    measured = s.Matrix.vstack(rho*(rows*generator)[:1, :], s.Matrix([[0, inertia*k/2, 0, j]]))
    exact_bracket = -rows*symplectic.inv()*rows.T
    checks.check("the actual current phase component exposes its precise configuration remainder",
                 s.simplify(exact_bracket[0, 1]-eps*k/(2*rho*om)) == 0)
    # The general physical matrices and bracket above remain fully symbolic.
    # Use an exact rational exposing instance for the infrastructure consumer,
    # avoiding expression swell in its second generic symplectic inversion.
    instance = {rho: 2, j: 3, inertia: 5, om: 7, k: s.Rational(1, 11)}
    chart = physical_configuration_chart(
        symplectic.subs(instance), generator.subs(instance), rows.subs(instance),
        configuration_rate=s.zeros(2, 4), configuration_acceleration=s.zeros(2, 4),
        generator_rate=s.zeros(4), momentum=measured.subs(instance),
    )
    checks.check("the configuration API reproduces the independently derived physical bracket",
                 zero(chart.configuration_bracket-exact_bracket.subs(instance)))
    checks.check("the ideal current is a true configuration map, not a changed momentum name",
                 chart.configuration_bracket.subs(eps, 0) == s.zeros(2)
                 and chart.momentum_difference != s.zeros(2, 4))
    print("Derived nodal susceptibility: lim h*mu_PP =", s.limit(h*mean2, h, 0))
    print("Derived physical mass second jet:", jet(mass))
    print("Derived hybrid configuration remainder:", exact_bracket[0, 1])
    print("Scope: affine Euler/Bloch identity, O(3) parity, fixed-radius material clock,")
    print("and full measured current pullback; no global EPS or all-k closure is inferred.")
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
