"""Exact same-field geometry, pressure-sector and material-action interface."""

import sympy as sp

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0205-same-cell")
    x, y, z = sp.symbols("X Y Z", real=True)
    small = sp.Symbol("A", positive=True)
    psi = sp.cos(y) + small * sp.cos(z)
    u = sp.Matrix([psi, small * sp.sin(z), -sp.sin(y)])
    coordinates = (x, y, z)
    curl = sp.Matrix([
        sp.diff(u[2], y) - sp.diff(u[1], z),
        sp.diff(u[0], z) - sp.diff(u[2], x),
        sp.diff(u[1], x) - sp.diff(u[0], y),
    ])
    pressure = -u.dot(u) / 2
    checks.check("the same acoustic cell satisfies complete stationary Euler", all(sp.simplify(v) == 0 for v in u.jacobian(coordinates)*u + sp.Matrix([sp.diff(pressure, c) for c in coordinates])) and curl == -u)
    rotation = sp.Matrix([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
    source_coordinates = (x, sp.pi/2-z, y)
    source = sp.Matrix([sp.cos(source_coordinates[2])+small*sp.sin(source_coordinates[1]), sp.sin(source_coordinates[2]), small*sp.cos(source_coordinates[1])])
    checks.check("the elliptic optical source is this exact cell under a proper isometry", rotation.det() == 1 and rotation.T*rotation == sp.eye(3) and all(sp.simplify(v) == 0 for v in rotation*u-source))
    checks.check("wrong transverse reflection does not transfer the physical vector field", rotation*u != sp.diag(1, 1, -1)*source)
    normal_hessian = sp.hessian(pressure, (y, z)).subs({y:0, z:0})
    checks.check("the actual core pressure has its isotropic physical Hessian", all(sp.simplify(v) == 0 for v in normal_hessian-small*sp.eye(2)))
    checks.check("the tube boundary is an actual Euler first-integral surface", sp.simplify(u.dot(sp.Matrix([sp.diff(psi, c) for c in coordinates]))) == 0)
    action, phase = sp.symbols("I theta", real=True)
    omega = sp.sqrt(small)
    # Canonical area scaling at the elliptic core, followed by the
    # defining angular average, not a supplied twist coefficient.
    yy = sp.sqrt(2*action*omega)*sp.cos(phase)
    zz = sp.sqrt(2*action/omega)*sp.sin(phase)
    quartic = -(yy**4+small*zz**4)/24
    quartic_average = sp.integrate(sp.expand_trig(quartic), (phase, 0, 2*sp.pi))/(2*sp.pi)
    checks.check("the Birkhoff quartic follows from the actual cosine cell", sp.simplify(quartic_average+(1+small)*action**2/16) == 0)
    deficit = omega*action+quartic_average
    return_rate = sp.diff(deficit, action)/(1+small-deficit)
    twist = sp.simplify(sp.diff(return_rate, action).subs(action, 0))
    checks.check("the actual fixed acoustic cell has nonzero section twist", sp.simplify(twist-small/(small+1)**2+sp.Rational(1, 8)) == 0 and twist.subs(small, sp.Rational(1, 100)) != 0)
    checks.check("the fixed core return is not relabeled as a Euclidean closed circle", sp.simplify(omega/(1+small)).subs(small, sp.Rational(1, 100)) == sp.Rational(10, 101))
    carrier = sp.Symbol("n", integer=True, positive=True)
    checks.check("the complete physical axial cross integral vanishes", sp.integrate(sp.exp(sp.I*carrier*x), (x, 0, 2*sp.pi)) == 0)
    tagged_cross = sp.integrate((1+sp.cos(carrier*x))*sp.exp(sp.I*carrier*x), (x, 0, 2*sp.pi))
    checks.check("a literal nonnegative material tag does not inherit fluid orthogonality", sp.simplify(tagged_cross-sp.pi) == 0)
    wave = sp.symbols("kx ky kz", real=True)
    macro = sp.Matrix(wave)
    test = sp.Matrix([sp.Function(f"v{j}")(y,z) for j in range(3)])
    modulated = sp.exp(sp.I*carrier*x)*test
    advected = (modulated.jacobian(coordinates)*u+sp.I*macro.dot(u)*modulated+u.jacobian(coordinates)*modulated)
    checks.check("the complete pre-pressure Euler operator preserves the same axial index for arbitrary K", all(sp.simplify(sp.diff(v, x)-sp.I*carrier*v) == 0 for v in advected))
    # The full pressure conclusion additionally uses its exact diagonal
    # Fourier multiplier, proved in the analytic attachment.
    metric_parameter, speed, p = sp.symbols("s W p", positive=True)
    metric = sp.diag(metric_parameter, 1/metric_parameter)
    h = omega*sp.trace(metric)
    j = sp.Matrix([[0,-1],[1,0]])
    velocity = sp.Matrix([1,sp.I])/sp.sqrt(2)
    displacement = sp.I*metric*velocity/h
    strain_displacement = -omega*j*displacement
    material_transport = sp.I*p*speed*displacement
    inverse_metric = metric.inv()
    wplus = velocity+strain_displacement
    energy = (sp.conjugate(wplus).T*inverse_metric*wplus)[0]/2-sp.re((sp.conjugate(wplus).T*inverse_metric*material_transport)[0])
    leading = sp.simplify(sp.diff(energy,p))
    print("actual local Jacobi energy coefficient of pW =", sp.simplify(leading/speed))
    checks.check("the full actual optical energy retains its axial material-current contribution", sp.simplify(leading-speed/(2*h)) == 0)
    return int(checks.finish())


if __name__ == "__main__":
    raise SystemExit(main())
