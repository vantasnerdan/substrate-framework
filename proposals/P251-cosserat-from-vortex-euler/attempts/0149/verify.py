"""Stationary square-pair Green symmetry, actual KKS and director action."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0149-stationary-pair")
    area, rho, circulation, length = s.symbols("A rho Gamma L", positive=True)
    gx, gy, gxx, gxy, gyy = s.symbols("gx gy gxx gxy gyy", real=True)
    quarter = s.Matrix([[0, -1], [1, 0]])
    offset = s.Matrix([length/2, length/2])
    checks.check("quarter turn fixes the antipode modulo the square period",
                 quarter*offset-offset == s.Matrix([-length, 0]))
    gradient = s.Matrix([gx, gy])
    stationary = s.solve(list(gradient+gradient), (gx, gy))
    checks.check("Green inversion symmetry forces exact point stationarity",
                 gradient.subs(stationary) == s.zeros(2, 1))
    hessian = s.Matrix([[gxx, gxy], [gxy, gyy]])
    conditions = list(quarter.T*hessian*quarter-hessian)
    conditions.append(s.trace(hessian)+1/area)
    values = s.solve(conditions, (gxx, gxy, gyy))
    actual = -rho*circulation**2*hessian.subs(values)
    checks.check("square symmetry and actual Poisson trace derive both positive Hessians",
                 actual == rho*circulation**2*s.eye(2)/(2*area))
    xx, yy, dx, dy = s.symbols("X Y rx ry", real=True)
    positions = s.Matrix([xx-dx/2, yy-dy/2, xx+dx/2, yy+dy/2])
    jacobian = positions.jacobian((xx, yy, dx, dy))
    skew = s.Matrix([[0, 1], [-1, 0]])
    kks = rho*circulation*s.diag(skew, skew)
    reduced = s.simplify(jacobian.T*kks*jacobian)
    checks.check("full point-Euler KKS splits into common and relative forms",
                 reduced == rho*circulation*s.diag(2*skew, skew/2))
    difference = s.Matrix([dx, dy]).jacobian((dx, dy))
    full_difference = s.Matrix([-s.eye(2), s.eye(2)]).T
    full_hessian = full_difference.T*actual*full_difference
    checks.check("full position Hessian has exactly the two translation zero modes",
                 full_hessian.rank() == 2
                 and full_hessian*s.Matrix([1, 0, 1, 0]) == s.zeros(4, 1)
                 and full_hessian*s.Matrix([0, 1, 0, 1]) == s.zeros(4, 1)
                 and difference == s.eye(2))

    radius, radial, angle, rate = s.symbols("d a theta theta_dot", real=True)
    b = rho*circulation*radius/2
    h = actual[0, 0]
    action = -b*radial*rate-h*(radial**2+radius**2*angle**2)/2
    reaction = s.solve(s.diff(action, radial), radial)[0]
    scalar = s.simplify(action.subs(radial, reaction))
    inertia = s.diff(scalar, rate, 2)
    stiffness = -s.diff(scalar, angle, 2)
    checks.check("radial reaction elimination derives positive geometric-angle inertia",
                 s.factor(inertia-rho*area*radius**2/2) == 0)
    checks.check("relative-angle locking is the actual periodic Green Hessian",
                 s.factor(stiffness-rho*circulation**2*radius**2/(2*area)) == 0)
    checks.check("actual pair optical frequency is fixed by action, not a fit",
                 s.factor(stiffness/inertia-circulation**2/area**2) == 0)
    momentum = s.diff(action, rate)
    checks.check("vorticity second-moment variation has the same initial spin row",
                 s.diff(-rho*circulation*(radius+radial)**2/4, radial).subs(radial, 0)
                 == s.diff(momentum, radial))
    tag_mass = s.symbols("M", positive=True)
    required = s.solve(tag_mass*radius**2/2-inertia, tag_mass)[0]
    checks.check("centroid-only tags expose missing ambient mass spin",
                 required == rho*area and 2*required > rho*area)

    center = s.symbols("Gamma_c", real=True)
    x = s.symbols("imbalance", real=True)
    central_energy = -3*rho*circulation*center*s.log(1-x*x)/(4*s.pi)
    ring_radial = 3*rho*circulation**2/(8*s.pi)
    joined = ring_radial+s.diff(central_energy, x, 2).subs(x, 0)
    critical_center = -s.Rational(5, 2)*circulation
    checks.check("stationarizing the dilute central polygon changes its radial sign",
                 s.factor(joined.subs(center, critical_center)
                          +27*rho*circulation**2/(8*s.pi)) == 0)
    # Remove the physical constant-vorticity term: the defining Poisson
    # trace is then zero and C4 symmetry falsely erases the positive mode.
    wrong_values = s.solve(conditions[:-1]+[s.trace(hessian)], (gxx, gxy, gyy))
    checks.check("omitting physical compensation is exposed by the mode Hessian",
                 hessian.subs(wrong_values) == s.zeros(2)
                 and actual != s.zeros(2))
    print("Scope: exact periodic point-Euler pair/action and analytic smooth-core bridge;")
    print("material-spin history, 3D coupled continuum and EPS join remain active")
    raise SystemExit(checks.finish())


if __name__ == "__main__":
    main()
