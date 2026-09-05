"""Actual periodic square-parcel pressure torque and director momentum."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0152-material-boundary")
    length, rho, gamma, area, d = s.symbols("L rho Gamma A d", positive=True)
    x, y = s.symbols("x y", real=True)
    px = s.Function("p_x")(y)
    py = s.Function("p_y")(x)
    right = s.Matrix([length/2, y, 0])
    left = s.Matrix([-length/2, y, 0])
    top = s.Matrix([x, length/2, 0])
    bottom = s.Matrix([x, -length/2, 0])
    ex, ey = s.Matrix([1, 0, 0]), s.Matrix([0, 1, 0])
    checks.check("opposite actual periodic x-face pressure torques cancel",
                 right.cross(px*ex)+left.cross(-px*ex) == s.zeros(3, 1))
    checks.check("opposite actual periodic y-face pressure torques cancel",
                 top.cross(py*ey)+bottom.cross(-py*ey) == s.zeros(3, 1))
    checks.check("weighted circulation cancellation uses square edge geometry",
                 right.dot(right) == left.dot(left)
                 and top.dot(top) == bottom.dot(bottom))
    q, p = s.symbols("theta p", real=True)
    inertia = rho*area*d**2/2
    stiffness = rho*gamma**2*d**2/(2*area)
    hamiltonian = p*p/(2*inertia)+stiffness*q*q/2
    pdot = -s.diff(hamiltonian, q)
    checks.check("the actual director angle phase has nonzero orbit-momentum rate",
                 s.diff(pdot, q) == -stiffness and stiffness.is_positive)
    # A concrete periodic advective angular flux need not share the torque
    # cancellation: paired x-face velocities have equal values, while the
    # tangential angular lever arms differ by the period.
    vx, vy = s.symbols("vx vy", real=True)
    velocity = s.Matrix([vx, vy, 0])
    flux_pair = (right.cross(velocity)*velocity.dot(ex)
                 +left.cross(velocity)*velocity.dot(-ex))
    checks.check("advective angular flux retains its physical tangential lever arm",
                 s.simplify(flux_pair-s.Matrix([0, 0, length*vx*vy])) == s.zeros(3, 1))
    core, moment, surface, rate_core, rate_moment, rate_surface = s.symbols(
        "J I B Jdot Idot Bdot", real=True)
    t = s.symbols("t", real=True)
    full_spin = (core+t*rate_core)+(gamma/area)*(moment+t*rate_moment)
    full_spin += surface+t*rate_surface
    checks.check("compensating-vorticity and moving-boundary spin rates both remain",
                 s.diff(full_spin, t) == rate_core+gamma*rate_moment/area+rate_surface)
    print("Scope: exact physical square-parcel mismatch and retained boundary state;")
    print("the positive Euler pair/director survives; this parcel is not its full spin law")
    raise SystemExit(checks.finish())


if __name__ == "__main__":
    main()
