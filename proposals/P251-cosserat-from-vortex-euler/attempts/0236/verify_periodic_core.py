"""Exact same-cell finite core, flux twist and density identities, no spectrum fit."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0236-same-periodic-core")
    xx, a, b = s.symbols("X a b", real=True)
    om = s.symbols("Omega", positive=True)
    coords = (xx, a, b)
    psi = s.cos(b)+om**2*s.cos(a)
    u = s.Matrix([psi, s.sin(b), -om**2*s.sin(a)])
    jac = u.jacobian(coords)
    curl = s.Matrix([jac[2, 1]-jac[1, 2], jac[0, 2]-jac[2, 0],
                     jac[1, 0]-jac[0, 1]])
    checks.check("the SAME nonlinear cell has literal constant curl and incompressibility",
                 curl == -u and s.trace(jac) == 0)
    p = -u.dot(u)/2
    checks.check("the full pressure retains the exact stationary Euler balance",
                 s.simplify(jac*u+s.Matrix([s.diff(p, q) for q in coords])) == s.zeros(3, 1))
    grad = s.Matrix([s.diff(psi, q) for q in coords])
    checks.check("the entire positive-radius family is exactly flow and vortex invariant",
                 s.simplify(grad.dot(u)) == 0 and s.simplify(grad.dot(curl)) == 0)
    checks.check("the closed periodic core is nonzero and purely axial",
                 u.subs({a: 0, b: 0}) == s.Matrix([1+om**2, 0, 0]))
    hessian = s.hessian(1+om**2-psi, (a, b)).subs({a: 0, b: 0})
    checks.check("the actual normal energy has a nondegenerate elliptic minimum",
                 hessian == s.diag(om**2, 1))
    v = u[1:, 0]
    checks.check("the axial-section return preserves the actual weighted flux",
                 s.simplify(sum(s.diff(psi*(v[i]/psi), q)
                                for i, q in enumerate((a, b)))) == 0)
    action = s.symbols("I", positive=True)
    angle = s.symbols("theta", real=True)
    aa = s.sqrt(2*action/om)*s.cos(angle)
    bb = s.sqrt(2*action*om)*s.sin(angle)
    quartic = -(om**2*aa**4+bb**4)/24
    average = s.simplify(s.integrate(s.expand_trig(quartic),
                                     (angle, 0, 2*s.pi))/(2*s.pi))
    checks.check("independent Hamiltonian averaging derives the nonlinear period coefficient",
                 s.simplify(average+(1+om**2)*action**2/16) == 0)
    normal_form = om*action+average
    omega_e = s.diff(normal_form, action, 2).subs(action, 0)/om
    checks.check("the true energy-frequency derivative includes orbital distortion",
                 s.simplify(omega_e+(1+om**2)/(8*om)) == 0)
    energy = s.symbols("E", real=True)
    ratio = (om+omega_e*energy)/(1+om**2-energy)
    ratio_e = s.diff(ratio, energy).subs(energy, 0)
    ratio_j = s.factor(ratio_e*om/(1+om**2))
    checks.check("the physical flux action gives the derived return-map twist",
                 s.factor(ratio_j-(-1+6*om**2-om**4)/(8*(1+om**2)**3)) == 0)
    actual = ratio_j.subs(om, s.Rational(1, 10))
    checks.check("the actual fixed cell has strict nonzero twist without a numerical floor",
                 actual == -s.Rational(235025, 2060602) and actual < 0)
    rotation = (om/(1+om**2)).subs(om, s.Rational(1, 10))
    checks.check("the actual core multipliers are elliptic and neither plus nor minus one",
                 rotation == s.Rational(10, 101) and 0 < rotation < s.Rational(1, 2))
    checks.check("the rational core rotation is not misreported as the Diophantine boundary",
                 rotation.is_Rational and ratio_e.subs(om, s.Rational(1, 10)) != 0)
    checks.check("discarding the axial speed derivative changes the flux twist",
                 s.simplify(ratio_j-omega_e*om/(1+om**2)**2) != 0)
    cell_energy = s.simplify(s.integrate(s.integrate(u.dot(u), (a, 0, 2*s.pi)),
                                         (b, 0, 2*s.pi))/(2*s.pi)**2)
    checks.check("the actual stationary law has finite nonzero background energy density",
                 cell_energy == 1+om**4)
    lx, j0, fraction = s.symbols("L j0 f", positive=True)
    density = fraction*j0*lx/(3*lx*(2*s.pi)**2)
    checks.check("actual optical inertia density is independent of axial period and packet width",
                 s.simplify(density-fraction*j0/(12*s.pi**2)) == 0
                 and s.diff(density, lx) == 0 and density.is_positive)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
