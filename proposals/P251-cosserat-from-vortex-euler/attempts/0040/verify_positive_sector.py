"""Full positive Euler orbit Hessian, KKS sign and physical angle normalization."""

import sympy as s

from fourier_orbit import (
    ZERO, actual_tube, core_cage_generators, cross, curl, derivative,
    divergence, inner, leray, mul, orbit_matrices, scale, trig,
)
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0040-positive-tube-orbit")
    a, b = s.symbols("a b", positive=True)
    background = actual_tube(a, b)
    angle, _ = core_cage_generators(1)
    f = mul(trig(0, kind="sin"), trig(1, kind="sin"))
    partner = (mul(derivative(f, 0), trig(2, kind="sin")),
               mul(derivative(f, 1), trig(2, kind="sin")),
               scale(mul(f, trig(2)), -2))
    generators = (angle, partner)
    tangents, hessian, kks = orbit_matrices(background, generators)
    ledger.check("actual tube background satisfies curl u0=u0", curl(background) == background)
    for index, (generator, tangent) in enumerate(zip(generators, tangents)):
        ledger.check(f"generator {index} is exactly volume preserving", divergence(generator) == {})
        ledger.check(f"isovortical tangent {index} is divergence free", divergence(tangent) == {})
        ledger.check(f"internal tangent {index} has zero harmonic velocity",
                     all(component.get(ZERO, 0) == 0 for component in tangent))
        ledger.check(f"generator {index} is not in the vorticity relabeling kernel",
                     any(curl(tangent)))
        unprojected = cross(generator, curl(background))
        ledger.check(f"Leray leaves the defining vorticity variation {index} intact",
                     curl(tangent) == curl(unprojected))

    expected_h = s.diag(-(17*a*a-60*a*b+17*b*b)/240, 5*(a*a+b*b)/96)
    expected_b = -(a-b)/4
    ledger.check("full two-coordinate Hessian derives with its zero mixed entry",
                 s.simplify(hessian-expected_h) == s.zeros(2))
    ledger.check("KKS pairing is nonzero for a noncircular core",
                 s.simplify(kks-s.Matrix([[0, expected_b], [-expected_b, 0]])) == s.zeros(2))
    # Different route: differentiate the physical energy along the pushed vorticity.
    for i in range(2):
        for j in range(i, 2):
            forward = cross(generators[i], curl(tangents[j]))
            reverse = cross(generators[j], curl(tangents[i]))
            direct = inner(tangents[i], tangents[j])+(inner(background, forward)
                                                      +inner(background, reverse))/2
            ledger.check(f"physical energy second derivative agrees at entry {i},{j}",
                         s.simplify(direct-hessian[i, j]) == 0)
    wrong_h = s.Matrix(2, 2, lambda i, j: inner(tangents[i], tangents[j])
                        +inner(tangents[i], curl(tangents[j])))
    ledger.check("wrong helicity sign changes the physical orbit Hessian", wrong_h != hessian)
    unprojected = cross(angle, curl(background))
    ledger.check("projection mutation exposes a compressible velocity tangent",
                 divergence(unprojected) != {})
    ledger.check("Leray projection is idempotent", leray(tangents[0]) == tangents[0])

    # KKS convention from Euler evolution, without assuming a canonical sign.
    uu = s.Matrix(s.symbols("u0:3"))
    ww = s.Matrix(s.symbols("omega0:3"))
    xx = s.Matrix(s.symbols("xi0:3"))
    d_energy = uu.dot(xx.cross(ww))
    omega_pairing = ww.dot(uu.cross(xx))
    ledger.check("i_X Omega=dE gives Euler generator X=u with positive KKS pairing",
                 s.expand(d_energy-omega_pairing) == 0)
    ledger.check("reversed KKS sign fails the directional energy identity",
                 s.expand(d_energy+omega_pairing) != 0)

    rho, ell = s.symbols("rho ell", positive=True)
    physical_h = rho*hessian.subs(a, 2*b)
    physical_b = rho*ell*kks[0, 1].subs(a, 2*b)
    ledger.check("explicit noncircular geometry has positive full energy Hessian",
                 physical_h == s.diag(7*rho*b*b/48, 25*rho*b*b/96)
                 and physical_h.det().is_positive)
    q, qdot, shape = s.symbols("q q_dot shape", real=True)
    action = physical_b*shape*qdot-(physical_h[0, 0]*q*q+physical_h[1, 1]*shape*shape)/2
    eliminated = s.solve(s.diff(action, shape), shape)[0]
    reduced = s.factor(action.subs(shape, eliminated))
    inertia = s.factor(physical_b**2/physical_h[1, 1])
    ledger.check("conjugate shape elimination gives positive angle inertia",
                 inertia == 6*rho*ell**2/25)
    ledger.check("same physical action gives kinetic minus restoring potential",
                 s.simplify(reduced-(inertia*qdot*qdot-physical_h[0, 0]*q*q)/2) == 0)
    frequency2 = s.factor(physical_h[0, 0]/inertia)
    ledger.check("one-action frequency is 175*b^2/(288*ell^2)",
                 frequency2 == 175*b*b/(288*ell**2))

    # Local core angle and conjugate return are geometrical, not merely variable names.
    x, y, z = s.symbols("x y z", real=True)
    angle_field = s.Matrix([-s.sin(y), s.sin(x), 0])*s.cos(z)
    partner_field = s.Matrix([s.cos(x)*s.sin(y)*s.sin(z),
                              s.sin(x)*s.cos(y)*s.sin(z),
                              -2*s.sin(x)*s.sin(y)*s.cos(z)])
    angle_jet = angle_field[:2, :].jacobian((x, y)).subs({x: 0, y: 0})
    ledger.check("opposite sections undergo opposite physical core rotations",
                 angle_jet.subs(z, 0) == s.Matrix([[0, -1], [1, 0]])
                 and angle_jet.subs(z, s.pi) == s.Matrix([[0, 1], [-1, 0]]))
    ledger.check("axial shape return closes incompressibility directly in real space",
                 s.simplify(sum(s.diff(partner_field[j], coord)
                                for j, coord in enumerate((x, y, z)))) == 0)
    ledger.check("physical relative angle chi=2q rescales both coefficients together",
                 inertia/4 == 3*rho*ell**2/50
                 and physical_h[0, 0]/4 == 7*rho*b*b/192)

    translations = [tuple({ZERO: s.Integer(1)} if j == i else {} for j in range(3))
                    for i in range(3)]
    _, extended_h, extended_kks = orbit_matrices(background, (*translations, *generators))
    ledger.check("uniform translations have zero orbit spring and spring cross terms",
                 extended_h[:3, :] == s.zeros(3, 5))
    ledger.check("uniform translation KKS cross terms vanish on this internal leaf",
                 extended_kks[:3, :] == s.zeros(3, 5))
    print("H per volume at a=2b:", physical_h)
    print("B per volume:", physical_b)
    print("I_q per volume:", inertia)
    print("I_relative per volume:", inertia/4)
    print("K_relative per volume:", physical_h[0, 0]/4)
    print("omega_squared:", frequency2)
    print("Scope: positive two-generator isovortical action of an actual smooth Beltrami tube.")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
