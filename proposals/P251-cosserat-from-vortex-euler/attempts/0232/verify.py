"""Exact pressure virtual work and null-Lagrangian current identities.

The smooth polynomial fields expose differential signs independently of
the distributional proof. The segment tests use general cubic virtual
fields restricted to the segment, not delta-function numerics. No Euler
mode, continuum modulus or free-boundary constitutive law is assumed.
"""

import sympy as s

from substrate_framework.micropolar import (
    MicropolarCoefficients,
    isotropic_micropolar_energy,
)
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0232-material-torque-boundary")
    x, y, z = xyz = s.symbols("x y z", real=True)
    point = s.Matrix(xyz)
    eye = s.eye(3)

    def grad(value):
        return s.Matrix([s.diff(value, item) for item in xyz])

    def div_matrix(value):
        return s.Matrix([sum(s.diff(value[i, j], xyz[j]) for j in range(3))
                         for i in range(3)]).applyfunc(s.expand)

    def curl(value):
        return s.Matrix([s.diff(value[2], y)-s.diff(value[1], z),
                         s.diff(value[0], z)-s.diff(value[2], x),
                         s.diff(value[1], x)-s.diff(value[0], y)])

    def zero(value):
        return all(s.expand(item) == 0 for item in value)

    xi = curl(s.Matrix([x*y*z, x**2*z+y*z**2, x*y**2+z**3]))
    dxi = xi.jacobian(xyz)
    p = x**3+y**2*z+2*x*z**2+x*y
    p1 = x*z+y**3+z**2
    center = s.Matrix(s.symbols("X1 X2 X3", real=True))
    dcenter = s.Matrix(s.symbols("dX1 dX2 dX3", real=True))
    radius = point-center
    traction = (p1+xi.dot(grad(p)))*eye-p*dxi.T
    force_volume = grad(p1)+s.hessian(p, xyz)*xi
    checks.check("the actual displacement fixture is divergence free",
                 s.trace(dxi) == 0 and xi != s.zeros(3, 1))
    checks.check("moved pressure and area vector give the exact volume force variation",
                 zero(div_matrix(traction)-force_volume))
    torque = s.Matrix.hstack(*[(xi-dcenter).cross(p*eye[:, j])
                              +radius.cross(traction[:, j]) for j in range(3)])
    torque_volume = (xi-dcenter).cross(grad(p))+radius.cross(force_volume)
    checks.check("lever arm and moved normal yield the exact complete torque variation",
                 zero(div_matrix(torque)-torque_volume))
    wrong_traction = (p1+xi.dot(grad(p)))*eye
    checks.check("dropping the moving-normal term exposes a nonzero pressure-force defect",
                 not zero(div_matrix(wrong_traction)-force_volume))
    wrong_torque = s.Matrix.hstack(*[radius.cross(traction[:, j]) for j in range(3)])
    checks.check("dropping the lever-arm variation exposes a nonzero torque defect",
                 not zero(div_matrix(wrong_torque)-torque_volume))
    u = curl(s.Matrix([y*z**2, x*y*z, x**2*y]))
    w = curl(s.Matrix([x*z, x*y**2, x*y*z]))
    convective = w.jacobian(xyz)*u+u.jacobian(xyz)*w
    div_convective = sum(s.diff(convective[i], xyz[i]) for i in range(3))
    checks.check("full linear Euler pressure retains both convective terms",
                 s.expand(div_convective-2*s.trace(u.jacobian(xyz)*w.jacobian(xyz))) == 0
                 and s.expand(div_convective) != 0)

    t = s.symbols("t", real=True)
    f = s.Matrix(s.symbols("F1 F2 F3", real=True))
    r = s.Matrix(s.symbols("r1 r2 r3", real=True))
    rows = s.symbols("v0:12", real=True)
    test = s.Matrix([sum(rows[4*i+j]*t**j for j in range(4)) for i in range(3)])
    bond_force = -s.integrate(f.dot(test.diff(t)), (t, 0, 1))
    checks.check("star-bond force has the true centroid and boundary endpoint signs",
                 s.expand(bond_force-f.dot(test.subs(t, 0)-test.subs(t, 1))) == 0)
    moment = r.cross(f)
    div_mu = -s.integrate((1-t)*moment.dot(test.diff(t)), (t, 0, 1))
    axial_sigma = -s.integrate(moment.dot(test), (t, 0, 1))
    checks.check("complete star couple flux minus axial stress is the centroid torque",
                 s.expand(div_mu-axial_sigma-moment.dot(test.subs(t, 0))) == 0)
    checks.check("reversing the axial-stress convention without its partner breaks balance",
                 s.expand(div_mu+axial_sigma-moment.dot(test.subs(t, 0))) != 0)

    entries = s.symbols("g0:9", real=True)
    g = s.Matrix(3, 3, entries)
    cs, ca, ct, e = s.symbols("cs ca ct e", real=True)
    coefficients = MicropolarCoefficients(0, 0, 0, ct, cs, ca)
    energy = isotropic_micropolar_energy(s.zeros(3), s.zeros(3, 1), g, coefficients)
    delta_energy = s.expand(energy.subs({cs: cs+e, ca: ca-e, ct: ct-e})-energy)
    delta_m = s.Matrix(3, 3, lambda i, j: s.diff(delta_energy, g[i, j]))
    checks.check("canonical energy variation derives the complete null-Lagrangian",
                 s.expand(delta_energy-e*(s.trace(g*g)-s.trace(g)**2)) == 0)
    checks.check("canonical couple-stress derivative has the required transpose and factor two",
                 zero(delta_m-2*e*(g.T-s.trace(g)*eye)))
    phi = s.Matrix([x*y+y*z**2, x**2*z+y**3, x*y*z+x**3])
    dphi = phi.jacobian(xyz)
    local_m = delta_m.subs(dict(zip(entries, list(dphi))))
    checks.check("compatible rotation gradients make the full stress improvement divergence free",
                 zero(div_matrix(local_m)))
    vector_flux = dphi*phi-phi*s.trace(dphi)
    div_flux = sum(s.diff(vector_flux[i], xyz[i]) for i in range(3))
    local_energy = delta_energy.subs(dict(zip(entries, list(dphi))))
    checks.check("the energy difference is exactly its retained boundary divergence",
                 s.expand(e*div_flux-local_energy) == 0)
    superpotential_div = s.Matrix(3, 3, lambda i, j: sum(
        s.diff(2*e*(eye[i, k]*phi[j]-eye[i, j]*phi[k]), xyz[k]) for k in range(3)))
    checks.check("the explicit antisymmetric superpotential reconstructs the current improvement",
                 zero(superpotential_div-local_m))
    eta = s.Matrix([x**2*y+z, y*z+x, z**2*x+y])
    virtual_power = sum(local_m[i, j]*s.diff(eta[i], xyz[j])
                        for i in range(3) for j in range(3))
    boundary_power = local_m.T*eta
    checks.check("boundary virtual work equals the full null-action variation pointwise",
                 s.expand(virtual_power-sum(s.diff(boundary_power[j], xyz[j])
                                            for j in range(3))) == 0
                 and s.expand(virtual_power) != 0)
    wrong_m = 2*e*(dphi-s.trace(dphi)*eye)
    checks.check("omitting the transpose creates a spurious bulk couple force",
                 not zero(div_matrix(wrong_m)))
    checks.check("identical bulk balances do not imply identical point boundary tractions",
                 not zero(local_m*eye[:, 0]))
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
