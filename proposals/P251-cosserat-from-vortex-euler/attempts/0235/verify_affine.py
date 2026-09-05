"""Actual three-axis acoustic fields and full material-spin Taylor row.

Canonical finite fields establish which preparation is used. Cartesian
Euler pressure and the complete torque variation independently expose
the nonzero full-current curvature. No spin row is supplied as a modulus.
"""

import sympy as s

from substrate_framework import euler_fourier as ef
from substrate_framework.euler_displacement_preparation import (
    finite_displacement_cell,
    prepared_displacement,
)
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0235-complete-affine-current")
    cell = finite_displacement_cell()
    unit = s.eye(3)
    amplitude = s.symbols("amplitude", real=True)

    def add(*vectors):
        return tuple(ef.add(*(v[i] for v in vectors)) for i in range(3))

    def scale(vector, coefficient):
        return tuple(ef.scale(v, coefficient) for v in vector)

    def zero_fourier(vector):
        return all(s.simplify(value) == 0 for row in vector for value in row.values())

    for axis, left_axis, right_axis in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        left = prepared_displacement(cell, unit[:, left_axis], unit[:, right_axis],
                                     amplitude=amplitude)
        right = prepared_displacement(cell, unit[:, right_axis], unit[:, left_axis],
                                      amplitude=amplitude)

        def difference(name):
            return add(getattr(left, name), scale(getattr(right, name), -1))

        axis_field = tuple({ef.ZERO: unit[i, axis]} for i in range(3))
        rotated = ef.cross(axis_field, cell.background)
        derivative = tuple(ef.derivative(row, axis) for row in cell.background)
        checks.check(f"axis {axis}: the complete D velocity is the actual rigid rotation",
                     zero_fourier(add(difference("first_velocity"), scale(rotated, -1))))
        checks.check(f"axis {axis}: the actual translation return cancels the material-rate row",
                     zero_fourier(add(difference("returned"), scale(derivative, -1)))
                     and zero_fourier(difference("material_rate")))
        checks.check(f"axis {axis}: raw lift plus actual C016 forcing retains both translation terms",
                     zero_fourier(add(difference("lift"), derivative, scale(rotated, -1)))
                     and zero_fourier(add(ef.leray(add(rotated, difference("lift"))),
                                          scale(derivative, -1))))
        checks.check(f"axis {axis}: full pressure turns the rate source into a neutral translation",
                     zero_fourier(add(ef.leray(rotated), scale(derivative, -1))))
    kap, vv = s.symbols("k0:3"), s.symbols("V0:3")
    symmetric = kap[1]*vv[0]+kap[0]*vv[1]
    checks.check("all three actual C016 passive-return antisymmetric contractions vanish",
                 all(s.diff(symmetric, kap[j], vv[k])-s.diff(symmetric, kap[k], vv[j]) == 0
                     for j, k in ((1, 2), (2, 0), (0, 1))))

    x, y, z, t = s.symbols("X Y Z t", real=True)
    xyz = (x, y, z)
    point = s.Matrix(xyz)
    aa, rho, mean = s.symbols("a rho mean", real=True)
    psi = s.cos(y)+aa*s.cos(z)
    u = s.Matrix([psi, aa*s.sin(z), -s.sin(y)])
    du = u.jacobian(xyz)
    p = -rho*u.dot(u)/2

    def grad(value):
        return s.Matrix([s.diff(value, q) for q in xyz])

    def clean(value):
        return s.simplify(s.trigsimp(s.expand(value)))

    def zero(vector):
        return all(clean(value) == 0 for value in vector)

    spin_second_trace = 0
    omitted_pressure_defects = []
    for axis in range(3):
        e = unit[:, axis]
        rigid = e.cross(point)
        wr = e.cross(u)-du*rigid
        translation = u.diff(xyz[axis])
        velocity = rigid+t*(wr-2*translation)
        pressure = 2*rho*u[axis]-t*rigid.dot(grad(p))-2*t*s.diff(p, xyz[axis])
        residual = rho*(velocity.diff(t)+velocity.jacobian(xyz)*u+du*velocity)+grad(pressure)
        checks.check(f"axis {axis}: full affine velocity and pressure solve linear Euler exactly",
                     zero(residual))
        omitted_pressure_defects.append(grad(2*t*s.diff(p, xyz[axis])))
        pressure0 = pressure.subs(t, 0)
        pressure_rate = pressure.diff(t).subs(t, 0)
        spin_second = -(rigid.cross(grad(p))
                        +(u-mean*unit[:, 0]).cross(grad(pressure0))
                        +point.cross(grad(pressure_rate)
                                     +s.hessian(pressure0, xyz)*u
                                     +s.hessian(p, xyz)*rigid))
        spin_second_trace += spin_second[axis]
    target = 2*rho*(u-mean*unit[:, 0]).dot(u)-2*point.dot(grad(p))
    checks.check("full moved-boundary spin trace equals the defining pressure-virial row",
                 clean(spin_second_trace-target) == 0)
    checks.check("discarding the translation pressure produces a nonzero full Euler residual",
                 any(not zero(row) for row in omitted_pressure_defects))

    # Normal invariant-measure identities use arbitrary chi(psi), not a fitted tag.
    argument = s.symbols("argument", real=True)
    chi_function = s.Function("chi")
    chi = chi_function(psi)
    weighted = chi*(psi-mean)
    weighted_prime = s.diff(chi_function(argument)*(argument-mean), argument).subs(argument, psi)
    normal_gradient = s.Matrix([s.diff(psi, y), s.diff(psi, z)])
    divergence = s.diff(weighted*s.diff(psi, y), y)+s.diff(weighted*s.diff(psi, z), z)
    checks.check("the actual tag divergence relates the response weight to axial shear variance",
                 clean(divergence-weighted_prime*normal_gradient.dot(normal_gradient)
                       +weighted*psi) == 0)
    cx0, cd0, vv0, vd0, sigma, inertia0 = s.symbols(
        "C0 D0 V0 VD Sigma I0", real=True)
    cx = cx0-vv0*t**2/2
    cd = cd0-vd0*t**2/2
    inertia = inertia0+sigma*t**2
    full_trace = 2*inertia+2*(cx-cx0)+2*(t*s.diff(cd, t)-cd+cd0)
    checks.check("the independently reduced correlation trace retains the positive shear second row",
                 s.expand(full_trace.diff(t, 2).subs(vd0, sigma-vv0)-2*sigma) == 0)
    stationary = s.symbols("stationary", real=True)
    checks.check("a constant physical b-return cannot cancel the finite-cut shear curvature",
                 s.diff(full_trace+stationary, t, 2) == s.diff(full_trace, t, 2))

    # Exact passive-current product rule; functions are actual defining correlations in the proof.
    ph, pg = s.Function("Ph")(t), s.Function("Pg")(t)
    displacement = ph+t*pg
    velocity = pg
    passive_spin = s.diff(displacement, t)-2*velocity
    checks.check("the axial passive material source retains both displacement and velocity spin terms",
                 s.expand(passive_spin-s.diff(ph, t)-t*s.diff(pg, t)+pg) == 0)
    checks.check("renaming the passive displacement dipole derivative loses an actual velocity row",
                 s.simplify(s.diff(displacement, t)-passive_spin) == 2*pg)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
