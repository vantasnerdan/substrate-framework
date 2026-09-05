"""Actual full Kelvin sector, retained shell and physical energy sign."""

import sympy as s

from substrate_framework import euler_fourier as ef
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0218-full-kelvin-sector")
    aa = s.Rational(1, 100)
    psi = ef.add(ef.trig(2), ef.scale(ef.trig(1), aa))
    u = (psi, ef.trig(2, kind="sin"), ef.scale(ef.trig(1, kind="sin"), -aa))
    streams = [ef.trig(1, 2, "sin"),
               ef.mul(ef.trig(1), ef.trig(2, kind="sin")),
               ef.mul(ef.trig(1, 2), ef.trig(2, kind="sin"))]

    def transport(field):
        return ef.transport(u, (field, {}, {}))[0]

    def inverse_h(field):
        return ef.add({k: value/s.Integer(sum(q*q for q in k))
                       for k, value in field.items() if k != ef.ZERO})

    def scalar_inner(left, right):
        return ef.inner((left, {}, {}), (right, {}, {}))

    xis = [( {}, ef.scale(ef.derivative(f, 2), -1), ef.derivative(f, 1)) for f in streams]
    velocities, hessian, phase = ef.coadjoint_matrices(u, xis, beltrami_eigenvalue=-1)
    pressure = ef.scale(ef.add(*(ef.mul(f, f) for f in u)), -s.Rational(1, 2))
    material_k = ef.material_jacobi_matrix(u, xis, pressure)
    rates = [ef.material_kelvin_operator(u, xi) for xi in xis]
    zetas = [transport(f) for f in streams]
    for index, (zeta, xi, velocity, rate) in enumerate(zip(zetas, xis, velocities, rates, strict=True)):
        phi = inverse_h(zeta)
        expected_velocity = (zeta, ef.scale(ef.derivative(phi, 2), -1), ef.derivative(phi, 1))
        ledger.check(f"full induced axial and planar velocity {index}", all(not ef.add(v, ef.scale(e, -1)) for v, e in zip(velocity, expected_velocity, strict=True)))
        sdot = ef.add(phi, ef.scale(zeta, -1))
        expected_rate = ({}, ef.scale(ef.derivative(sdot, 2), -1), ef.derivative(sdot, 1))
        ledger.check(f"full Kelvin reconstruction has zero axial displacement rate {index}", all(not ef.add(v, ef.scale(e, -1)) for v, e in zip(rate, expected_rate, strict=True)))
        energy_weight = scalar_inner(zeta, ef.add(zeta, ef.scale(phi, -1)))
        ledger.check(f"positive weight is NEGATIVE physical coadjoint Hessian {index}", s.simplify(hessian[index, index]+energy_weight) == 0)
        ledger.check(f"physical full Jacobi Hessian agrees with coadjoint form {index}", s.simplify(material_k[index, index]+ef.inner(rate, rate)-hessian[index, index]) == 0)
        ledger.check(f"nontrivial tested sector has strict negative physical energy {index}", hessian[index, index].is_negative is True)
    expected_phase = s.Matrix(3, 3, lambda i, j: scalar_inner(streams[i], zetas[j]))
    ledger.check("complete Kelvin phase has opposite sign from homogeneous-label phase", s.simplify(phase-expected_phase) == s.zeros(3))
    shell = {k: value for k, value in zetas[1].items() if sum(q*q for q in k) == 1}
    ledger.check("actual high-stream data can force a nonzero first vorticity shell", bool(shell))
    ledger.check("first shell has zero positive energy weight rather than being inverted", not ef.add(shell, ef.scale(inverse_h(shell), -1)))
    field = zetas[1]
    bz = ef.add(field, ef.scale(inverse_h(field), -1))
    lz = ef.scale(transport(bz), -1)
    wrong = ef.scale(transport(field), -1)
    ledger.check("dropping full elliptic reaction changes actual Euler generator", bool(ef.add(lz, ef.scale(wrong, -1))))
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
