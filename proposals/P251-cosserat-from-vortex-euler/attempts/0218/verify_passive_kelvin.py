"""Field-changing positive Kelvin lift, with its complete initial forms."""

import sympy as s

from substrate_framework import euler_fourier as ef
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0218-positive-passive-kelvin")
    aa = s.Rational(1, 100)
    psi = ef.add(ef.trig(2), ef.scale(ef.trig(1), aa))
    u = (psi, ef.trig(2, kind="sin"), ef.scale(ef.trig(1, kind="sin"), -aa))
    streams = [ef.trig(1, 2, "sin"),
               ef.mul(ef.trig(1), ef.trig(2, kind="sin")),
               ef.mul(ef.trig(1, 2), ef.trig(2, kind="sin"))]

    def transport(field):
        return ef.transport(u, (field, {}, {}))[0]

    def scalar_inner(left, right):
        return ef.inner((left, {}, {}), (right, {}, {}))

    def add_vectors(left, right, factor=1):
        return tuple(ef.add(a, ef.scale(b, factor)) for a, b in zip(left, right, strict=True))

    xis = [(ef.scale(f, -1), ef.scale(ef.derivative(f, 2), -1), ef.derivative(f, 1)) for f in streams]
    velocities, hessian, phase = ef.coadjoint_matrices(u, xis, beltrami_eigenvalue=-1)
    pressure = ef.scale(ef.add(*(ef.mul(f, f) for f in u)), -s.Rational(1, 2))
    material_k = ef.material_jacobi_matrix(u, xis, pressure)
    rates = [ef.material_kelvin_operator(u, xi) for xi in xis]
    gs = [transport(f) for f in streams]
    for index, (g, xi, velocity, rate) in enumerate(zip(gs, xis, velocities, rates, strict=True)):
        ledger.check(f"actual coadjoint velocity is nonzero axial transport {index}", velocity == (g, {}, {}) and bool(g))
        expected_rate = (g, ef.derivative(g, 2), ef.scale(ef.derivative(g, 1), -1))
        ledger.check(f"full material Kelvin evolution preserves the passive lift {index}", not any(add_vectors(rate, expected_rate, -1)))
        ledger.check(f"full physical coadjoint Hessian is the positive axial norm {index}", s.simplify(hessian[index, index]-scalar_inner(g, g)) == 0)
        ledger.check(f"full Jacobi stiffness plus rate energy has the same sign {index}", s.simplify(material_k[index, index]+ef.inner(rate, rate)-hessian[index, index]) == 0)
        euler_rate = tuple(ef.scale(f, -1) for f in ef.leray(add_vectors(ef.transport(u, velocity), ef.transport(velocity, u))))
        ledger.check(f"actual full Euler pressure gives exactly passive velocity evolution {index}", euler_rate == (ef.scale(transport(g), -1), {}, {}))
    expected_phase = s.Matrix(3, 3, lambda i, j: -scalar_inner(streams[i], gs[j]))
    ledger.check("full Kelvin phase has the derived positive-clock orientation", s.simplify(phase-expected_phase) == s.zeros(3))
    ledger.check("same physical phase does not imply homogeneous-label energy", hessian[0, 0] != 2*scalar_inner(gs[0], gs[0]))
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
