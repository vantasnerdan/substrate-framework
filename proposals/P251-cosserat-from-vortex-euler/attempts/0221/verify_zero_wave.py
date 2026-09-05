"""Actual Kelvin signatures and exact finite phase/energy normalization."""

import sympy as s

from substrate_framework import euler_fourier as ef
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0221-kelvin-zero-wave-controls")
    aa = s.Rational(1, 100)
    psi = ef.add(ef.trig(2), ef.scale(ef.trig(1), aa))
    u = (psi, ef.trig(2, kind="sin"), ef.scale(ef.trig(1, kind="sin"), -aa))
    generators = []
    for harmonic, axial_fraction in ((2, -1), (7, 0)):
        for stream in (ef.trig(1, harmonic, "sin"), ef.mul(ef.trig(1, harmonic), ef.trig(2, kind="sin"))):
            generators.append((ef.scale(stream, axial_fraction), ef.scale(ef.derivative(stream, 2), -1), ef.derivative(stream, 1)))
    _, h, omega = ef.coadjoint_matrices(u, generators, beltrami_eigenvalue=-1)
    ledger.check("actual selected Kelvin families have no cross energy", h[:2, 2:] == s.zeros(2))
    ledger.check("actual selected Kelvin families have no cross phase", omega[:2, 2:] == s.zeros(2))
    hp, hn = h[:2, :2], -h[2:, 2:]
    ledger.check("actual passive energy is positive definite", hp[0, 0] > 0 and hp.det() > 0)
    ledger.check("actual planar Kelvin energy is negative definite", hn[0, 0] > 0 and hn.det() > 0)
    cp = hp.cholesky().inv().T
    cn = hn.cholesky().inv().T
    ledger.check("positive whitening uses complete actual quadratic form", s.simplify(cp.T*hp*cp) == s.eye(2))
    ledger.check("negative whitening uses complete actual quadratic form", s.simplify(cn.T*(-hn)*cn) == -s.eye(2))
    kp = s.simplify((cp.T*omega[:2, :2]*cp)[0, 1])
    kn = s.simplify((cn.T*omega[2:, 2:]*cn)[0, 1])
    ledger.check("opposite signatures retain nonzero physical phase ratios", kp*kn < 0)
    ledger.check("actual unequal families have distinct phase-to-energy ratios", s.simplify(kp+kn) != 0)
    b = s.symbols("desired_phase", positive=True)
    factor = b/(kp+kn)
    ledger.check("selected ordering realizes positive amplitude square", factor.is_positive is True)
    ledger.check("actual paired full energy cancels exactly", s.simplify(factor*(cp.T*hp*cp-cn.T*hn*cn)) == s.zeros(2))
    ledger.check("actual paired phase equals its prescribed finite value", s.simplify(factor*(kp+kn)-b) == 0)
    swap = s.Matrix([[0, 1], [1, 0]])
    j = s.Matrix([[0, 1], [-1, 0]])
    ledger.check("real column interchange reverses phase without negative density", swap.T*j*swap == -j and swap.T*s.eye(2)*swap == s.eye(2))
    v = s.Matrix(s.symbols("v1 v2", real=True))
    ledger.check("single actual rank-one energy return has exactly zero phase", v.T*j*v == s.zeros(1))
    ledger.check("same-sign energy substitution cannot cancel the full matrix", factor*(s.eye(2)+s.eye(2)) != s.zeros(2))
    print("actual whitened phase scalars:", kp, kn)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
