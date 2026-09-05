"""Exact phase mass and full kinetic-Gram centering identities."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0083-phase-centering")
    rho, fraction = s.symbols("rho fraction", positive=True)
    velocity = s.Matrix(s.symbols("v0:3", real=True))
    phase_energy = (rho*fraction+rho*(1-fraction))*velocity.dot(velocity)/2
    ledger.check("tube and ambient mass give the same total translation metric",
                 s.simplify(phase_energy-rho*velocity.dot(velocity)/2) == 0)

    # A four-dimensional mass-weighted velocity space; two macro momenta
    # and three reaction columns. The rectangular maps have nonzero crosses.
    b = s.Matrix([[1, 2], [0, 1], [1, 0], [2, -1]])
    a = s.Matrix([[1, 2, -1], [3, 0, 1], [-1, 1, 2], [2, 1, 0]])
    m, g, c = b.T*b, a.T*a, b.T*a
    projector = b*m.inv()*b.T
    residual = g-c.T*m.inv()*c
    ledger.check("macro kinetic projection is self-adjoint and idempotent",
                 projector.T == projector and projector**2 == projector)
    ledger.check("full noncommuting Schur complement is the orthogonal residual Gram",
                 residual == a.T*(s.eye(4)-projector)*a)
    # Exact principal minors establish PSD without rounded small eigenvalues.
    minors = [residual.extract(indices, indices).det()
              for indices in ((0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2))]
    ledger.check("the residual kinetic Gram is positive semidefinite including its null",
                 all(value >= 0 for value in minors) and residual.det() == 0)
    local = s.Matrix([[5, 1, 0], [1, 4, 1], [0, 1, 3]])
    centered = g+local-c.T*m.inv()*c
    ledger.check("positive local helicity survives full ambient-inclusive centering",
                 all(centered[:n, :n].det() > 0 for n in (1, 2, 3)))
    ledger.check("the centered Hessian keeps precisely the local plus residual terms",
                 centered == local+residual)
    p = s.Matrix(s.symbols("p0:2", real=True))
    reaction = s.Matrix(s.symbols("s0:3", real=True))
    before = ((b*p+a*reaction).dot(b*p+a*reaction)
              +(reaction.T*local*reaction)[0])/2
    shift = p+m.inv()*c*reaction
    after = ((shift.T*m*shift)[0]+(reaction.T*centered*reaction)[0])/2
    ledger.check("one complete kinetic square supplies the exact momentum shift",
                 s.expand(before-after) == 0)
    # With A=B the true residual is zero; double subtraction is negative.
    wrong = m-2*m*m.inv()*m
    ledger.check("over-subtracting the same kinetic cross is exposed analytically",
                 wrong == -m and wrong[0, 0] < 0)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
