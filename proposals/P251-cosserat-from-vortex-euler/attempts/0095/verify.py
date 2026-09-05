"""Exact complement dynamics, Schur action, and physical observation repair."""

import sympy as s

from substrate_framework.euler_orbit import hermitian_schur_jet
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0095-Kelvin-complement")
    omega = s.Symbol("omega", real=True)
    j2 = s.Matrix([[0, 1], [-1, 0]])
    symplectic = s.diag(j2, j2)
    hessian = s.Matrix([[3, 0, 0, 0], [0, 2, 0, 1],
                        [0, 0, 7, 0], [0, 1, 0, 5]])
    ledger.check("exposing full Hamiltonian has a positive definite Hessian",
                 all(hessian[:n, :n].det() > 0 for n in range(1, 5)))
    evolution = -symplectic.inv()*hessian
    embedding = s.eye(4)[:, :2]
    omega_e = embedding.T*symplectic*embedding
    projection = omega_e.inv()*embedding.T*symplectic
    a_e = projection*evolution*embedding
    residual = embedding*a_e-evolution*embedding
    ledger.check("finite symplectic projection preserves retained coordinates",
                 projection*embedding == s.eye(2))
    ledger.check("restricted reconstruction residual is symplectically orthogonal but nonzero",
                 embedding.T*symplectic*residual == s.zeros(2)
                 and residual == s.Matrix([[0, 0], [0, 0], [0, -1], [0, 0]]))
    ledger.check("residual norm is explicit rather than assumed negligible",
                 residual.T*residual == s.diag(0, 1))

    pencil = hessian-s.I*omega*symplectic
    retained, coupling, reaction = pencil[:2, :2], pencil[2:, :2], pencil[2:, 2:]
    reduced = s.simplify(retained-coupling.conjugate().T*reaction.inv()*coupling)
    p_eff = s.factor(reduced[1, 1])
    j_eff = s.factor(1/p_eff)
    ledger.check("complete reaction elimination retains the actual frequency dependence",
                 s.factor(p_eff-(2-7/(35-omega**2))) == 0)
    ledger.check("physical retained momentum has the repaired inertia, not the frozen value",
                 j_eff == (omega**2-35)/(2*omega**2-63)
                 and j_eff.subs(omega, 0) == s.Rational(5, 9)
                 and j_eff.subs(omega, 0) != s.Rational(1, 2))
    characteristic = s.factor(pencil.det())
    ledger.check("full characteristic polynomial is recovered by the same Schur action",
                 characteristic == omega**4-41*omega**2+189
                 and s.factor(reduced.det()*reaction.det()-characteristic) == 0)
    ledger.check("frozen single-pair pole does not solve the full retained dynamics",
                 characteristic.subs(omega**2, 6) != 0)

    jet = hermitian_schur_jet(
        (hessian[2:, 2:], -s.I*j2, s.zeros(2)),
        (hessian[2:, :2], s.zeros(2), s.zeros(2)),
        (hessian[:2, :2], -s.I*j2, s.zeros(2)))
    ledger.check("importable noncommuting Schur jets agree with the complete exact pencil",
                 all(s.simplify(jet.reduced[n]-reduced.diff(omega, n).subs(omega, 0)
                                /s.factorial(n)) == s.zeros(2) for n in range(3)))
    ledger.check("first dynamic inertia correction is derived rather than discarded",
                 s.series(j_eff, omega, 0, 4).removeO()
                 == s.Rational(5, 9)+omega**2/s.Integer(567))

    obs_e, obs_r = s.Matrix([[1, 0]]), s.Matrix([[2, 0]])
    obs_eff = s.simplify(obs_e-obs_r*reaction.inv()*coupling)
    full_state = s.Matrix.vstack(s.eye(2), -reaction.inv()*coupling)
    ledger.check("physical field map is transformed together with the eliminated action",
                 s.simplify(s.Matrix.hstack(obs_e, obs_r)*full_state-obs_eff) == s.zeros(1, 2)
                 and s.simplify(obs_eff-s.Matrix([[1, -2*s.I*omega/(35-omega**2)]]))
                 == s.zeros(1, 2))
    ledger.check("keeping the old physical field row after reaction elimination is exposed",
                 obs_eff != obs_e)
    q0, p0 = s.symbols("q0 p0", real=True)
    initial = s.Matrix([q0, p0, 0, 0])
    exact_third = (evolution**3*initial)[0]
    frozen_third = (a_e**3*s.Matrix([q0, p0]))[0]
    ledger.check("zero initial complement still produces a computed later discrepancy",
                 s.expand(exact_third-frozen_third) == -7*p0)
    ledger.check("eliminated state solves both exact reaction rows",
                 s.simplify(pencil[2:, :]*full_state) == s.zeros(2))
    print("EXACT full characteristic:", characteristic)
    print("EXACT reduced inertia:", j_eff)
    print("EXACT physical observation map:", obs_eff)
    print("EXACT third-derivative discrepancy:", exact_third-frozen_third)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
