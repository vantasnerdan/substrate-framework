"""Exact joint-jet triangularity with all intrinsic control rows retained.

Actual Kelvin blocks are supplied by reviewed0221. This verifier checks
the new polynomial input-map implication, not their Euler existence.
Phase is skew Hermitian, energy Hermitian, with real Bloch conjugation.
"""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0228-joint-kelvin-jet-controls")
    k = s.symbols("k", real=True)
    p, q, r, z, a, b, c, d = s.symbols("p q r z a b c d", real=True)
    unit = s.Matrix([[0, 1], [-1, 0]])
    phase = unit+s.I*k*s.Matrix([[p, q], [q, r]])+k**2*z*unit
    energy = s.I*k*a*unit+k**2*s.Matrix([[b, c], [c, d]])
    checks.check("complete phase-block jets have the required two reality symmetries",
                 phase.conjugate().T == -phase and phase.subs(k, -k) == phase.conjugate())
    checks.check("complete energy-block jets retain allowed odd Hermitian cross terms",
                 energy.conjugate().T == energy and energy.subs(k, -k) == energy.conjugate()
                 and energy.diff(k).subs(k, 0) != s.zeros(2))
    v = s.Matrix(s.symbols("v1 v2", real=True))
    w = s.Matrix(s.symbols("w1 w2", real=True))
    plus = v.T+s.I*k*w.T
    minus = v.T
    hp, hm, sp, sm = s.symbols("hp hm sp sm", real=True)
    e_first = plus.conjugate().T*(1+k**2*hp)*plus+minus.T*(-1+k**2*hm)*minus
    o_first = plus.conjugate().T*(s.I*k*sp)*plus+minus.T*(s.I*k*sm)*minus
    checks.check("actual opposite-energy rows cancel their full order-zero matrix",
                 e_first.subs(k, 0) == s.zeros(2))
    checks.check("their first energy is the requested imaginary skew input wedge",
                 s.simplify(e_first.diff(k).subs(k, 0)-s.I*(v*w.T-w*v.T)) == s.zeros(2))
    checks.check("their actual induced first phase is retained as a symmetric rank-one row",
                 s.simplify(o_first.diff(k).subs(k, 0)-s.I*(sp+sm)*v*v.T) == s.zeros(2))
    target = s.symbols("target", real=True)
    amplitude = (target-p)/2
    return_map = s.Matrix.vstack(v.T, s.I*k*amplitude*v.T)
    repaired_phase = s.simplify(return_map.conjugate().T*phase*return_map)
    repaired_energy = s.simplify(return_map.conjugate().T*energy*return_map)
    checks.check("first-phase correction leaves both constant forms zero",
                 repaired_phase.subs(k, 0) == s.zeros(2)
                 and repaired_energy.subs(k, 0) == s.zeros(2))
    checks.check("the intrinsic diagonal phase is included in the exact correction",
                 s.simplify(repaired_phase.diff(k).subs(k, 0)-s.I*target*v*v.T) == s.zeros(2))
    checks.check("first-phase correction does not undo the first-energy matching",
                 repaired_energy.diff(k).subs(k, 0) == s.zeros(2))
    checks.check("the intrinsic odd energy still contributes to the later second row",
                 s.simplify(repaired_energy.diff(k, 2).subs(k, 0)/2
                            -(b-a*(target-p))*v*v.T) == s.zeros(2))
    ka, kb = s.symbols("Ka Kb", real=True)
    second_map = s.I*k*s.Matrix.vstack(ka*v.T, kb*w.T)
    second_phase = second_map.conjugate().T*phase*second_map
    second_energy = second_map.conjugate().T*energy*second_map
    checks.check("linear real-Bloch amplitude rows supply the exact mixed quadratic phase",
                 s.simplify(second_phase.diff(k, 2).subs(k, 0)/2
                            -ka*kb*(v*w.T-w*v.T)) == s.zeros(2))
    checks.check("that quadratic phase return has no energy through degree two",
                 all(second_energy.diff(k, order).subs(k, 0) == s.zeros(2)
                     for order in range(3)))
    scalar_map = s.I*k*ka*v.T
    scalar_phase = scalar_map.conjugate().T*(s.I*k*sp)*scalar_map
    scalar_energy = scalar_map.conjugate().T*(1+k**2*hp)*scalar_map
    checks.check("linear scalar-energy return supplies its exact quadratic rank-one energy",
                 s.simplify(scalar_energy.diff(k, 2).subs(k, 0)/2-ka**2*v*v.T) == s.zeros(2))
    checks.check("the scalar-energy return leaves phase zero through degree two",
                 all(scalar_phase.diff(k, order).subs(k, 0) == s.zeros(2)
                     for order in range(3)))
    checks.check("signed actual energy squares span mixed spatial monomials",
                 s.expand(((ka+kb)**2-(ka-kb)**2)/4-ka*kb) == 0)
    checks.check("omitting the intrinsic phase row produces a nonzero normalization defect",
                 s.simplify(p+2*(target/2)-target) == p and p != 0)
    checks.check("dropping the opposite-energy companion leaves a lower-order debt",
                 (plus.conjugate().T*plus).subs(k, 0) == v*v.T
                 and v*v.T != s.zeros(2))
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
