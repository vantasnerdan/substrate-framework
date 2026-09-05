"""Retain ambient transport zero modes while deriving their leading mean decoupling."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0150-retained-ambient-kernel")
    eps, z, omega = s.symbols("epsilon z omega", positive=True)
    pi = s.diag(1, 0)
    transport = s.diag(0, s.I*omega)
    resolvent = (eps*z*s.eye(2)+transport).inv()
    scaled_limit = (eps*resolvent).applyfunc(lambda value: s.limit(value, eps, 0))
    checks.check("the ambient orbit-average projection is retained in the acoustic limit",
                 scaled_limit == pi/z)
    a, b, c, d = s.symbols("a b c d")
    coupling = s.Matrix([[a, b], [c, d]])
    feedback = eps**2*resolvent*coupling*resolvent
    feedback_limit = feedback.applyfunc(lambda value: s.limit(value, eps, 0))
    checks.check("both compact feedback resolvents keep the actual averaged block",
                 feedback_limit == pi*coupling*pi/z**2)
    primitive = s.Matrix([0, 1/(s.I*omega)])
    velocity = transport*primitive
    checks.check("bounded streamline coordinates force zero orbit-average velocity",
                 pi*velocity == s.zeros(2, 1) and velocity == s.Matrix([0, 1]))
    dressed = (s.eye(2)-feedback_limit).inv()*velocity
    checks.check("the full averaged response is invisible to the leading physical mean",
                 s.simplify((velocity.T*dressed)[0]-1) == 0
                 and feedback_limit != s.zeros(2))
    row = s.Matrix([[2, 3]])*transport
    checks.check("the actual translation primitive annihilates ambient zero modes",
                 row*pi == s.zeros(1, 2))
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
