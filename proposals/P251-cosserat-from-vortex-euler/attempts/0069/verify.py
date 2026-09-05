"""Exact first-moment repair and full curl/double-divergence Leray symbol."""

import sympy as s

from substrate_framework.verification import CheckLedger

ledger = CheckLedger("P251/0069")
x, y, z = s.symbols("px py pz", real=True)
p = s.Matrix([x, y, z])
square = p.dot(p)
cross = s.Matrix([[0, -z, y], [z, 0, -x], [-y, x, 0]])
curl = s.I*cross
projector = s.eye(3)-p*p.T/square
a = s.Matrix(s.symbols("a0:3", real=True))
moment = s.Matrix(3, 3, lambda i, j:
                  -sum(s.LeviCivita(i, j, m)*a[m] for m in range(3)))
c = s.Symbol("c", real=True)
stf = s.Matrix([[x, y, z], [y, -x, 1], [z, 1, 0]])
ledger.check("compact curl moment is antisymmetric and carries all three components",
             moment+moment.T == s.zeros(3)
             and moment[0, 1] == -a[2] and moment[1, 2] == -a[0]
             and moment[2, 0] == -a[1])
ledger.check("scalar plus curl moments are orthogonal to symmetric tracefree strain",
             s.expand(sum((c*s.eye(3)+moment)[i, j]*stf[i, j]
                          for i in range(3) for j in range(3))) == 0)
ledger.check("retaining only translation constraints would leave STF moments",
             sum(stf[i, j]**2 for i in range(3) for j in range(3)) != 0)
ledger.check("Leray fixes the actual compact-curl spectral summand",
             s.simplify(projector*curl-curl) == s.zeros(3)
             and s.simplify(curl.H*projector-curl.H) == s.zeros(3))
d = s.zeros(3, 27)
for i in range(3):
    for j in range(3):
        for ell in range(3):
            d[i, 9*i+3*j+ell] = -p[j]*p[ell]
ledger.check("complete mixed kinetic block is polynomial despite the Leray denominator",
             s.simplify(curl.H*projector*d-curl.H*d) == s.zeros(3, 27))
ledger.check("curl kinetic block has the exact transverse polynomial symbol",
             s.expand(curl.H*curl-(square*s.eye(3)-p*p.T)) == s.zeros(3))
# One nonpolynomial entry exposes the actual remaining homogeneous degree.
entry = (d.T*projector*d)[0, 0]
t = s.Symbol("t", positive=True)
ledger.check("double-divergence block is homogeneous degree four, not degree two",
             s.simplify(entry.subs({x: t*x, y: t*y, z: t*z}, simultaneous=True)
                        -t**4*entry) == 0
             and s.simplify(entry-t**2*entry) != 0)
raise SystemExit(ledger.finish())
