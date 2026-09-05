"""Finite exact isotropization of the source Beltrami spectral covariance."""

import sympy as s

from substrate_framework.verification import CheckLedger

ledger = CheckLedger("P251/0071")
gold = (1+s.sqrt(5))/2
scale = s.sqrt(1+gold**2)
vertices = []
for left in (-1, 1):
    for right in (-1, 1):
        vertices.extend([s.Matrix([0, left, right*gold])/scale,
                         s.Matrix([left, right*gold, 0])/scale,
                         s.Matrix([right*gold, 0, left])/scale])
n = s.Matrix(s.symbols("x y z", real=True))
square = n.dot(n)
second = s.simplify(sum(v.dot(n)**2 for v in vertices)/12)
fourth = s.simplify(s.expand(sum(v.dot(n)**4 for v in vertices)/12))
ledger.check("actual twelve vertices have unit norm",
             len(vertices) == 12 and all(s.simplify(v.dot(v)-1) == 0 for v in vertices))
ledger.check("second directional moment is exactly isotropic",
             s.simplify(second-square/3) == 0)
ledger.check("fourth directional moment is exactly isotropic",
             s.simplify(fourth-square**2/5) == 0)
c0, c2, c4 = s.symbols("c0 c2 c4", real=True)
ledger.check("every source quartic spectral weight becomes radial",
             s.simplify(c0+c2*second+c4*fourth
                        -(c0+c2*square/3+c4*square**2/5)) == 0)
cross = s.Matrix([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
projector = (s.eye(3)-n*n.T+s.I*cross)/2
unit_relation = s.groebner([square-1], *n, domain=s.QQ)


def sphere_zero(value):
    real, imaginary = s.expand_complex(s.expand(value)).as_real_imag()
    return (unit_relation.reduce(s.expand(real))[1] == 0
            and unit_relation.reduce(s.expand(imaginary))[1] == 0)


ledger.check("helicity projector is idempotent with eigenvalue plus one",
             all(sphere_zero(v) for v in projector*projector-projector)
             and all(sphere_zero(v) for v in s.I*cross*projector-projector))
raw = s.Matrix([n[0]**2-1, n[0]*n[1]-s.I*n[2], n[0]*n[2]+s.I*n[1]])
ledger.check("actual source vector outer product has only its scalar axial weight",
             all(sphere_zero(v) for v in raw*raw.H-2*(1-n[0]**2)*projector))
ledger.check("using one axis fails isotropization despite a valid Beltrami field",
             s.expand(n[0]**4-square**2/5) != 0)
raise SystemExit(ledger.finish())
