"""Exact audit of the full-history acoustic gain singularity in P251/0265.

This is distinct from P251/0250's bounded inverse of the already divided
quadratic coefficient q_N=Gamma_A/|K|^2.  It checks only the finite-dimensional
asymptotic mechanism and does not replace the physical source construction.
"""

import sympy as sp


k = sp.symbols("k", real=True)
h2, bo = sp.symbols("h2 bo", nonzero=True)
b1 = sp.symbols("b1")

# Scalar acoustic/optical representative of the triangular leading-carrier
# map in P251/0265.  The full proof has two transverse acoustic components and
# a three-component optical block; either acoustic singular value has this
# same exact order.
gain = sp.Matrix([[k**2 * h2, k * b1], [0, bo]])
inverse = sp.simplify(gain.inv())

assert sp.simplify(gain.det() - k**2 * h2 * bo) == 0
assert sp.simplify(k**2 * inverse[0, 0] - 1 / h2) == 0
assert sp.simplify(k**3 * sp.diff(inverse[0, 0], k) + 2 / h2) == 0

print("det(M) =", sp.factor(gain.det()))
print("M^{-1}_{AA} =", inverse[0, 0])
print("d/dk M^{-1}_{AA} =", sp.diff(inverse[0, 0], k))
print("VERDICT: no bounded C^3 gain inverse on a ball containing k=0")
print("REPAIR: freeze each nonzero k before the carrier/error diagonal")
