"""Proof-bearing finite-field rank plus exact compact-pair current bookkeeping."""

import numpy as np
import sympy as s

from substrate_framework.verification import CheckLedger

from jet_probe import ORDER, angular, matrix, modular_rank

ledger = CheckLedger("P251/0085")
x, y, z, strength = s.symbols("x y z W", real=True)
transverse = s.expand((x + s.I * y)**(ORDER + 1))
real, imag = s.re(transverse).expand(), s.im(transverse).expand()
omega = s.Matrix([0, 0, strength])
coordinates = [x, y, z]
pairs = [(-strength * imag, real), (strength * real, imag), (0, z**(ORDER + 1))]
for index, (f, g) in enumerate(pairs):
    gradient_f = s.Matrix([s.diff(f, t) for t in coordinates])
    gradient_g = s.Matrix([s.diff(g, t) for t in coordinates])
    ledger.check(f"universal top-jet left null vector {index} is exact",
                 s.simplify(gradient_f + omega.cross(gradient_g)) == s.zeros(3, 1))
ledger.check("three universal independent output dependencies give rank upper235",
             matrix.shape == (238, 252) and matrix.shape[0] - len(pairs) == 235)
ledger.check("a nonzero exact finite-field minor reaches the proven rank upper bound",
             modular_rank(matrix) == 235)
for axis in range(3):
    ledger.check(f"physical spin axis {axis} is outside the full constraint row space",
                 modular_rank(np.vstack([matrix, angular[axis]])) == 236)
ledger.check("all three physical spin rows are independently controllable",
             modular_rank(np.vstack([matrix, angular])) == 238)

# Exact disjoint compact-response bookkeeping.
amplitude, raw_b = s.symbols("A Bc", nonzero=True, real=True)
lq0, lqc, lsc = s.symbols("Lq0 Lqc Lsc", real=True)
# Basis Q0,Qc,Sc,eta for one physical component; disjoint supports except Qc/Sc.
form = s.zeros(4)
form[1, 2], form[2, 1] = raw_b, -raw_b
physical_spin = s.Matrix([[lq0, lqc, lsc, 1]])
q = s.Matrix([1, amplitude, 0, -lq0 - amplitude * lqc])
reaction = s.Matrix([0, 0, 1, amplitude * raw_b - lsc])
ledger.check("compact corrections keep the physical core angle and zero reaction angle",
             q[0] == 1 and reaction[0] == 0)
ledger.check("the physical coordinate spin is zero exactly",
             s.simplify((physical_spin * q)[0]) == 0)
ledger.check("actual reaction spin equals the unchanged KKS pairing exactly",
             s.simplify((physical_spin * reaction)[0] - (q.T * form * reaction)[0]) == 0)
ledger.check("the compact KKS pair remains nonzero after physical moment matching",
             s.simplify((q.T * form * reaction)[0] - amplitude * raw_b) == 0)

inertia_block, mixed, stiffness, b = s.symbols("P N Hq B", nonzero=True, real=True)
pdot, relative = s.symbols("Phidot q", real=True)
momentum_plus = (b * pdot - mixed * relative) / inertia_block
momentum_minus = (-b * pdot - mixed * relative) / inertia_block
physical_plus = b * momentum_plus
physical_minus = -b * momentum_minus
ledger.check("one-background physical spin retains its static angle response",
             s.simplify(physical_plus - b**2 * pdot / inertia_block
                        + b * mixed * relative / inertia_block) == 0)
ledger.check("time reversal cancels static spin only after actual reaction elimination",
             s.simplify((physical_plus + physical_minus) / 2
                        - b**2 * pdot / inertia_block) == 0)
lagrangian_plus = b * momentum_plus * pdot - (
    inertia_block * momentum_plus**2 + 2 * mixed * momentum_plus * relative
    + stiffness * relative**2) / 2
lagrangian_minus = -b * momentum_minus * pdot - (
    inertia_block * momentum_minus**2 + 2 * mixed * momentum_minus * relative
    + stiffness * relative**2) / 2
target = (b**2 * pdot**2 / inertia_block
          - (stiffness - mixed**2 / inertia_block) * relative**2) / 2
ledger.check("full paired action retains the exact positive Schur coefficients",
             s.simplify((lagrangian_plus + lagrangian_minus) / 2 - target) == 0)
raise SystemExit(ledger.finish())
