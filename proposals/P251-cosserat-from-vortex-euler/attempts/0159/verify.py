"""Thin replay of the new packet APIs at the frozen C-CST-011 parameters."""

import sympy as s

from substrate_framework.euler_core_packet import (
    common_circle_angular_rule,
    common_circle_moment_weights,
    gaussian_carrier_filter,
    laguerre_packet_angle,
    packet_material_moment_rows,
)
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0159-packet-api")
    w, x = s.symbols("w x", real=True)
    p, lam, om, delta, length, c = s.symbols("p lambda Omega delta L c", positive=True)
    angle = laguerre_packet_angle(8, w)
    gamma = -2*om+(18+s.diff(angle, w).subs(w, 1))*s.sqrt(2)*om*s.sqrt(lam/p)
    curvature = s.simplify((p*p*s.diff(gamma*gamma, p, 2)).subs(p, lam/delta**2))
    ledger.check("physical curvature is derived by actual carrier differentiation",
                 s.simplify(s.limit(curvature/delta, delta, 0)-s.sqrt(2)*om**2/3) == 0)
    rows = packet_material_moment_rows(8, 7, x)
    jets = s.Matrix([[s.diff(row, x, j).subs(x, 0) for j in range(len(rows))]
                     for row in rows])
    determinant = jets.det(method="domain-ge")
    ledger.check("full twenty-row radial/reference construction has nonzero minor", determinant != 0)
    filtered = gaussian_carrier_filter(length, c*length, p, p)
    ledger.check("actual finite marker filter retains the nonunit curvature factor",
                 s.simplify(filtered.carrier_factor-1/(1+c*c)) == 0
                 and s.simplify(filtered.carrier_factor-1) != 0)
    ledger.check("whole finite-packet measure is not an assigned fiber length",
                 filtered.plancherel_weight == s.sqrt(s.pi)*length)
    rule = common_circle_angular_rule(3)
    actual_weights = common_circle_moment_weights(rule.nodes)
    ledger.check("common-circle actual moment solve preserves the positive rule",
                 actual_weights.weights == rule.weights
                 and all(weight.is_positive for weight in actual_weights.weights))
    ledger.check("six angular moments agree with the exact arcsine measure",
                 all(s.simplify(sum(a*b**j for a, b in zip(rule.weights, rule.nodes))
                                -s.binomial(2*j, j)/4**j) == 0 for j in range(6)))
    print(f"actual radial/reference determinant={determinant}")
    print(f"derived physical curvature={curvature}")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
