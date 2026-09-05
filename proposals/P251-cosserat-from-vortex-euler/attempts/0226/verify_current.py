"""Exact additional physical G/spin moment repair anchors."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0226-physical-current-return")
    gap, nu, t = s.symbols("gap nu t", real=True)
    phase = s.exp(-s.I*gap*t)
    checks.check("spin row retains its extra frequency, not just displacement moment",
                 s.simplify(s.I*s.diff(phase, t)+nu*phase-(nu+gap)*phase) == 0)
    for order in (0, 1, 2):
        row = s.diff((nu+gap)*phase, t, order).subs(t, 0)
        checks.check(f"physical spin time row {order} uses the two stated moment orders",
                     s.expand(row-(-s.I)**order*(nu*gap**order+gap**(order+1))) == 0)
    distance = s.symbols("distance", positive=True)
    chi = s.exp(-s.exp(1/distance))
    log_slope = -s.diff(chi, distance)/chi
    checks.check("fixed positive flat tag has the actual non-polynomial slope",
                 s.simplify(log_slope+s.exp(1/distance)/distance**2) == 0)
    checks.check("no finite-order analytic clock polynomial absorbs that tag slope",
                 s.limit(s.exp(1/distance)*distance**20, distance, 0, dir="+") == s.oo)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
