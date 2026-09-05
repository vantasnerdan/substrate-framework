"""Continue the first actual stress jet to the fourth time derivative."""

import sympy as s
from verify import equal, mean_product, vector_add, vector_scale

from substrate_framework.euler_fourier import (
    add, curl, derivative, divergence, leray, mul, scale, transport, trig,
)
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0175-higher-time-continuation")
    # Exact rational amplitudes; all product harmonics are retained.
    u = (add(trig(2), trig(1, kind="sin")), trig(2, kind="sin"), trig(1))
    p = scale(add(*(mul(v, v) for v in u)), -s.Rational(1, 2))
    dp = tuple(derivative(p, j) for j in range(3))
    ledger.check("same stationary constant-curl field", equal(curl(u), vector_scale(u, -1))
                 and not any(vector_add(transport(u, u), dp)))

    def c(field):
        return vector_add(vector_scale(transport(u, field), -1), transport(field, u))

    def linear_euler(field):
        return vector_scale(leray(vector_add(transport(u, field), transport(field, u))), -1)

    # The exact isotropic contractions reduce to Z/5-Y/15 for symmetric
    # first-cell rows chi^{ml}=chi^{lm}; pressure has the same contraction.
    z = [s.S.Zero for _ in range(5)]
    y = [s.S.Zero for _ in range(5)]
    trace_cells = [({}, {}, {}) for _ in range(6)]
    for m in range(3):
        for ell in range(m, 3):
            seed = tuple(add(u[m] if i == ell else {}, u[ell] if i == m else {})
                         for i in range(3))
            w = vector_scale(leray(seed), -1)
            forcing = leray(tuple(add(dp[m] if i == ell else {}, dp[ell] if i == m else {})
                                  for i in range(3)))
            cells = [({}, {}, {}), w]
            for n in range(1, 5):
                w = vector_add(linear_euler(w), forcing) if n == 1 else linear_euler(w)
                cells.append(vector_add(c(cells[-1]), w))
            ledger.check(f"full Fourier Lin derivatives stay solenoidal ({m},{ell})",
                         all(not divergence(field) for field in cells))
            for order in range(5):
                chi, rate = cells[order], cells[order+1]
                z[order] += mean_product(u[m], rate[ell]) + mean_product(dp[m], chi[ell])
                if m != ell:
                    z[order] += mean_product(u[ell], rate[m]) + mean_product(dp[ell], chi[m])
                if m == ell:
                    y[order] += sum(mean_product(u[i], rate[i]) + mean_product(dp[i], chi[i])
                                    for i in range(3))
            if m == ell:
                trace_cells = [vector_add(trace_cells[n], cells[n]) for n in range(6)]
            print("Completed full derivatives for cell", m, ell, flush=True)
    ledger.check("exact trace first cell is minus twice t*u",
                 equal(trace_cells[1], vector_scale(u, -2))
                 and all(not any(field) for field in trace_cells[2:]))
    energy = sum(mean_product(v, v) for v in u)
    observed = [s.factor(z[n]/5-y[n]/15+(energy/3 if n == 0 else 0)) for n in range(5)]
    print("R_D^(n)(0), n=0..4:", observed, flush=True)
    ledger.check("first independent contraction agrees with earlier tensor route",
                 observed[:3] == [-s.Rational(4, 15), 0, 0])
    print("Exact temporal cancellation:", "survives through fourth jet" if all(v == 0 for v in observed[1:])
          else "refuted for this two-wave field at the displayed nonzero jet")
    raise SystemExit(ledger.finish())


if __name__ == "__main__":
    main()
