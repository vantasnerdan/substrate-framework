"""Exact full-Fourier initial-time stress jets, without spectral truncation."""

import sympy as s

from substrate_framework.euler_fourier import (
    ZERO, add, curl, derivative, divergence, leray, mul, scale, transport, trig,
)
from substrate_framework.verification import CheckLedger


def vector_add(*vectors):
    return tuple(add(*(v[i] for v in vectors)) for i in range(3))


def vector_scale(vector, coefficient):
    return tuple(scale(v, coefficient) for v in vector)


def equal(left, right):
    return not any(vector_add(left, vector_scale(right, -1)))


def mean_product(left, right):
    return mul(left, right).get(ZERO, 0)


def main():
    ledger = CheckLedger("P251-0175-observed-stress-time-jets")
    d = s.Symbol("d", real=True)
    # Actual two-wave member of the negative-helicity ABC family.
    u = (add(trig(2), scale(trig(1, kind="sin"), d)),
         trig(2, kind="sin"), scale(trig(1), d))
    p = scale(add(*(mul(v, v) for v in u)), -s.Rational(1, 2))
    dp = tuple(derivative(p, j) for j in range(3))
    ledger.check("actual background is constant-curl", equal(curl(u), vector_scale(u, -1)))
    ledger.check("actual stationary Euler pressure", not any(vector_add(transport(u, u), dp)))

    def c(field):
        return vector_add(vector_scale(transport(u, field), -1), transport(field, u))

    def linear_euler(field):
        return vector_scale(leray(vector_add(transport(u, field), transport(field, u))), -1)

    def jacobi_potential(field):
        hp = tuple(add(*(mul(derivative(dp[i], j), field[j]) for j in range(3)))
                   for i in range(3))
        return leray(vector_add(transport(u, transport(u, field)), hp))

    # T[i,l,j,m] multiplies kappa_j*kappa_m*D_l before the slow projector.
    tensors = [s.MutableDenseNDimArray.zeros(3, 3, 3, 3) for _ in range(3)]
    for m in range(3):
        for ell in range(3):
            seed = tuple(add(u[m] if i == ell else {}, u[ell] if i == m else {})
                         for i in range(3))
            chi1 = vector_scale(leray(seed), -1)
            forcing = leray(tuple(add(dp[m] if i == ell else {}, dp[ell] if i == m else {})
                                  for i in range(3)))
            w1 = vector_add(linear_euler(chi1), forcing)
            chi2 = vector_add(c(chi1), w1)
            w2 = linear_euler(w1)
            chi3 = vector_add(c(chi2), w2)
            direct2 = vector_add(vector_scale(leray(transport(u, chi1)), -2), forcing)
            direct3 = vector_add(vector_scale(leray(transport(u, chi2)), -2),
                                 vector_scale(jacobi_potential(chi1), -1))
            ledger.check(f"Euler/Lin equals Jacobi for cell ({m},{ell})",
                         equal(chi2, direct2) and equal(chi3, direct3)
                         and not divergence(chi2) and not divergence(chi3))
            zero = ({}, {}, {})
            for order, (chi, rate) in enumerate(((zero, chi1), (chi1, chi2), (chi2, chi3))):
                for i in range(3):
                    for j in range(3):
                        value = (mean_product(u[j], rate[i]) + mean_product(u[i], rate[j])
                                 + mean_product(dp[j], chi[i]) + mean_product(dp[i], chi[j]))
                        if order == 0 and i == ell:
                            value += mean_product(u[j], u[m])
                        tensors[order][i, ell, j, m] = s.expand(value)

    def isotropic_transverse(tensor):
        # First average the transverse polarization, then the common direction.
        # E k_j k_m=delta_jm/3; E k_i k_l k_j k_m=(three pairings)/15.
        diagonal = sum(tensor[i, i, j, j] for i in range(3) for j in range(3))
        crossed = sum(tensor[i, ell, i, ell] + tensor[i, ell, ell, i]
                      for i in range(3) for ell in range(3))
        return s.factor(s.Rational(2, 15)*diagonal - s.Rational(1, 30)*crossed)

    observed = [isotropic_transverse(t) for t in tensors]
    print("Whole-field isotropic physical R_D(0), R_D'(0), R_D''(0):", observed, flush=True)
    ledger.check("one-wave reference initial acceleration", observed[0].subs(d, 0) == -s.Rational(2, 15))
    ledger.check("one-wave reference has no temporal stress derivatives",
                 all(value.subs(d, 0) == 0 for value in observed[1:]))
    ledger.check("time-reversal odd derivative vanishes in this paired family", observed[1] == 0)
    print("Exact cancellation candidate:", "survives this jet" if observed[2] == 0 else "refuted by nonzero second time derivative")
    print("No statement about all stationary backgrounds, higher jets, or parent completion.")
    raise SystemExit(ledger.finish())


if __name__ == "__main__":
    main()
