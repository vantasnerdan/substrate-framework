"""Actual full Fourier Kelvin phase/energy versus integrated transport rows."""

import sympy as sp

from substrate_framework import euler_fourier as ef
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0223-complete-action")
    k = sp.Symbol("k", real=True, nonzero=True)
    psi = ef.add(ef.trig(2), ef.scale(ef.trig(1), sp.Rational(1, 100)))
    u = (psi, ef.trig(2, kind="sin"), ef.scale(ef.trig(1, kind="sin"), -sp.Rational(1, 100)))
    def trans(a):
        return ef.transport(u, (a, {}, {}))[0]

    def integ(a):
        return sp.cancel(a.get(ef.ZERO, 0))

    def pair(a, b):
        return integ(ef.mul(a, b))

    def hinv(a):
        return {q: c / (q[1]**2 + q[2]**2) for q, c in a.items() if q[1]**2 + q[2]**2}

    def grad_pair(a, b):
        return ef.add(*(ef.mul(ef.derivative(a, j), ef.derivative(b, j)) for j in (1, 2)))
    h_l = ef.trig(1, 2)
    s_l = h_l
    odd_product = ef.mul(ef.trig(1, 2, kind="sin"), ef.trig(2, kind="sin"))
    even_product = ef.mul(ef.trig(1, 2), ef.trig(2))
    h_r = ef.add(odd_product, even_product)
    s_r = ef.add(ef.scale(odd_product, 2), even_product)
    f_l, f_r = hinv(h_l), hinv(h_r)

    def generator(h, s, sign):
        def lift(field):
            return {(sign*k, q[1], q[2]): c for q, c in field.items()}
        f = hinv(h)
        return (lift(h),
                lift(ef.add(ef.scale(ef.derivative(s, 2), -1), ef.scale(ef.derivative(f, 1), sign*sp.I*k))),
                lift(ef.add(ef.derivative(s, 1), ef.scale(ef.derivative(f, 2), sign*sp.I*k))))

    xi_l, xi_r = generator(h_l, s_l, -1), generator(h_r, s_r, 1)
    velocities, energy, phase = ef.coadjoint_matrices(u, [xi_l, xi_r], beltrami_eigenvalue=-1)
    omega0 = pair(s_l, trans(s_r)) + pair(h_l, trans(s_r)) + pair(s_l, trans(h_r))
    omega1 = integ(ef.add(ef.mul(psi, ef.add(ef.mul(s_l, h_r), ef.mul(s_r, h_l))),
                         ef.scale(ef.mul(ef.add(s_l, h_l), grad_pair(psi, f_r)), -1),
                         ef.scale(ef.mul(ef.add(s_r, h_r), grad_pair(psi, f_l)), -1)))
    omega2 = pair(f_l, trans(f_r))
    ledger.check("the COMPLETE canonical KKS form equals all three derived transport coefficients",
                 sp.cancel(phase[0, 1] - omega0 - sp.I*k*omega1-k*k*omega2) == 0)
    ledger.check("the selected real columns expose a genuinely nonzero odd phase coefficient", omega1 != 0)
    ledger.check("the selected columns expose a genuinely nonzero second phase coefficient", omega2 != 0)
    left, right = velocities
    curl_l, curl_r = ef.curl(left)[0], ef.curl(right)[0]
    normal_j = ({}, ef.scale(right[2], -1), right[1])
    full_integrated = ef.inner(left, right) + pair(left[0], curl_r) + pair(curl_l, right[0]) + sp.I*k*ef.inner(left, normal_j)
    ledger.check("the full canonical energy equals the retained helicity/current formula at arbitrary k",
                 sp.cancel(energy[0, 1]-full_integrated) == 0)

    for sign, (h, s, f, velocity) in zip((-1, 1), ((h_l, s_l, f_l, left), (h_r, s_r, f_r, right)), strict=True):
        weighted = ({}, ef.mul(psi, ef.derivative(f, 1)), ef.mul(psi, ef.derivative(f, 2)))
        expected = ef.add(ef.scale(trans(ef.add(s, h)), -1), ef.scale(ef.divergence(weighted), sign*sp.I*k))
        lifted = {(sign*k, q[1], q[2]): c for q, c in expected.items()}
        difference = ef.add(ef.curl(velocity)[0], ef.scale(lifted, -1))
        ledger.check(f"full pressure normal curl uses only the true scalar transport rows {sign}",
                     all(sp.cancel(c) == 0 for c in difference.values()))

    expected_zero = -pair(trans(s_l), trans(s_r)) - pair(trans(s_l), trans(h_r)) - pair(trans(h_l), trans(s_r))
    expected_zero += pair(trans(ef.add(s_l, h_l)), hinv(trans(ef.add(s_r, h_r))))
    ledger.check("the zero-wave complete energy is the exact indefinite two-sector form",
                 sp.cancel(energy[0, 1].limit(k, 0)-expected_zero) == 0)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
