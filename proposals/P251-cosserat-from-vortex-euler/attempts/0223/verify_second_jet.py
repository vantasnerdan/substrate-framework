"""Full canonical scalar Kelvin second jet and local physical-current oracle."""

import sympy as sp

from substrate_framework import euler_fourier as ef
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0223-second-jet")
    k = sp.Symbol("k", real=True, nonzero=True)
    psi = ef.add(ef.trig(2), ef.scale(ef.trig(1), sp.Rational(1, 100)))
    u = (psi, ef.trig(2, kind="sin"), ef.scale(ef.trig(1, kind="sin"), -sp.Rational(1, 100)))

    def clean(field):
        return {q: c for q, value in field.items() if (c := sp.cancel(value)) != 0}

    def hinv(field):
        return {q: c / (q[1] ** 2 + q[2] ** 2) for q, c in field.items() if q[1] ** 2 + q[2] ** 2}

    def subtract(left, right):
        return clean(ef.add(left, ef.scale(right, -1)))

    def mf(field):
        return {q: c for q, c in field.items() if q != ef.ZERO}

    def trans(field):
        return ef.transport(u, (field, {}, {}))[0]

    def lift(field):
        return {(k, q[1], q[2]): c for q, c in field.items()}

    def jet(field, order):
        return clean({(0, q[1], q[2]): sp.diff(c, k, order).limit(k, 0) / sp.factorial(order) for q, c in field.items()})

    def weighted_grad(field):
        return ({}, ef.mul(psi, ef.derivative(field, 1)), ef.mul(psi, ef.derivative(field, 2)))

    pairs = [(ef.trig(1), ef.trig(2, 2)),
             (ef.mul(ef.trig(1, 2), ef.trig(2, 2)), ef.trig(1, 2))]
    for index, (h, s) in enumerate(pairs):
        f = hinv(h)
        xi = (lift(h),
              lift(ef.add(ef.scale(ef.derivative(s, 2), -1), ef.scale(ef.derivative(f, 1), sp.I * k))),
              lift(ef.add(ef.derivative(s, 1), ef.scale(ef.derivative(f, 2), sp.I * k))))
        actual = ef.material_kelvin_operator(u, xi)
        c = ef.divergence(({}, ef.mul(ef.add(s, h), ef.derivative(psi, 1)),
                           ef.mul(ef.add(s, h), ef.derivative(psi, 2))))
        sym_grad = ({}, *tuple(ef.add(*(ef.mul(ef.add(ef.derivative(u[i], j), ef.derivative(u[j], i)),
                                                ef.derivative(f, j)) for j in (1, 2))) for i in (1, 2)))
        csym = ef.curl(sym_grad)[0]
        expected_h = [ef.scale(trans(h), -1),
                      ef.scale(mf(ef.add(ef.mul(psi, subtract(s, h)), hinv(c))), sp.I),
                      ef.scale(hinv(trans(subtract(s, f))), -1)]
        expected_s = [subtract(hinv(trans(ef.add(s, h))), trans(s)),
                      ef.scale(hinv(ef.add(ef.divergence(weighted_grad(subtract(f, s))), csym)), -sp.I),
                      ef.scale(hinv(trans(f)), -1)]
        for order in range(3):
            got_h = jet(actual[0], order)
            normal = ({}, jet(actual[1], order), jet(actual[2], order))
            got_s = ef.scale(hinv(ef.curl(normal)[0]), -1)
            ledger.check(f"full pressure gives exact h Taylor coefficient {order}, preparation {index}",
                         not subtract(got_h, expected_h[order]))
            ledger.check(f"full normal Kelvin curl gives exact s Taylor coefficient {order}, preparation {index}",
                         not subtract(got_s, expected_s[order]))
        ledger.check(f"the actual complete material rate remains Bloch solenoidal {index}",
                     not clean(ef.divergence(actual)))

    # Exposing pointwise identity: no integration or fictitious periodic
    # coordinate weight is used here. A compact physical tag subsequently
    # removes the proved divergence in the source.
    a, b = sp.symbols("a b", real=True)
    ps = sp.cos(b) + sp.cos(a) / 100
    vv = sp.Matrix([-sp.diff(ps, b), sp.diff(ps, a)])
    ss = sp.cos(2*a) * sp.cos(2*b)
    ff = sp.cos(a)
    hh = -sp.diff(ff, a, 2) - sp.diff(ff, b, 2)
    def grad(q):
        return sp.Matrix([sp.diff(q, a), sp.diff(q, b)])

    def div(vec):
        return sp.diff(vec[0], a) + sp.diff(vec[1], b)

    def cross(left, right):
        return left[0]*right[1] - left[1]*right[0]
    rr = sp.Matrix([a, b])
    xx = sp.Matrix([-sp.diff(ss, b), sp.diff(ss, a)]) + sp.I*k*grad(ff)
    chi = ps**2
    primitive = ps**3/3
    dp = rr.dot(grad(ps))
    ww = sp.Matrix(sp.symbols("w_a w_b"))
    rate = ww - xx.jacobian([a, b])*vv + vv.jacobian([a, b])*xx - sp.I*k*ps*xx
    delta = -xx.dot(grad(chi))
    sdens = chi*cross(rr, ww) + delta*dp
    gdens_rate = chi*cross(rr, rate)
    aux = 2*primitive - chi*dp
    flux = sp.I*k*(chi*ps*cross(rr, xx) + aux*hh)
    boundary = div(chi*vv*cross(rr, xx) + aux*xx)
    ledger.check("the complete moving spin/current identity holds POINTWISE including axial flux",
                 sp.expand(sdens-gdens_rate-flux-boundary) == 0)
    ledger.check("the current uses the true nonzero normal Bloch divergence",
                 sp.expand(div(xx)+sp.I*k*hh) == 0)
    ledger.check("discarding the retained axial flux changes this physical local identity",
                 sp.expand(flux) != 0)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
