"""New exact Langreth identity and discrete two-sided invariant, not old proofs."""
from pathlib import Path
import json
import sys
import sympy as s


def blocks():
    # Exact noncommuting finite matrices encode the retarded and advanced blocks.
    g = s.Matrix([[-s.I, s.Rational(1, 5)], [0, -s.I]])
    sigma = s.Matrix([[s.Rational(1, 7), s.I/9], [s.Rational(1, 11), -s.Rational(1, 8)]])
    q = s.I*s.Matrix([[s.Rational(2, 5), s.Rational(1, 6)+s.I/7],
                      [s.Rational(1, 6)-s.I/7, s.Rational(3, 5)]])
    k = s.I*s.Matrix([[s.Rational(1, 8), s.I/13], [-s.I/13, s.Rational(1, 3)]])
    one = s.eye(2)
    r = (one-g*sigma).inv()*g  # never inverse g, including singular boson g0
    a = r*sigma
    dressed = (one+a)*q*(one+a.adjoint())+r*k*r.adjoint()
    residual = (one-g*sigma)*dressed-q*(one+sigma.adjoint()*r.adjoint())-g*k*r.adjoint()
    langreth = all(s.cancel(v)==0 for v in residual)
    # Real row quadrature weights enter independently on the two sides.
    w = s.diag(s.Rational(1, 2), s.Rational(-1, 12))
    t = one+a*w
    discrete = t*q*t.adjoint()+r*w*k*w*r.adjoint()
    exchange = all(s.cancel(v)==0 for v in discrete+discrete.adjoint())
    wrong = (one+w*a)*q*t.adjoint()+r*w*k*w*r.adjoint()
    order_sensitive = any(s.cancel(v)!=0 for v in wrong+wrong.adjoint())
    omitted_initial_sensitive = any(s.cancel(v)!=0 for v in t*q*t.adjoint())
    omitted_noise_sensitive = any(s.cancel(v)!=0 for v in r*w*k*w*r.adjoint())
    assert langreth and exchange and order_sensitive and omitted_initial_sensitive and omitted_noise_sensitive
    return {'langreth_left_equation_exact':langreth, 'real_quadrature_two_sided_exchange_exact':exchange,
            'wrong_noncommuting_order_detected':order_sensitive,
            'initial_source_omission_detected':omitted_initial_sensitive,
            'noise_omission_detected':omitted_noise_sensitive,
            'statistics':'Same real-time Langreth plus sign for fermion and Cartesian boson; statistics enters mixed branch only, identically zero here.',
            'scope':'Exact encoded block algebra only; no automatic equality of nested NESSi finite-step Dyson and two-sided discretizations.'}

if __name__ == '__main__':
    Path(sys.argv[1]).write_text(json.dumps(blocks(),indent=2)+'\n')
