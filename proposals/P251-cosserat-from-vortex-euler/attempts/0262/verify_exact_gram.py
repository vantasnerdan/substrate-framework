"""Exact oracle for the simultaneous Gram map; no PDE claim is tested here."""

import sympy as sp

from substrate_framework.euler_joint_normalizer import common_null_gram_lift


def main():
    n = 2
    eye = sp.eye(n)
    zero = sp.zeros(n)
    h_form = sp.diag(eye, -eye, zero, zero)
    o_form = sp.BlockMatrix(
        [[zero, zero, zero, zero], [zero, zero, zero, zero],
         [zero, zero, zero, eye], [zero, zero, -eye, zero]]
    ).as_explicit()
    v = sp.Matrix.vstack(eye, eye, eye, zero)
    w = sp.Matrix.vstack(eye / 2, -eye / 2, zero, zero)
    q = sp.Matrix.vstack(zero, zero, zero, eye)
    a, b, c, d, e, f, g, h = sp.symbols('a b c d e f g h', real=True)
    target_h = sp.Matrix([[a, b + sp.I*c], [b - sp.I*c, d]])
    target_o = sp.Matrix([[sp.I*e, f + sp.I*g], [-f + sp.I*g, sp.I*h]])
    scale = sp.symbols('L', positive=True)
    y = common_null_gram_lift(target_h, target_o, scale)
    assert y == scale*v + (w*target_h + q*target_o)/(2*scale)
    result_h = sp.simplify(y.conjugate().T*h_form*y)
    result_o = sp.simplify(y.conjugate().T*o_form*y)
    assert result_h == target_h
    assert result_o == target_o
    dual = w.row_join(q)
    assert dual.T*h_form*dual == sp.zeros(2*n)
    assert dual.T*o_form*dual == sp.zeros(2*n)
    # Removing either energy signature destroys the null-base mechanism.
    wrong_h = sp.diag(eye, eye, zero, zero)
    assert sp.simplify(y.conjugate().T*wrong_h*y-target_h) != sp.zeros(n)
    # The phase sign and half-normalization are independently load bearing.
    assert sp.simplify(y.conjugate().T*(-o_form)*y-target_o) != sp.zeros(n)
    wrong_y = v+w*target_h+q*target_o
    assert sp.simplify(wrong_y.conjugate().T*h_form*wrong_y) == 2*target_h
    # Imaginary diagonal rows are tested, not silently reduced to real skew.
    assert result_o[0, 0] == sp.I*e
    print('Exact complex joint Gram identities, common-null dual span, '
          'scale invariance and sign/signature/normalization mutations: PASS')


if __name__ == '__main__':
    main()
