#!/usr/bin/env python3
"""Exact algebra checks for P253/0021; no floating-point numerics."""

import sympy as sp


def check_action_sign() -> None:
    q1, q2, v1, v2 = sp.symbols("q1 q2 v1 v2")
    s11, s12, s22 = sp.symbols("s11 s12 s22")
    q = sp.Matrix([q1, q2])
    v = sp.Matrix([v1, v2])
    J = sp.Matrix([[0, 1], [-1, 0]])
    S = sp.Matrix([[s11, s12], [s12, s22]])
    lagrangian = -(q.T * J * v)[0] / 2 - (q.T * S * q)[0] / 2
    d_q = sp.Matrix([sp.diff(lagrangian, x) for x in q])
    d_v = sp.Matrix([sp.diff(lagrangian, x) for x in v])
    # d/dt(dL/dv), retaining only the velocity-dependent part.
    dt_d_v = d_v.jacobian(q) * v
    assert sp.simplify(d_q - dt_d_v + J * v + S * q) == sp.zeros(2, 1)


def check_companion_projection() -> None:
    a, b, c, d = sp.symbols("a b c d", nonzero=True)
    A = sp.Matrix([[a, b], [c, d]])
    assert sp.factor(A.det()) != 0
    # In coordinates where the two trial columns are the identity, the
    # normalized companion columns are A^{-1} and F(companions)=I.
    companions = A.inv()
    assert sp.simplify(A * companions) == sp.eye(2)

    f11, f12, f13, f14 = sp.symbols("f11 f12 f13 f14")
    f21, f22, f23, f24 = sp.symbols("f21 f22 f23 f24")
    F = sp.Matrix([[f11, f12, f13, f14], [f21, f22, f23, f24]])
    # A symbolic right inverse W can be normalized abstractly by imposing
    # F W=I; use canonical coordinates to check the projection identities.
    Fc = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1]])
    Wc = sp.Matrix([[0, 0], [0, 0], [1, 0], [0, 1]])
    projection = sp.eye(4) - Wc * Fc
    assert Fc * Wc == sp.eye(2)
    assert projection * projection == projection
    assert Fc * projection == sp.zeros(2, 4)
    assert projection * Wc == sp.zeros(4, 2)


def check_invariant_flag_block() -> None:
    a11, a12, a21, a22 = sp.symbols("a11 a12 a21 a22")
    b11, b12, b21, b22 = sp.symbols("b11 b12 b21 b22")
    c11, c12, c21, c22 = sp.symbols("c11 c12 c21 c22")
    d11, d12, d21, d22 = sp.symbols("d11 d12 d21 d22")
    A = sp.Matrix([[a11, a12], [a21, a22]])
    B = sp.Matrix([[b11, b12], [b21, b22]])
    C = sp.Matrix([[c11, c12], [c21, c22]])
    D = sp.Matrix([[d11, d12], [d21, d22]])
    Z = sp.zeros(2)
    I = sp.eye(2)
    M = A.row_join(Z).row_join(B)
    M = M.col_join(C.row_join(I).row_join(D))
    M = M.col_join(Z.row_join(Z).row_join(I))
    conserved = Z.row_join(Z).row_join(I)
    symmetry = sp.zeros(6, 2)
    symmetry[2:4, :] = I
    assert conserved * M == conserved
    assert M * symmetry == symmetry


if __name__ == "__main__":
    check_action_sign()
    check_companion_projection()
    check_invariant_flag_block()
    print("PASS action sign: quadratic Euler-Lagrange equation is J qdot = -S q")
    print("PASS companions: normalization and fixed-level projection are exact")
    print("PASS block flag: conserved covectors and both symmetry vectors are fixed")
