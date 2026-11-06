"""Product-rule connection reconstruction and exact square/rectangular ambiguities."""
import sympy as s
from .exact import Evidence, native


def square_product_rule(out):
    r = s.symbols("r", real=True)
    E = s.diag(-2, 3)
    M = s.Matrix([[1, 2+s.I], [-1+3*s.I, -2]])
    V = s.Matrix([[3, -2+4*s.I], [-2-4*s.I, -1]])
    T, h = s.eye(2)+r*M, E+r*V
    S, H = T.H*T, T.H*h*T
    B, A = S.diff(r).subs(r, 0), H.diff(r).subs(r, 0)
    D = native(A, B, E)
    out.identity("square_product_rule_Sprime", B-M.H-M)
    out.identity("square_product_rule_Hprime", A-V-M.H*E-E*M)
    out.identity("square_connection_reconstruction", D-(E*M-M*E)-V)
    out.exposed("native_equals_laboratory_derivative_mutation", D-V)
    out.exposed("Hermitian_part_equals_laboratory_derivative_mutation", (D+D.H)/2-V)
    out.identity("Hermitian_part_requires_connection_skew", (D+D.H)/2-V-(E*(M-M.H)-(M-M.H)*E)/2)
    symmetric = (M+M.H)/2
    Dsym = native(V+symmetric.H*E+E*symmetric, symmetric+symmetric.H, E)
    out.identity("symmetric_connection_suffices_for_Hermitian_part", (Dsym+Dsym.H)/2-V)
    return r, E, M, V, T, h, S, H


def square_nonuniqueness(out, data):
    r, E, M, V, T, h, S, H = data
    J = s.Matrix([[s.I, 1+2*s.I], [-1+2*s.I, -s.I]])
    U = (s.eye(2)+r*J/2)*(s.eye(2)-r*J/2).inv()
    out.identity("Cayley_unitary_exact_curve", U.H*U-s.eye(2))
    T2, h2 = U*T, U*h*U.H
    out.identity("two_embeddings_same_overlap_entire_curve", T2.H*T2-S)
    out.identity("two_embeddings_same_H_entire_curve", T2.H*h2*T2-H)
    V2 = h2.diff(r).subs(r, 0)
    out.identity("embedding_ambiguity_commutator", V2-V-(J*E-E*J))
    out.exposed("same_HS_and_all_derivatives_different_hprime", V2-V)
    out.identity("connection_resolves_square_ambiguity", V2+(E*(M+J)-(M+J)*E)-(V+E*M-M*E))
    return {"T1": str(T), "h1": str(h), "J": str(J),
            "T2_definition": "U*T1, U=(I+rJ/2)*(I-rJ/2)^-1", "h2_definition": "U*h1*Udagger",
            "hprime1": str(V), "hprime2": str(V2), "same_curves": "H and S exactly equal for real r"}


def rectangular_nonuniqueness(out):
    r = s.symbols("r", real=True)
    E = s.diag(-2, 3)
    T = s.Matrix([[1, 0], [0, 1], [r, s.I*r]])
    n = s.Matrix([-r, s.I*r, 1])
    v = s.Matrix([1+s.I, 2-s.I, 0])
    h1 = s.diag(-2+3*r, 3-r, 7)
    h2 = h1+n*v.H+v*n.H
    out.identity("normal_vector_orthogonal_entire_curve", T.H*n)
    out.identity("rectangular_same_T_same_H_entire_curve", T.H*h2*T-T.H*h1*T)
    Phi, Y = T.subs(r, 0), T.diff(r).subs(r, 0)
    V1 = Phi.H*h1.diff(r).subs(r, 0)*Phi
    V2 = Phi.H*h2.diff(r).subs(r, 0)*Phi
    out.exposed("even_known_T_without_ambient_h_not_sufficient", V2-V1)
    h = h2.subs(r, 0)
    W, Z = h*Phi-Phi*E, Phi.H*Y
    out.identity("generalized_eigenproblem_only_projected", Phi.H*W)
    out.exposed("ambient_eigenvector_assumption_false", W)
    H, S = T.H*h2*T, T.H*T
    D = native(H.diff(r).subs(r, 0), S.diff(r).subs(r, 0), E)
    residual_term = Y.H*W+W.H*Y
    out.identity("rectangular_full_product_rule", D-V2-(E*Z-Z*E)-residual_term)
    out.exposed("square_reconstruction_fails_for_rectangular_leakage", D-(E*Z-Z*E)-V2)
    out.exposed("normalization_metric_does_not_fix_leakage", residual_term)
    return {"T": str(T), "h1": str(h1), "h2": str(h2), "V1": str(V1), "V2": str(V2),
            "Z": str(Z), "W": str(W), "leakage_correction": str(residual_term)}


def run():
    out = Evidence()
    square = square_nonuniqueness(out, square_product_rule(out))
    rectangular = rectangular_nonuniqueness(out)
    return {"square_counterexample": square, "rectangular_counterexample": rectangular, "checks": out.rows}
