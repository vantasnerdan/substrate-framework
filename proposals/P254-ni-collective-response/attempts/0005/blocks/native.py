"""Derive native adjoint and eigenvector-tangent constraints in an S-orthonormal frame."""
import sympy as s
from .exact import Evidence, native


def symbolic_identities(out):
    e0, e1, a0, a1, b0, b1, ar, ai, br, bi = s.symbols(
        "e0 e1 a0 a1 b0 b1 ar ai br bi", real=True)
    E = s.diag(e0, e1)
    A = s.Matrix([[a0, ar+s.I*ai], [ar-s.I*ai, a1]])
    B = s.Matrix([[b0, br+s.I*bi], [br-s.I*bi, b1]])
    D = native(A, B, E)
    out.identity("native_adjoint_D_minus_Ddag_equals_EB_minus_BE", D-D.H-(E*B-B*E))
    F = A-(E*B+B*E)/2
    out.identity("Hermitian_part_is_arithmetic_energy_symmetrization", (D+D.H)/2-F)
    out.identity("Hermitian_part_adjoint", F-F.H)
    out.identity("native_diagonal_real", D[0, 0]-s.conjugate(D[0, 0]))
    K = s.Matrix([[-b0/2, D[0, 1]/(e1-e0)],
                  [D[1, 0]/(e0-e1), -b1/2]])
    Ep = s.diag(D[0, 0], D[1, 1])
    out.identity("differentiated_generalized_eigenproblem", A+E*K-B*E-K*E-Ep)
    out.identity("differentiated_S_normalization", K+K.H+B)
    out.identity("degenerate_energy_block_native_Hermitian", (D-D.H).subs(e1, e0))
    return {"D": str(D), "Hermitian_part": str(F), "K_non_degenerate": str(K),
            "assumptions": "real eigenvalues; e0 != e1 only for K division; arbitrary imaginary diagonal gauge may be added to K"}


def asymmetric_example(out):
    C = s.Matrix([[1, s.I], [1+s.I, 2]])
    inverse = C.inv()
    E = s.diag(-2, 3)
    S = inverse.H*inverse
    H = inverse.H*E*inverse
    A = s.Matrix([[2, 1+3*s.I], [1-3*s.I, -1]])
    B = s.Matrix([[1, 2-s.I], [2+s.I, -2]])
    Hp, Sp = inverse.H*A*inverse, inverse.H*B*inverse
    out.identity("nonorthogonal_reference_normalized", C.H*S*C-s.eye(2))
    out.identity("generalized_reference_eigenproblem", H*C-S*C*E)
    out.exposed("ordinary_normalization_is_not_S_normalization", C.H*C-s.eye(2))
    D = native(C.H*Hp*C, C.H*Sp*C, E)
    out.exposed("initial_to_final_energy_mutation", D-(A-E*B))
    out.exposed("omitted_overlap_mutation", D-A)
    out.exposed("native_is_not_generally_Hermitian", D-D.H)
    out.exposed("transpose_is_not_adjoint", D.T-D.H)
    out.exposed("noncommuting_overlap_energy", E*B-B*E)
    return {"C": str(C), "S": str(S), "H": str(H), "A": str(A), "B": str(B),
            "E": str(E), "D": str(D), "S_positive_reason": "C invertible, S=C^-dagger C^-1"}


def run():
    out = Evidence()
    identities = symbolic_identities(out)
    example = asymmetric_example(out)
    return {"identities": identities, "example": example, "checks": out.rows}
