"""Check energy, gauge, normalization, and moving coordinate-frame transformations."""
import sympy as s
from .exact import Evidence, native


def energy_origin(out, E, A, B):
    shift, slope = s.symbols("shift slope", real=True)
    D = native(A, B, E)
    shifted = native(A+shift*B+slope*s.eye(2), B, E+shift*s.eye(2))
    out.identity("energy_origin_shift_derivative", shifted-D-slope*s.eye(2))
    out.exposed("wrong_energy_shift_missing_shift_Sprime", native(A, B, E+5*s.eye(2))-D)
    out.exposed("wrong_energy_shift_sign", native(A-5*B, B, E+5*s.eye(2))-D)
    C = s.Matrix([[1, s.I], [1+s.I, 2]])
    S = C.inv().H*C.inv()
    H = C.inv().H*E*C.inv()
    out.identity("generalized_shift_is_H_plus_shift_S", (H+5*S)*C-S*C*(E+5*s.eye(2)))
    out.exposed("H_plus_shift_identity_is_not_generalized_shift", (H+5*s.eye(2))*C-S*C*(E+5*s.eye(2)))


def reference_gauges(out, E, A, B):
    D = native(A, B, E)
    P = s.diag((3+4*s.I)/5, s.I)
    out.identity("phase_gauge_unitarity", P.H*P-s.eye(2))
    out.identity("phase_covariance", native(P.H*A*P, P.H*B*P, E)-P.H*D*P)
    out.exposed("missing_bra_conjugation", P.T*D*P-P.H*D*P)
    out.exposed("entries_are_not_gauge_invariants", P.H*D*P-D)
    out.identity("entry_absolute_square_gauge_invariant", s.Abs((P.H*D*P)[0, 1])**2-s.Abs(D[0, 1])**2)
    U = s.Matrix([[1, s.I], [s.I, 1]])/s.sqrt(2)
    Ed = 3*s.eye(2)
    out.identity("degenerate_subspace_unitary_covariance", native(U.H*A*U, U.H*B*U, Ed)-U.H*native(A, B, Ed)*U)
    L = s.diag(2, 1+s.I)
    out.identity("unnormalized_raw_congruence", native(L.H*A*L, L.H*B*L, E)-L.H*D*L)
    out.exposed("normalization_changes_raw_entries", L.H*D*L-D)
    out.identity("normalization_recovery", L.inv().H*(L.H*D*L)*L.inv()-D)


def moving_frames(out, E, A, B):
    r = s.symbols("r", real=True)
    G = s.Matrix([[1, 2+s.I], [-1+3*s.I, -2]])
    U = s.eye(2)+r*G
    H = E+r*A
    S = s.eye(2)+r*B
    Hp = (U.H*H*U).diff(r).subs(r, 0)
    Sp = (U.H*S*U).diff(r).subs(r, 0)
    D, transformed = native(A, B, E), native(Hp, Sp, E)
    out.identity("moving_frame_native_commutator", transformed-D-(E*G-G*E))
    out.exposed("native_not_invariant_under_moving_frame", transformed-D)
    F = (D+D.H)/2
    Fnew = (transformed+transformed.H)/2
    out.identity("Hermitian_part_connection_dependence", Fnew-F-(E*(G-G.H)-(G-G.H)*E)/2)
    out.exposed("arbitrary_Hermitization_does_not_remove_frame", Fnew-F)
    return {"G": str(G), "moving_native_change": str(transformed-D),
            "Hermitian_part_change": str(Fnew-F)}


def orthonormalization_tangent(out, E, A, B):
    r = s.symbols("r", real=True)
    O = s.eye(2)-r*B/2
    H, S = E+r*A, s.eye(2)+r*B
    out.identity("symmetric_orthonormalization_metric_first_jet", (O.H*S*O).diff(r).subs(r, 0))
    F = (native(A, B, E)+native(A, B, E).H)/2
    out.identity("symmetric_orthonormalization_operator_derivative", (O.H*H*O).diff(r).subs(r, 0)-F)
    J = s.Matrix([[s.I, 1+2*s.I], [-1+2*s.I, -s.I]])
    O2 = O*(s.eye(2)+r*J)
    out.identity("other_orthonormalization_same_metric_first_jet", (O2.H*S*O2).diff(r).subs(r, 0))
    derivative2 = (O2.H*H*O2).diff(r).subs(r, 0)
    out.identity("orthonormalization_unitary_connection_shift", derivative2-F-(E*J-J*E))
    out.exposed("orthonormal_metric_does_not_select_operator_derivative", derivative2-F)


def run():
    out = Evidence()
    E = s.diag(-2, 3)
    A = s.Matrix([[2, 1+3*s.I], [1-3*s.I, -1]])
    B = s.Matrix([[1, 2-s.I], [2+s.I, -2]])
    energy_origin(out, E, A, B)
    reference_gauges(out, E, A, B)
    frames = moving_frames(out, E, A, B)
    orthonormalization_tangent(out, E, A, B)
    return {"moving_frame": frames, "checks": out.rows}
