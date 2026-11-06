"""Use independent generalized eigensolves and centered refinement, not formula echoes."""
import numpy as np
from . import numeric_model as model


def maxnorm(matrix):
    return float(np.max(np.abs(matrix)))


def reference():
    T, h, H, S = model.curves(0.)
    e, C, Phi = model.solve(0.)
    Tp, Hp, Sp = model.analytic_jets()
    assert np.min(np.linalg.eigvalsh(S)) > .1
    assert np.min(np.diff(e)) > .1
    assert maxnorm(H@C-S@C@np.diag(e)) < 1.e-12
    assert maxnorm(C.conj().T@S@C-np.eye(2)) < 1.e-12
    A, B = C.conj().T@Hp@C, C.conj().T@Sp@C
    D = A-B@np.diag(e)
    Z = Phi.conj().T@Tp@C
    reconstructed = D-(np.diag(e)@Z-Z@np.diag(e))
    direct = Phi.conj().T@model.V@Phi
    assert maxnorm(reconstructed-direct) < 1.e-12
    return e, C, Phi, S, A, B, D, reconstructed


def finite_difference(step, ref):
    e, C, Phi, S, A, B, D, reconstructed = ref
    ep, Cp, _ = model.solve(step)
    em, Cm, _ = model.solve(-step)
    Cp, Cm = model.align(Cp, C, S), model.align(Cm, C, S)
    K = C.conj().T@S@((Cp-Cm)/(2*step))
    ep_fd = (ep-em)/(2*step)
    native_fd = np.diag(ep_fd).astype(complex)
    for j in range(2):
        for i in range(2):
            if i != j:
                native_fd[j, i] = (e[i]-e[j])*K[j, i]
    Tp, hp, Hp, Sp = model.curves(step)
    Tm, hm, Hm, Sm = model.curves(-step)
    physical_K = Phi.conj().T@((Tp@Cp-Tm@Cm)/(2*step))
    physical_fd = np.diag(ep_fd).astype(complex)
    for j in range(2):
        for i in range(2):
            if i != j:
                physical_fd[j, i] = (e[i]-e[j])*physical_K[j, i]
    direct_h_fd = Phi.conj().T@((hp-hm)/(2*step))@Phi
    return {
        "step": step,
        "native_from_eigenvectors_error": maxnorm(native_fd-D),
        "eigenvalue_derivative_error": maxnorm(ep_fd-np.real(np.diag(D))),
        "normalization_derivative_error": maxnorm(K+K.conj().T+B),
        "laboratory_from_physical_eigenvectors_error": maxnorm(physical_fd-reconstructed),
        "laboratory_from_direct_curve_error": maxnorm(direct_h_fd-reconstructed),
        "Hprime_direct_curve_error": maxnorm(C.conj().T@((Hp-Hm)/(2*step))@C-A),
        "Sprime_direct_curve_error": maxnorm(C.conj().T@((Sp-Sm)/(2*step))@C-B),
        "mutant_final_energy_error": maxnorm(native_fd-(A-np.diag(e)@B)),
        "mutant_omitted_overlap_error": maxnorm(native_fd-A),
        "mutant_omitted_connection_error": maxnorm(physical_fd-D),
        "mutant_arbitrary_Hermitization_error": maxnorm(physical_fd-(D+D.conj().T)/2),
    }


def check_refinement(rows):
    metrics = [key for key in rows[0] if key.endswith("_error") and not key.startswith("mutant")]
    orders = {}
    for key in metrics:
        errors = [row[key] for row in rows]
        assert all(np.isfinite(errors)), (key, errors)
        assert all(a > b for a, b in zip(errors, errors[1:])), (key, errors)
        orders[key] = [float(np.log2(a/b)) for a, b in zip(errors, errors[1:])]
        assert all(1.7 < order < 2.3 for order in orders[key][-2:]), (key, orders[key])
        assert errors[-1] < .002, (key, errors[-1])
    for key in rows[-1]:
        if key.startswith("mutant"):
            assert rows[-1][key] > .01, (key, rows[-1][key])
    return orders


def run():
    ref = reference()
    rows = [finite_difference(step, ref) for step in (.04, .02, .01, .005)]
    orders = check_refinement(rows)
    return {"scope": "float64 numeric_evidence, not exact proof or rigorous enclosure",
            "solver": "scipy.linalg.eigh generalized Hermitian driver gvd; scipy.linalg.expm; centered differences",
            "norm": "maximum absolute matrix entry; dimensionless",
            "inputs": model.inputs(), "reference_energies": ref[0].tolist(),
            "refinement": rows, "observed_orders": orders}
