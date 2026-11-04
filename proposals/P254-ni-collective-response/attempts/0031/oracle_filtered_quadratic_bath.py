"""Owning A31 exact-small-block/SciPy oracle. Run only when authorized by parent.

No material data, subspace module, stationary state or FDT enters this oracle.
Truth routes differ at the load-bearing integration risk: adaptive scalar direct
quadrature, and Gauss-Legendre direct matrix exponentials with refinement, versus
the production spectral/Bessel/Legendre integral. Matrix exponentials are cached
per quadrature node and reused across ALL modes/probes and mutation checks.
Independent 16-dimensional CAR/Fock traces test ordered nonstationary covariance.
Exact here describes the finite four-orbital model, not rigorous roundoff bounds.
The program makes no claim of material convergence or full moving/RPA closure.
"""

from __future__ import annotations

import json

import numpy as np
from numpy.polynomial.chebyshev import chebval
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.linalg import expm, solve

from filtered_quadratic_bath import integrate_filtered_operators


TOLERANCE = 3e-10


def _error(actual, expected):
    return float(np.linalg.norm(actual - expected) / max(1.0, np.linalg.norm(expected)))


def _agree(name, actual, expected, records, tolerance=TOLERANCE):
    error = _error(actual, expected)
    records[name] = error
    if not np.isfinite(error) or error > tolerance:
        raise AssertionError(f"{name}: relative/absolute scaled error {error} > {tolerance}")


def _exposes(name, wrong, expected, records):
    error = _error(wrong, expected)
    records[name] = error
    if not np.isfinite(error) or error < 1e-5:
        raise AssertionError(f"insensitive fixture for {name}: scaled difference {error}")


def _complex_quad(function, left, right):
    real, real_error = quad(lambda t: float(np.real(function(t))), left, right,
                            epsabs=2e-13, epsrel=2e-13, limit=100)
    imag, imag_error = quad(lambda t: float(np.imag(function(t))), left, right,
                            epsabs=2e-13, epsrel=2e-13, limit=100)
    return real + 1j * imag, real_error + imag_error


def _scalar_gap_oracle(records):
    # Exact degeneracy, signed gaps, and gaps far below inverse-gap-recurrence
    # stability, with every degree separately visible and complex coefficients.
    energies = np.array([0.0, 0.0, 1e-13, 2.3])
    bounds = np.array([[-0.8, 0.35], [0.35, 1.7], [1.7, 1.83]])
    coefficients = np.zeros((3, 8, 1, 8), complex)
    for s in range(3):
        for degree in range(8):
            coefficients[s, degree, 0, degree] = (1 + 0.17 * s + 0.03 * degree) + 1j * (0.2 - 0.09 * s + 0.07 * degree)
    forces = np.ones((1, 4, 4), complex)
    force_adjoints = (0.73 + 0.19j) * forces
    result = integrate_filtered_operators(np.diag(energies), forces, force_adjoints,
                                          bounds, coefficients)
    expected = np.zeros((8, 4, 4), complex)
    expected_adj = np.zeros_like(expected)
    wrong_negative_parity = np.zeros_like(expected)
    maximum_quad_error = 0.0
    # Cache repeated scalar gaps (the exact degenerate block repeats entries).
    cache = {}
    for i, ei in enumerate(energies):
        for j, ej in enumerate(energies):
            gap = float(ei - ej)
            if gap not in cache:
                values = np.zeros(8, complex)
                adj_values = np.zeros(8, complex)
                bad_values = np.zeros(8, complex)
                for s, (left, right) in enumerate(bounds):
                    half = (right - left) / 2
                    midpoint = left + half
                    for probe in range(8):
                        poly = coefficients[s, :, 0, probe]
                        def value(t):
                            return chebval((t - midpoint) / half, poly)
                        integral, estimate = _complex_quad(
                            lambda t: value(t) * np.exp(1j * gap * t), left, right)
                        adj_integral, adj_estimate = _complex_quad(
                            lambda t: value(t).conjugate() * np.exp(1j * gap * t), left, right)
                        # This is exactly the missing j_n(-z)=(-1)^n*j_n(z)
                        # parity defect, while retaining the correct midpoint phase.
                        bad_integral, _ = _complex_quad(
                            lambda t: value(t) * np.exp(1j * (gap * midpoint + abs(gap) * (t - midpoint))), left, right)
                        values[probe] += integral
                        adj_values[probe] += adj_integral * (0.73 + 0.19j)
                        bad_values[probe] += bad_integral
                        maximum_quad_error = max(maximum_quad_error, estimate, adj_estimate)
                cache[gap] = values, adj_values, bad_values
            expected[:, i, j], expected_adj[:, i, j], wrong_negative_parity[:, i, j] = cache[gap]
    _agree("scalar_all_degrees_zero_negative_close_gaps", result["operators"], expected, records)
    _agree("scalar_complex_supplied_adjoint_weights", result["adjoint_operators"], expected_adj, records)
    _exposes("mutation_negative_Bessel_parity", wrong_negative_parity, expected, records)
    zero_gap_wrong = expected.copy()
    zero_gap_wrong[:, energies[:, None] == energies[None, :]] = 0
    _exposes("mutation_zero_gap_dropped", zero_gap_wrong, expected, records)
    # Analytic T_0 zero-gap normalization: integral is interval length, not 2*length.
    analytic_zero = sum((right - left) * coefficients[s, 0, 0, 0]
                        for s, (left, right) in enumerate(bounds))
    _agree("exact_constant_zero_gap_limit", result["operators"][0, 0, 1], analytic_zero, records)
    records["scalar_quad_maximum_reported_error"] = maximum_quad_error
    if result["diagnostics"]["exact_zero_energy_gap_count"] != 6:
        raise AssertionError("raw exact degeneracy was not retained")


def _fixture():
    index = np.arange(4)
    # An explicit complex unitary, with no production eigensolver dependence.
    unitary = np.exp(2j * np.pi * index[:, None] * index[None, :] / 4) / 2
    energies = np.array([-0.7, -0.7, -0.7 + 1e-11, 1.3])
    h = (unitary * energies) @ unitary.conj().T
    row, col = np.indices((4, 4))
    f = np.stack([
        (0.2 + 0.13 * row - 0.09 * col) + 1j * (0.11 + 0.07 * row * col - 0.16 * col),
        (0.17 - 0.08 * row * col + 0.03 * col) + 1j * (-0.12 + 0.19 * row + 0.04 * col),
    ])
    fa = f.conj().transpose(0, 2, 1).copy()
    # A deliberately independent supplied-adjoint input exposes its replacement.
    fa[0, 0, 2] += 0.14 - 0.08j
    fa[1, 3, 1] -= 0.09 + 0.12j
    bounds = np.array([[-0.6, 0.2], [0.2, 3.4], [3.4, 3.57]])
    s, d, m, p = np.indices((3, 8, 2, 3))
    coefficients = ((0.21 + 0.06 * s + 0.07 * m - 0.03 * p) * (-1.) ** (d + p)
                    + 1j * (0.12 - 0.03 * s + 0.05 * m + 0.08 * p)) / (d + 1.)
    coefficients[:, 7] *= 3  # Do not let a truncated degree-six solver pass.
    return h, f, fa, bounds, coefficients


def _matrix_quadrature(h, forces, supplied_adjoints, bounds, coefficients, order,
                       mutations=False):
    """One expm per unique node; all probes/modes share that exact evolution.

    Inverse evolution is obtained by a solve, not a second expm or a Hermiticity
    assumption. Thus the raw non-Hermitian fixture is tested as well.
    """
    nodes, weights = leggauss(order)
    probes = coefficients.shape[3]
    output = np.zeros((probes, h.shape[0], h.shape[0]), complex)
    adjoint = np.zeros_like(output)
    mutants = {name: np.zeros_like(output) for name in (
        "gap_sign", "missing_interval_phase", "missing_jacobian",
        "unconjugated_adjoint_weights", "manufactured_force_adjoints", "last_interval_dropped",
        "degree_seven_dropped")}
    identity = np.eye(h.shape[0])
    for segment, (left, right) in enumerate(bounds):
        half = (right - left) / 2
        midpoint = left + half
        if mutations:
            undo_midpoint = expm(-1j * h * midpoint)
            redo_midpoint = solve(undo_midpoint, identity)
        for x, weight in zip(nodes, weights):
            t = midpoint + half * x
            forward = expm(1j * h * t)
            backward = solve(forward, identity)
            evolved = forward @ forces @ backward
            adjoint_evolved = backward.conj().T @ supplied_adjoints @ forward.conj().T
            scalar = chebval(x, coefficients[segment])  # [mode,probe], Clenshaw
            value = np.einsum("mp,mij->pij", scalar, evolved)
            adj_value = np.einsum("mp,mij->pij", scalar.conj(), adjoint_evolved)
            output += (half * weight) * value
            adjoint += (half * weight) * adj_value
            if mutations:
                mutants["gap_sign"] += half * weight * np.einsum(
                    "mp,mij->pij", scalar, backward @ forces @ forward)
                centered_forward = forward @ undo_midpoint
                centered_backward = redo_midpoint @ backward
                mutants["missing_interval_phase"] += half * weight * np.einsum(
                    "mp,mij->pij", scalar, centered_forward @ forces @ centered_backward)
                mutants["missing_jacobian"] += weight * value
                mutants["unconjugated_adjoint_weights"] += half * weight * np.einsum(
                    "mp,mij->pij", scalar, adjoint_evolved)
                mutants["manufactured_force_adjoints"] += half * weight * np.einsum(
                    "mp,mij->pij", scalar.conj(),
                    backward.conj().T @ forces.conj().transpose(0, 2, 1) @ forward.conj().T)
                if segment != len(bounds) - 1:
                    mutants["last_interval_dropped"] += half * weight * value
                low_degree_scalar = chebval(x, coefficients[segment, :-1])
                mutants["degree_seven_dropped"] += half * weight * np.einsum(
                    "mp,mij->pij", low_degree_scalar, evolved)
    return output, adjoint, mutants


def _annihilators(size):
    """Entrywise CAR matrices with explicit occupation-bit/Jordan-Wigner sign."""
    operators = np.zeros((size, 2 ** size, 2 ** size), complex)
    for mode in range(size):
        for ket in range(2 ** size):
            if ket & (1 << mode):
                bra = ket ^ (1 << mode)
                sign = (-1) ** ((ket & ((1 << mode) - 1)).bit_count())
                operators[mode, bra, ket] = sign
    return operators


def _lift(matrix, annihilators):
    result = np.zeros(annihilators.shape[1:], complex)
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            result += matrix[i, j] * (annihilators[i].conj().T @ annihilators[j])
    return result


def _ordered_car_oracle(h, operators, adjoints, reference, reference_adjoints,
                        records):
    # Fixed physical Gaussian state, not diagonal in H. Its eigenvalues are
    # fixture occupation probabilities, not a temperature or a fitted refill.
    state_rotation = expm(1j * np.array([
        [0.2, 0.31 + 0.13j, -0.09j, 0.17],
        [0.31 - 0.13j, -0.1, 0.23, 0.11j],
        [0.09j, 0.23, 0.07, -0.19 + 0.08j],
        [0.17, -0.11j, -0.19 - 0.08j, 0.3],
    ]))
    occupations = np.array([0.12, 0.35, 0.68, 0.9])
    state = (state_rotation * occupations) @ state_rotation.conj().T
    commutator = float(np.linalg.norm(h @ state - state @ h))
    records["noncommuting_state_commutator_fro"] = commutator
    if commutator < 0.05:
        raise AssertionError("state fixture accidentally stationary")
    annihilators = _annihilators(4)
    identity = np.eye(16)
    for i in range(4):
        for j in range(4):
            _agree(f"CAR_{i}_{j}", annihilators[i] @ annihilators[j].conj().T
                   + annihilators[j].conj().T @ annihilators[i],
                   identity if i == j else np.zeros_like(identity), records)
    creators = np.einsum("ik,ifg->kfg", state_rotation,
                         annihilators.conj().transpose(0, 2, 1))
    density = np.zeros((16, 16), complex)
    for mask in range(16):
        ket = np.zeros(16, complex)
        ket[0] = 1
        probability = 1.0
        for mode, occupation in enumerate(occupations):
            occupied = bool(mask & (1 << mode))
            probability *= occupation if occupied else 1 - occupation
            if occupied:
                ket = creators[mode] @ ket
        density += probability * np.outer(ket, ket.conj())
    _agree("Fock_density_trace", np.trace(density), 1.0, records)
    reconstructed = np.array([[np.trace(density @ annihilators[j].conj().T @ annihilators[i])
                               for j in range(4)] for i in range(4)])
    _agree("Fock_original_one_particle_state", reconstructed, state, records)
    lifted = [_lift(a, annihilators) for a in reference]
    lifted_adj = [_lift(a, annihilators) for a in reference_adjoints]
    expected = np.empty((len(operators), len(operators)), complex)
    expected_reverse = np.empty_like(expected)
    calculated = np.empty_like(expected)
    reversed_calculated = np.empty_like(expected)
    for p in range(len(operators)):
        for q in range(len(operators)):
            means = np.trace(density @ lifted[p]) * np.trace(density @ lifted_adj[q])
            expected[p, q] = np.trace(density @ lifted[p] @ lifted_adj[q]) - means
            expected_reverse[p, q] = np.trace(density @ lifted_adj[q] @ lifted[p]) - means
            calculated[p, q] = np.trace(state @ operators[p] @ (np.eye(4) - state) @ adjoints[q])
            reversed_calculated[p, q] = np.trace(state @ adjoints[q] @ (np.eye(4) - state) @ operators[p])
    _agree("ordered_CAR_covariance_direct_Fock_trace", calculated, expected, records)
    _agree("reverse_ordered_CAR_covariance_direct_Fock_trace", reversed_calculated, expected_reverse, records)
    _exposes("mutation_covariance_order_reversed", reversed_calculated, expected, records)
    # The tempting spectral-population shortcut changes the original state.
    energies, vectors = np.linalg.eigh(h)  # ONLY oracle mutation, never production H
    pinched_state = (vectors * np.diag(vectors.conj().T @ state @ vectors)) @ vectors.conj().T
    stationary_wrong = np.array([[np.trace(pinched_state @ a @ (np.eye(4) - pinched_state) @ b)
                                 for b in adjoints] for a in operators])
    _exposes("mutation_stationary_population_replacement", stationary_wrong, expected, records)


def run_oracle():
    """Execute only tiny synthetic fixtures; raises on agreement OR sensitivity failure."""
    records = {}
    _scalar_gap_oracle(records)
    h, forces, supplied_adjoints, bounds, coefficients = _fixture()
    result = integrate_filtered_operators(h, forces, supplied_adjoints, bounds, coefficients)
    coarse, coarse_adj, _ = _matrix_quadrature(h, forces, supplied_adjoints, bounds, coefficients, 24)
    reference, reference_adj, mutants = _matrix_quadrature(
        h, forces, supplied_adjoints, bounds, coefficients, 48, mutations=True)
    _agree("matrix_quadrature_refinement", coarse, reference, records)
    _agree("matrix_adjoint_quadrature_refinement", coarse_adj, reference_adj, records)
    _agree("complex_H_direct_expm_integral", result["operators"], reference, records)
    _agree("independent_supplied_adjoint_direct_expm_integral", result["adjoint_operators"], reference_adj, records)
    for name, wrong in mutants.items():
        target = reference_adj if name in ("unconjugated_adjoint_weights", "manufactured_force_adjoints") else reference
        _exposes("mutation_" + name, wrong, target, records)
    # Physical adjoints recover an actual adjoint, without manufacturing the
    # deliberately distinct input in the preceding test.
    true_adjoints = forces.conj().transpose(0, 2, 1)
    physical = integrate_filtered_operators(h, forces, true_adjoints, bounds, coefficients)
    _agree("true_adjoint_integral", physical["adjoint_operators"],
           physical["operators"].conj().transpose(0, 2, 1), records)
    # A separate direct-quadrature physical-adjoint leg is implied entrywise by
    # adjointing its independent direct-expm operator, not the production result.
    _ordered_car_oracle(h, physical["operators"], physical["adjoint_operators"],
                        reference, reference.conj().transpose(0, 2, 1), records)
    raw_h = h.copy()
    raw_h[0, 1] += 0.025 + 0.017j
    raw_h[2, 2] += 0.006j
    raw = integrate_filtered_operators(raw_h, forces, supplied_adjoints, bounds, coefficients)
    raw_coarse, raw_coarse_adj, _ = _matrix_quadrature(raw_h, forces, supplied_adjoints, bounds, coefficients, 24)
    raw_reference, raw_reference_adj, _ = _matrix_quadrature(raw_h, forces, supplied_adjoints, bounds, coefficients, 48)
    _agree("raw_H_matrix_quadrature_refinement", raw_coarse, raw_reference, records)
    _agree("raw_H_adjoint_quadrature_refinement", raw_coarse_adj, raw_reference_adj, records)
    _agree("raw_nonhermitian_H_not_repaired", raw["operators"], raw_reference, records)
    _agree("raw_Hdagger_supplied_adjoint_evolution", raw["adjoint_operators"], raw_reference_adj, records)
    _agree("raw_H_hermiticity_diagnostic", raw["diagnostics"]["hamiltonian_hermiticity_absolute_fro"],
           np.linalg.norm(raw_h - raw_h.conj().T), records)
    repaired = integrate_filtered_operators((raw_h + raw_h.conj().T) / 2, forces,
                                            supplied_adjoints, bounds, coefficients)
    _exposes("mutation_silent_Hermitization", repaired["operators"], raw_reference, records)
    if not raw["diagnostics"]["all_outputs_finite"]:
        raise AssertionError("nonfinite raw-H outputs")
    return {
        "schema": "P254-A31-filtered-quadratic-bath-owning-oracle-v1",
        "scope": "audited synthetic finite-block integration and ordered nonstationary CAR covariance only",
        "approximation": "frozen-initial-H selected-mode finite-basis comparator, NOT moving/RPA closure",
        "truth_routes": ["adaptive scalar direct quadrature", "refined Gauss-Legendre direct matrix exponential quadrature", "independent occupation-bit CAR/Fock traces"],
        "tolerance": TOLERANCE,
        "checks": records,
        "raw_H_diagnostics": raw["diagnostics"],
        "material_data_read": False,
        "not_established": ["rigorous enclosure", "material basis convergence", "physical loss", "full campaign"],
    }


if __name__ == "__main__":
    print(json.dumps(run_oracle(), indent=2, sort_keys=True, allow_nan=False))
