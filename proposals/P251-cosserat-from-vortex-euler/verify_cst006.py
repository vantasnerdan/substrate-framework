"""C-CST-006 verifier: orientation-ergodic contrast closure.

The exact result is about the *coherent signed response*: under an independent
uniform frame phase, the first-order frame map and hence any response linear in
that map average to zero. The quadratic fluctuation energy does not average to
zero, because ``<L.T L> = I``. Keeping both facts explicit prevents a vanishing
mean from being misreported as vanishing microscopic fluctuation energy.

The no-Cosserat macroscopic limit therefore additionally assumes the declared
ergodic closure: no coherent microrotation field is retained after the signed
response is averaged. Under that explicit closure the remaining constitutive
state is the translational Navier--Cauchy sector. The Monte Carlo check measures
the same signed first-moment observable as the exact calculation.
"""

import sys

import numpy as np
import sympy as sp

from substrate_framework.verification import CheckLedger


def phase_rotation(phi: sp.Expr) -> sp.Matrix:
    return sp.Matrix([[sp.cos(phi), sp.sin(phi)], [-sp.sin(phi), sp.cos(phi)]])


def check_phase_average(ledger: CheckLedger) -> None:
    phi = sp.Symbol("phi", real=True)
    rotation = phase_rotation(phi)
    average = sp.simplify(sp.integrate(rotation, (phi, 0, 2 * sp.pi)) / (2 * sp.pi))
    average_gram = sp.simplify(
        sp.integrate(rotation.T * rotation, (phi, 0, 2 * sp.pi)) / (2 * sp.pi)
    )
    ledger.check(
        "uniform phase has zero coherent frame map",
        average == sp.zeros(2, 2),
        "<L(phi)> = 0",
    )
    ledger.check(
        "uniform phase retains quadratic frame fluctuations",
        average_gram == sp.eye(2),
        "<L(phi).T L(phi)> = I, so zero mean alone does not erase energy",
    )


def check_factorization(ledger: CheckLedger) -> None:
    theta, phi = sp.symbols("theta phi", real=True)
    n = sp.Matrix(
        [sp.sin(theta) * sp.cos(phi), sp.sin(theta) * sp.sin(phi), sp.cos(theta)]
    )
    mean_n = sp.Matrix(
        [
            sp.simplify(
                sp.integrate(
                    sp.integrate(component * sp.sin(theta), (phi, 0, 2 * sp.pi)),
                    (theta, 0, sp.pi),
                )
                / (4 * sp.pi)
            )
            for component in n
        ]
    )
    mean_t = mean_n.copy()
    joint_independent = sp.simplify(mean_n * mean_t.T)
    ledger.check(
        "isotropic vector means vanish by exact sphere integration",
        mean_n == sp.zeros(3, 1) and mean_t == sp.zeros(3, 1),
        "<n> = <t> = 0",
    )
    ledger.check(
        "independent joint first moment factorizes to zero",
        joint_independent == sp.zeros(3, 3),
        "<n_i t_j> = <n_i><t_j> = 0",
    )


def check_coherent_response_closure(ledger: CheckLedger) -> None:
    phi = sp.Symbol("phi", real=True)
    response = sp.Matrix(sp.symbols("r1:3", real=True))
    average_map = sp.simplify(
        sp.integrate(phase_rotation(phi), (phi, 0, 2 * sp.pi)) / (2 * sp.pi)
    )
    mean_response = sp.simplify(average_map * response)
    ledger.check(
        "ergodic signed-response closure has zero coherent couple response",
        mean_response == sp.zeros(2, 1),
        "the macroscopic closure retains no coherent Phi response",
    )


def sample_signed_response(
    seed: int, count: int, bias: float = 0.0
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    gaussian = rng.standard_normal((count, 3))
    directions = gaussian / np.linalg.norm(gaussian, axis=1, keepdims=True)
    phase = rng.uniform(0.0, 2.0 * np.pi, count)
    weight = directions[:, 0] ** 2 - directions[:, 1] ** 2
    samples = np.column_stack((weight * np.cos(phase) + bias, weight * np.sin(phase)))
    return samples, weight


def check_simulation(ledger: CheckLedger) -> None:
    count = 200_000
    samples, _weight = sample_signed_response(20260904, count)
    means = np.mean(samples, axis=0)
    standard_errors = np.std(samples, axis=0, ddof=1) / np.sqrt(count)
    z_scores = np.abs(means) / standard_errors
    print(
        "MC_EVIDENCE "
        f"N={count} seed=20260904 means={means.tolist()} "
        f"standard_errors={standard_errors.tolist()} z_scores={z_scores.tolist()}"
    )
    ledger.check(
        "Monte Carlo signed response agrees with the zero first moment",
        bool(np.all(z_scores <= 5.0)),
        f"N={count}, seed=20260904, z_scores={z_scores.tolist()}, float64",
    )


def check_mutations(ledger: CheckLedger) -> None:
    count = 200_000
    baseline, _weight = sample_signed_response(20260904, count)
    standard_error = float(np.std(baseline[:, 0], ddof=1) / np.sqrt(count))
    biased, _weight = sample_signed_response(20260904, count, bias=7.0 * standard_error)
    biased_z = abs(float(np.mean(biased[:, 0]))) / float(
        np.std(biased[:, 0], ddof=1) / np.sqrt(count)
    )
    ledger.check(
        "M1 biased signed response is rejected",
        biased_z > 5.0,
        f"biased z-score={biased_z:.6f}",
    )

    phi, epsilon = sp.symbols("phi epsilon", real=True)
    density = (1 + 2 * epsilon * sp.cos(phi)) / (2 * sp.pi)
    biased_map = sp.simplify(
        sp.integrate(phase_rotation(phi) * density, (phi, 0, 2 * sp.pi))
    )
    ledger.check(
        "M2 phase-locked distribution produces a nonzero coherent map",
        biased_map != sp.zeros(2, 2),
        f"<L>_biased = {biased_map}",
    )


def main() -> int:
    ledger = CheckLedger("C-CST-006")
    check_phase_average(ledger)
    check_factorization(ledger)
    check_coherent_response_closure(ledger)
    check_simulation(ledger)
    check_mutations(ledger)
    return int(ledger.finish())


if __name__ == "__main__":
    sys.exit(main())
