"""Exact exterior-degenerate quadratic census for P249 attempt 0002."""

from __future__ import annotations

from itertools import combinations_with_replacement

import sympy as sp

from substrate_framework.m5_covariant_action import (
    MINKOWSKI_MOSTLY_PLUS,
    projected_spatial_ldg_potential,
)
from substrate_framework.m5_fluctuation_spectrum import (
    certify_projector_variation,
    timelike_rotation_kinetic_metric,
)
from substrate_framework.verification import CheckLedger


PAIRS = tuple(combinations_with_replacement(range(4), 2))


def _basis() -> tuple[sp.Matrix, ...]:
    result = []
    for row, column in PAIRS:
        matrix = sp.zeros(4)
        matrix[row, column] = 1
        matrix[column, row] = 1
        result.append(sp.Matrix(matrix))
    return tuple(result)


def _coordinates(matrix: sp.Matrix) -> sp.Matrix:
    return sp.Matrix([matrix[row, column] for row, column in PAIRS])


def main() -> int:
    ledger = CheckLedger("P249-O2-0002")
    eta = sp.Matrix(MINKOWSKI_MOSTLY_PLUS)
    vacuum = sp.diag(-4, 1, 0, 0)
    projector_t = sp.diag(1, 0, 0, 0)
    epsilon = sp.Symbol("epsilon", real=True)
    kappa = sp.Symbol("kappa", positive=True)
    basis = _basis()
    ledger.check("symmetric channel order exact", PAIRS == ((0, 0), (0, 1), (0, 2), (0, 3), (1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3)))

    # Potential orbit invariance makes the boost rows/columns exactly zero.
    # Every remaining channel stays block diagonal with the same spectral Pt,
    # so its complete quadratic form is evaluated directly from the canonical
    # projected M5.17 potential.
    boost_indices = {1, 2, 3}

    def second_coefficient(variation: sp.Matrix) -> sp.Expr:
        potential = projected_spatial_ldg_potential(
            vacuum + epsilon * variation,
            projector_t,
            beta=1,
            scale=1,
            timelike_target_eigenvalue=4,
            timelike_stiffness=1,
            metric_covariant=eta,
        )
        return sp.expand(potential).coeff(epsilon, 2)

    stiffness = sp.zeros(10)
    rest_indices = [index for index in range(10) if index not in boost_indices]
    for index in rest_indices:
        stiffness[index, index] = sp.cancel(second_coefficient(basis[index]))
    for left_offset, left in enumerate(rest_indices):
        for right in rest_indices[left_offset + 1 :]:
            cross = (
                second_coefficient(basis[left] + basis[right])
                - stiffness[left, left]
                - stiffness[right, right]
            ) / 2
            stiffness[left, right] = sp.cancel(cross)
            stiffness[right, left] = sp.cancel(cross)
    expected_k = sp.diag(1, 0, 0, 0, sp.Rational(5, 2), 0, 0, sp.Rational(3, 2), 3, sp.Rational(3, 2))
    ledger.check("full projected-potential stiffness exact", stiffness == expected_k)
    ledger.check("stiffness rank exact", stiffness.rank() == 5)
    ledger.check("stiffness is positive semidefinite", all(value >= 0 for value in stiffness.diagonal()))

    canonical_basis, kinetic, _ = timelike_rotation_kinetic_metric(
        targets=(4, 1, 0, 0),
        projector_stiffness=kappa,
        basis=basis,
    )
    ledger.check("canonical kinetic basis retained", canonical_basis == basis)
    expected_g = sp.diag(0, kappa / 9, kappa / 16, kappa / 16, 0, 0, 0, 0, 0, 0)
    ledger.check("projector-current kinetic metric exact", kinetic == expected_g)
    ledger.check("kinetic rank exact", kinetic.rank() == 3)
    ledger.check(
        "kinetic image positive for kappa positive",
        all(kinetic[index, index].is_positive for index in (1, 2, 3)),
    )
    ledger.check(
        "degenerate spatial targets retain unique timelike-projector variation",
        certify_projector_variation(targets=(4, 1, 0, 0), basis=basis) == [0] * 10,
    )

    common_kernel = stiffness.nullspace()
    kinetic_kernel = kinetic.nullspace()
    common = sp.Matrix.vstack(stiffness, kinetic).nullspace()
    ledger.check("stiffness nullity exact", len(common_kernel) == 5)
    ledger.check("kinetic nullity exact", len(kinetic_kernel) == 7)
    ledger.check("common inert kernel dimension exact", len(common) == 2)
    expected_common = [sp.eye(10)[:, 5], sp.eye(10)[:, 6]]
    ledger.check(
        "common kernel is director shear pair",
        common == expected_common,
    )

    generator = sp.zeros(4)
    generator[2, 3] = -1
    generator[3, 2] = 1
    representation = sp.zeros(10)
    for column, element in enumerate(basis):
        representation[:, column] = _coordinates(generator * element + element * generator.T)
    ledger.check("SO(2) preserves stiffness", representation.T * stiffness + stiffness * representation == sp.zeros(10))
    ledger.check("SO(2) preserves kinetic metric", representation.T * kinetic + kinetic * representation == sp.zeros(10))

    split_diagonal = sp.eye(10)[:, 7] - sp.eye(10)[:, 9]
    split_off_diagonal = sp.eye(10)[:, 8]
    ledger.check(
        "clock split is weight two",
        representation * split_diagonal == 2 * split_off_diagonal
        and representation * split_off_diagonal == -2 * split_diagonal,
    )
    ledger.check(
        "clock split stiffness isotropic",
        (split_diagonal.T * stiffness * split_diagonal)[0] == 3
        and (split_off_diagonal.T * stiffness * split_off_diagonal)[0] == 3,
    )
    ledger.check(
        "clock split has zero quadratic kinetic metric",
        kinetic * split_diagonal == sp.zeros(10, 1)
        and kinetic * split_off_diagonal == sp.zeros(10, 1),
    )

    wave_number_squared, frequency_squared = sp.symbols("k2 omega2", nonnegative=True)
    principal = stiffness + wave_number_squared * kinetic
    for index in (1, 2, 3):
        ledger.check(
            f"physical branch {PAIRS[index]} is exactly massless",
            sp.factor(principal[index, index] - frequency_squared * kinetic[index, index]).subs(frequency_squared, wave_number_squared) == 0,
        )
    ledger.check("radiation threshold is exactly zero", sp.solve(frequency_squared - wave_number_squared, frequency_squared) == [wave_number_squared])

    # Exact perturbative order check for a noncommuting derivative pair:
    # D_mu=epsilon W_mu makes F=epsilon^2[W_mu,W_nu], hence F^2 starts at 4.
    w0 = basis[4]
    w1 = basis[5]
    commutator = (epsilon * w0) * eta * (epsilon * w1) - (epsilon * w1) * eta * (epsilon * w0)
    curvature_norm = sp.expand(sp.trace(commutator.T * commutator))
    ledger.check("curvature quadratic coefficient vanishes exactly", curvature_norm.coeff(epsilon, 2) == 0)
    ledger.check("curvature first appears at fourth order", curvature_norm.coeff(epsilon, 4) != 0)

    mutated_targets = (4, 1, sp.Rational(1, 10), 0)
    _, mutated_g, _ = timelike_rotation_kinetic_metric(
        targets=mutated_targets, projector_stiffness=kappa, basis=basis
    )
    ledger.check("tangent-degeneracy mutation changes the kinetic record", mutated_g != kinetic)

    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
