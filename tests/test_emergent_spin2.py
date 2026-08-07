import pytest
import sympy as sp

from substrate_framework.dimensional_sine_gordon import dimensional_sine_gordon_coefficients
from substrate_framework.emergent_spin2 import (
    collective_tensor_constraint_matrix,
    fierz_pauli_gauge_ledger,
    fierz_pauli_gauge_residual,
    sine_gordon_spin_two_coupling_ledger,
    sine_gordon_spin_two_ledger,
    spin_two_mode_ledger,
    tensor_mode_count_from_constraints,
)


def test_gauge_invariance_solves_a_unique_fierz_pauli_ray() -> None:
    ledger = fierz_pauli_gauge_ledger()
    assert ledger.sampled_constraint_rank == 3
    assert ledger.allowed_coefficient_dimension == 1
    assert ledger.normalized_coefficient_ray == (1, -2, 2, -1)
    assert ledger.symbolic_gauge_residual == sp.zeros(10, 4)


def test_fierz_pauli_coefficient_mutation_breaks_the_gauge_kernel() -> None:
    momentum = (1, 2, 3, 5)
    accepted = fierz_pauli_gauge_residual((1, -2, 2, -1), momentum)
    mutated = fierz_pauli_gauge_residual(
        (1, -2, 2, sp.Rational(-9, 10)),
        momentum,
    )
    assert accepted == sp.zeros(10, 4)
    assert mutated.rank() > 0


@pytest.mark.parametrize("wavevector", [(0, 0, 3), (1, 2, 3), (-2, 5, 1)])
def test_volume_and_force_balance_derive_two_positive_modes(
    wavevector: tuple[int, int, int],
) -> None:
    ledger = spin_two_mode_ledger(wavevector, 2, sp.Rational(2, 3))
    assert ledger.constraint_rank == 4
    assert ledger.physical_mode_count == 2
    assert ledger.admissible_basis.shape == (6, 2)
    assert ledger.projected_kinetic_rank == 2
    assert ledger.projected_frobenius_metric.is_positive_definite is True


def test_dropping_volume_balance_exposes_the_extra_scalar_mode() -> None:
    full = collective_tensor_constraint_matrix((1, 2, 3))
    without_volume = collective_tensor_constraint_matrix(
        (1, 2, 3),
        include_volume_constraint=False,
    )
    assert tensor_mode_count_from_constraints(full) == 2
    assert tensor_mode_count_from_constraints(without_volume) == 3


def test_local_lattice_symbol_has_the_declared_relativistic_limit_and_correction() -> None:
    spacing = sp.Symbol("ell", positive=True)
    ledger = spin_two_mode_ledger((3, 0, 0), 2, spacing)
    assert ledger.continuum_angular_frequency_squared == 36
    assert sp.limit(ledger.lattice_angular_frequency_squared, spacing, 0) == 36
    assert sp.limit(
        (ledger.lattice_angular_frequency_squared - ledger.continuum_angular_frequency_squared)
        / spacing**2,
        spacing,
        0,
    ) == ledger.leading_lattice_correction / spacing**2
    assert ledger.relative_lattice_correction == -sp.Rational(3, 4) * spacing**2


def test_spin_two_composer_uses_the_wall_normalization_on_the_gauge_ray() -> None:
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    ledger = sine_gordon_spin_two_ledger(coefficients, (1, 2, 3))
    assert ledger.coupling.normalized_fierz_pauli_coefficients == (
        -sp.Rational(288, 5),
        sp.Rational(576, 5),
        -sp.Rational(576, 5),
        sp.Rational(288, 5),
    )
    assert ledger.coupling.trace_reversed_source_coefficient == -sp.Rational(5, 1152)
    assert ledger.coupling.newton_constant == sp.Rational(5, 1152) / sp.pi


def test_common_sine_gordon_scale_changes_the_derived_coupling() -> None:
    baseline = sine_gordon_spin_two_coupling_ledger(
        dimensional_sine_gordon_coefficients(2, 8, 18)
    )
    scaled = sine_gordon_spin_two_coupling_ledger(
        dimensional_sine_gordon_coefficients(14, 56, 126)
    )
    assert scaled.wall_scales.signal_speed == baseline.wall_scales.signal_speed
    assert scaled.wall_scales.profile_length == baseline.wall_scales.profile_length
    assert (
        scaled.wall_scales.spin2_spacetime_normalization
        == 7 * baseline.wall_scales.spin2_spacetime_normalization
    )
    assert scaled.newton_constant == baseline.newton_constant / 7


def test_mode_ledger_rejects_zero_or_inexact_wavevectors() -> None:
    with pytest.raises(ValueError, match="nonzero"):
        spin_two_mode_ledger((0, 0, 0), 1, 1)
    with pytest.raises(ValueError, match="exact"):
        spin_two_mode_ledger((1.0, 0, 0), 1, 1)
