"""Direct exact tests of the actual controller and retained physical action."""

import pytest
import sympy as sp

from substrate_framework import euler_fourier as ef
from substrate_framework.euler_displacement_preparation import finite_displacement_cell
from substrate_framework.euler_passive_control import (
    acoustic_action_jet,
    diagonal_wave_number,
    passive_packet,
    phase_matched_controls,
    wrapped_streamline,
)


def test_wrapped_streamline_actual_period_coarea_and_regular_domain():
    branch = wrapped_streamline(sp.Rational(1, 100), 1, sp.Rational(1, 2))
    z = branch.coordinate
    assert branch.speed_squared == 1-(sp.Rational(1, 2)-sp.cos(z)/100)**2
    assert sp.simplify(branch.frequency*branch.period-2*sp.pi) == 0
    assert sp.simplify(branch.coarea_density-branch.period/(4*sp.pi**2)) == 0
    assert branch.domain_condition is sp.true
    assert isinstance(branch.period, sp.Integral)


def test_period_derivative_is_the_actual_saddle_frequency_source():
    aa, bb, cc = sp.symbols("A B c", positive=True)
    branch = wrapped_streamline(aa, bb, cc)
    z = branch.coordinate
    expected = (cc-aa*sp.cos(z))/branch.speed_squared**sp.Rational(3, 2)
    assert sp.simplify(sp.diff(branch.period.function, cc)-expected) == 0
    assert branch.domain_condition.has(sp.Gt(bb-aa-cc, 0))


def test_passive_sector_on_the_same_canonical_euler_cell_has_zero_pressure():
    cell = finite_displacement_cell()
    scalar = ef.mul(ef.trig(1, 3), ef.trig(2, 2, kind="sin"))
    perturbation = (scalar, {}, {})
    assert not any(ef.transport(perturbation, cell.background))
    transport = ef.transport(cell.background, perturbation)
    assert not transport[1] and not transport[2]
    assert not ef.divergence(transport)
    assert ef.leray(transport) == transport


def test_actual_packet_euler_and_lin_histories_and_initial_configuration():
    theta, time = sp.symbols("theta t", real=True)
    omega, amplitude = sp.symbols("omega G", positive=True)
    packet = passive_packet(amplitude, omega, theta, time)
    assert sp.simplify(sp.diff(packet.euler_velocity, time)
                       +omega*sp.diff(packet.euler_velocity, theta)) == 0
    assert sp.simplify(packet.coordinate_rate+omega*sp.diff(packet.lin_displacement, theta)
                       -packet.euler_velocity) == 0
    assert packet.lin_displacement.subs(time, 0) == packet.initial_configuration
    assert packet.coordinate_rate.subs(time, 0) == packet.initial_velocity/2


def test_half_inverse_retains_full_energy_and_exposes_wrong_fraction():
    theta, time = sp.symbols("theta t", real=True)
    packets = [passive_packet(3, 2, theta, time, configuration_fraction=alpha)
               for alpha in (0, sp.Rational(1, 2), 1)]
    assert [p.initial_energy_angle_average for p in packets] == [sp.Rational(9, 4), 0,
                                                                              -sp.Rational(9, 4)]
    assert sp.simplify(packets[1].initial_energy_density) == 0
    assert packets[1].initial_configuration != 0


def test_exact_phase_and_output_moment_inversion_includes_the_extra_row():
    epsilon = sp.Symbol("epsilon", positive=True)
    matrix = sp.Matrix([[1/epsilon, 1/(2*epsilon)], [epsilon, 2*epsilon]])
    control = phase_matched_controls(matrix, [1])
    assert control.exponents == (-1, 1)
    assert control.weights == sp.Matrix([-1/(3*epsilon), 2/(3*epsilon)])
    assert control.phase_residual == 0 and control.output_residuals == (0,)
    assert control.domain_condition is sp.true


def test_positive_width_band_moments_are_not_replaced_by_point_carriers():
    # Exact uniform-band algebra fixture; the theorem's actual Euler bands
    # are smooth. This tests the moment API, not an Euler existence inference.
    omega = sp.Symbol("omega", positive=True)
    intervals = [(1, 2), (3, 4)]
    matrix = sp.Matrix([[sp.integrate(omega**power, (omega, lo, hi))
                         for lo, hi in intervals] for power in (-1, 1)])
    control = phase_matched_controls(matrix, [1])
    assert control.phase_residual == 0 and control.output_residuals == (0,)
    assert sp.simplify((matrix*control.weights)[1]-1) == 0
    assert matrix[0, 0] == sp.log(2)


def test_higher_odd_polynomial_coefficients_have_correct_factorials_and_signs():
    frequencies = (1, 2, 3)
    matrix = sp.Matrix([[sp.Rational(1, q) if power == -1 else q**power
                         for q in frequencies] for power in (-1, 1, 3)])
    control = phase_matched_controls(matrix, [1, 2])
    time = sp.Symbol("t", real=True)
    response = sum(w*sp.sin(q*time) for w, q in zip(control.weights, frequencies))
    assert sp.series(response, time, 0, 5).removeO() == time+2*time**3
    assert sum(w/q for w, q in zip(control.weights, frequencies)) == 0


@pytest.fixture(scope="module")
def acoustic():
    time = sp.Symbol("t", real=True)
    return time, acoustic_action_jet(2, 3, time, time**3/6)


def test_physical_acoustic_action_keeps_the_actual_moving_mass(acoustic):
    time, result = acoustic
    q = result.spatial_square
    assert sp.simplify(result.wronskian-1-q*time**2/2) == 0
    assert sp.simplify(result.mass-3*(1-q*time**2/2)) == 0
    assert sp.simplify(result.mass_rate+3*q*time) == 0
    assert result.stiffness == 6*q


def test_conserved_energy_is_distinct_from_moving_mechanical_energy(acoustic):
    time, result = acoustic
    q = result.spatial_square
    assert result.conserved_energy == sp.diag(6*q, 3*(1-q*time**2))
    assert sp.simplify(result.conserved_energy[1, 1]-result.mass) != 0
    exact = acoustic_action_jet(2, 3, time, 0)
    assert exact.mass == 3 and exact.mass_rate == 0
    assert exact.conserved_energy == sp.diag(6*exact.spatial_square, 3)


def test_actual_diagonal_scale_suppresses_arbitrarily_large_finite_norms():
    accuracy = sp.Rational(1, 1000)
    norm, constant = sp.Integer(10)**20, sp.Integer(10)**40
    chosen = diagonal_wave_number(accuracy, norm, constant, previous=sp.Rational(1, 100))
    assert 0 < chosen < sp.Rational(1, 200)
    assert chosen*norm < accuracy and chosen*constant < accuracy


def test_explicit_invalid_domains_are_rejected_without_numeric_tolerances():
    theta, time = sp.symbols("theta t", real=True)
    calls = [lambda: wrapped_streamline(1, 1, 0),
             lambda: wrapped_streamline(1, 3, 2),
             lambda: passive_packet(1, 0, theta, time),
             lambda: passive_packet(1, 2, theta, theta),
             lambda: passive_packet(sp.cos(theta), 2, theta, time),
             lambda: phase_matched_controls([[1, 1], [1, 1]], [1]),
             lambda: phase_matched_controls([[1, 2], [3, 4]], []),
             lambda: phase_matched_controls([[1]], [1]),
             lambda: phase_matched_controls([[1, -2], [3, 4]], [1]),
             lambda: acoustic_action_jet(2, 3, time, 1),
             lambda: acoustic_action_jet(2, 3, time, time),
             lambda: diagonal_wave_number(sp.nan, 1, 1),
             lambda: diagonal_wave_number(1, sp.oo, 1)]
    for call in calls:
        with pytest.raises(ValueError):
            call()
