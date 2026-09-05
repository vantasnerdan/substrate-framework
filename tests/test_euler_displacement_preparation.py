"""Direct exact-field tests of the displacement-only construction."""

import pytest
import sympy as sp

from substrate_framework import euler_fourier as ef
from substrate_framework.euler_displacement_preparation import (
    finite_displacement_cell,
    negative_helicity_shell,
    prepared_displacement,
    transverse_pair_average,
)


@pytest.fixture(scope="module")
def cell():
    return finite_displacement_cell()


@pytest.fixture(scope="module")
def rows(cell):
    kap = sp.Matrix(sp.symbols("kx ky kz", real=True))
    disp = sp.Matrix(sp.symbols("Dx Dy Dz", real=True))
    amplitude = sp.Symbol("amplitude", real=True)
    return kap, disp, amplitude, prepared_displacement(cell, kap, disp, amplitude=amplitude)


def test_actual_background_complete_euler_and_helicity(cell):
    u = cell.background
    for i, component in enumerate(ef.transport(u, u)):
        assert not ef.add(component, ef.derivative(cell.pressure, i))
        assert not ef.add(ef.curl(u)[i], u[i])
    assert not ef.divergence(u)


def test_full_forced_corrector_and_untruncated_kernel_inverse(cell):
    z, u = cell.velocity_corrector, cell.background
    left, right = ef.transport(u, z), ef.transport(z, u)
    generator = ef.leray(tuple(ef.add(left[i], right[i]) for i in range(3)))
    assert all(not ef.add(generator[i], ef.scale(cell.forcing[i], -1)) for i in range(3))
    assert not ef.divergence(z)
    assert all(sum(q*q for q in wave) != 1 for wave in cell.scalar_corrector)
    assert max(sum(abs(q) for q in wave) for wave in cell.scalar_corrector) == 7
    advected = ef.transport(u, (cell.first_integral, {}, {}))[0]
    assert not advected


def test_actual_range_energy_positive_matching_margin(cell):
    assert cell.range_energy == sp.Rational(
        636165619333494275181329617802709521563769,
        112592022217345440000000000000000000000000000)
    assert 0 < cell.range_energy < 13*cell.background_energy/1280
    assert cell.background_energy/240 < cell.restoring_coefficient < cell.background_energy/120


def test_return_is_stationary_full_helicity_and_has_zero_phase_mean(rows):
    returned = rows[3].returned
    assert not ef.divergence(returned)
    assert all(not ef.add(ef.curl(returned)[i], returned[i]) for i in range(3))
    assert all(component.get(ef.ZERO, 0) == 0 for component in returned)
    projected = negative_helicity_shell(returned)
    assert all(not ef.add(projected[i], ef.scale(returned[i], -1)) for i in range(3))


def test_independent_fourth_tensor_matches_general_pair_integral(rows):
    kap, disp = rows[:2]
    for i, j, a, b in ((0, 0, 0, 0), (0, 0, 1, 1), (0, 1, 0, 1), (2, 1, 0, 0)):
        expected = sp.Rational(2, 15)*int(i == j and a == b)
        expected -= sp.Rational(1, 30)*(int(i == a and j == b)+int(i == b and j == a))
        assert transverse_pair_average(kap[i]*kap[j]*disp[a]*disp[b], kap, disp) == expected
    assert transverse_pair_average(disp.dot(disp), kap, disp) == 1
    assert transverse_pair_average(kap.dot(disp)**2, kap, disp) == 0
    assert transverse_pair_average(0, kap, disp) == 0


def test_complete_field_energy_and_actual_current_match_independent_constants(cell, rows):
    kap, disp, amplitude, prepared = rows
    energy = transverse_pair_average(prepared.energy_coefficient, kap, disp)
    acceleration = transverse_pair_average(prepared.acceleration_contraction, kap, disp)
    assert sp.factor(energy-cell.background_energy*(amplitude**2/15-sp.Rational(47, 240))
                     -cell.range_energy) == 0
    assert sp.factor(acceleration-cell.background_energy*(8*amplitude+13)/120) == 0
    assert sp.simplify(energy.subs(amplitude, cell.matching_amplitude)
                       -cell.restoring_coefficient) == 0
    assert sp.simplify(acceleration.subs(amplitude, cell.matching_amplitude)
                       +cell.restoring_coefficient) == 0


def test_exact_finite_k_lift_has_actual_bloch_divergence(cell):
    kap, disp = sp.Matrix([0, 1, 0]), sp.Matrix([0, 0, 1])
    prepared = prepared_displacement(cell, kap, disp, amplitude=0)
    wave_number = sp.Symbol("k", real=True)
    for wave in set().union(*(component.keys() for component in prepared.translation)):
        full_wave = sp.Matrix(wave)+wave_number*kap
        full_value = sp.Matrix([prepared.translation[i].get(wave, 0)
                               +sp.I*wave_number*prepared.lift[i].get(wave, 0) for i in range(3)])
        assert sp.expand(full_wave.dot(full_value)) == 0


def test_omitted_actual_euler_subtraction_is_exposed(rows):
    kap, disp, _, prepared = rows
    omitted = ef.inner(prepared.material_rate, prepared.material_rate)-prepared.energy_coefficient
    assert transverse_pair_average(omitted, kap, disp) > 0


def test_pair_average_rejects_wrong_degree_and_repeated_coordinates(rows):
    kap, disp = rows[:2]
    with pytest.raises(ValueError, match="quadratic"):
        transverse_pair_average(disp[0], kap, disp)
    with pytest.raises(ValueError, match="disjoint"):
        transverse_pair_average(disp[0]**2, disp, disp)
    with pytest.raises(ValueError, match="disjoint"):
        transverse_pair_average(disp[0]**2, [0, 1, 2], disp)


def test_preparation_input_shapes_and_real_amplitude(cell):
    with pytest.raises(ValueError, match="three"):
        prepared_displacement(cell, [1, 0], [0, 1, 0])
    with pytest.raises(ValueError, match="finite and real"):
        prepared_displacement(cell, [1, 0, 0], [0, 1, 0], amplitude=sp.I)
    with pytest.raises(ValueError, match="finite and real"):
        prepared_displacement(cell, [1, 0, 0], [0, 1, 0], amplitude=sp.nan)
    with pytest.raises(ValueError, match="three"):
        negative_helicity_shell(({}, {}))
