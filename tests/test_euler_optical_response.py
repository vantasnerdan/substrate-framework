"""Independent physical-row/energy differentiation for the conditional optical API."""

import pytest
import sympy as sp

from substrate_framework.euler_optical_response import (
    OpticalModeJet,
    correlated_optical_preparation,
    positive_optical_curvature,
    transverse_optical_action,
    transverse_tilt_tensor,
)
from substrate_framework.micropolar import micropolar_fourier_stiffness


def example():
    return correlated_optical_preparation(
        (OpticalModeJet(1, 2, 0, 1), OpticalModeJet(1, 1, 0, 1)),
        (sp.Rational(1, 2), sp.Rational(1, 2)), 1,
    )


def test_linear_system_and_signed_probability_distinction():
    result = example()
    assert result.amplitudes == (-sp.Rational(2, 3), sp.Rational(8, 3))
    assert result.slopes == (sp.Rational(9, 20), sp.Rational(1, 10))
    assert (result.J0, result.J2) == (sp.Rational(34, 9), sp.Rational(17, 80))
    assert result.slope_matrix*sp.Matrix(result.slopes) == result.slope_rhs
    assert result.domain_condition is sp.true
    assert all(r == 0 for r in (result.normalization_residual, result.variance_residual,
                                result.curvature_residual, result.energy_residual))
    # Taking absolute amplitudes destroys the physical variance cancellation.
    assert sum(w*abs(a)*m.frequency_first**2
               for w, a, m in zip(result.weights, result.amplitudes, result.modes)) != 0
    # Signed probabilities/phase weights would instead return a different mass.
    assert sum(w*a*m.phase_mass
               for w, a, m in zip(result.weights, result.amplitudes, result.modes)) != result.J0


@pytest.mark.parametrize("mass_jets", [(1, 0, 0, 1, 0, 0), (2, 3, -1, 5, -2, 4)])
def test_actual_observed_rows_and_energy_by_carrier_differentiation(mass_jets):
    j1, jp1, jpp1, j2, jp2, jpp2 = mass_jets
    modes = (OpticalModeJet(3, 2, -4, j1, jp1, jpp1),
             OpticalModeJet(3, -1, 5, j2, jp2, jpp2))
    result = correlated_optical_preparation(modes, (sp.Rational(1, 3), sp.Rational(2, 3)), 7)
    p, time, q0, rate0, frequency = sp.symbols("p time q0 rate0 frequency", real=True)
    observed = [sp.S.Zero, sp.S.Zero]
    energy, initial_phase = sp.S.Zero, sp.S.Zero
    for w, m, a, d in zip(result.weights, modes, result.amplitudes, result.slopes):
        nu = m.frequency+m.frequency_first*p+m.frequency_second*p**2/2
        mass = m.phase_mass+m.phase_mass_first*p+m.phase_mass_second*p**2/2
        amplitude = a+d*p
        observed[0] += w*amplitude*sp.cos(nu*time)
        observed[1] += w*amplitude*sp.sin(nu*time)/nu
        initial_phase += w*mass*amplitude**2
        energy += w*mass*amplitude**2*(nu**2*q0**2+rate0**2)/2
    reference = (sp.cos(frequency*time), sp.sin(frequency*time)/frequency)
    for row, ref in zip(observed, reference):
        assert sp.simplify(row.subs(p, 0)-ref.subs(frequency, 3)) == 0
        # Identity for every time, not just an initial time Taylor coefficient.
        expected = 7*sp.diff(ref, frequency).subs(frequency, 3)
        assert sp.simplify(sp.diff(row, p, 2).subs(p, 0)-expected) == 0
    assert sp.simplify(initial_phase.subs(p, 0)-result.J0) == 0
    assert sp.simplify(sp.diff(initial_phase, p, 2).subs(p, 0)-result.J2) == 0
    # The coefficient of a symmetric signed carrier band is E''/2.
    expected = (result.J2*(rate0**2+9*q0**2)+2*result.J0*3*7*q0**2)/4
    assert sp.simplify(sp.diff(energy, p, 2).subs(p, 0)/2-expected) == 0
    if jp1:
        omitted_mass_derivative = sum(
            2*w*3*m.phase_mass_first*a**2*m.frequency_first
            for w, m, a in zip(result.weights, modes, result.amplitudes))
        assert omitted_mass_derivative != 0


def test_full_transverse_initial_phase_and_raw_normalization():
    raw, nu, amplitude = sp.symbols("raw nu amplitude", positive=True)
    j = sp.Matrix([[0, 1], [-1, 0]])
    for sign in (-1, 1):
        phase_map = sp.diag(amplitude, sign*amplitude/nu)
        pulled = phase_map.T*(sign*raw*nu*j)*phase_map
        assert sp.simplify(pulled-raw*amplitude**2*j) == sp.zeros(2)
        energy = phase_map.T*(raw*nu**2*sp.eye(2))*phase_map
        assert sp.simplify(energy-raw*amplitude**2*sp.diag(nu**2, 1)) == sp.zeros(2)
    mode = OpticalModeJet.from_raw_tilt(2, -3, 4, 15, 6, -9)
    assert (mode.phase_mass, mode.phase_mass_first, mode.phase_mass_second) == (5, 2, -3)


def test_actual_haar_frame_integral_and_rotation_covariance():
    polar, azimuth, spin = sp.symbols("polar azimuth spin", real=True)
    # Explicit orthonormal frame; carrier t is NOT the observed tilt n.
    t_z = sp.cos(polar)
    n_x = sp.cos(polar)*sp.cos(azimuth)*sp.cos(spin)-sp.sin(azimuth)*sp.sin(spin)
    n_z = -sp.sin(polar)*sp.cos(spin)
    values = []
    for component in (n_x, n_z):
        integral = 3*component**2*t_z**2*sp.sin(polar)/(8*sp.pi**2)
        for coordinate, end in ((spin, 2*sp.pi), (azimuth, 2*sp.pi), (polar, sp.pi)):
            integral = sp.integrate(integral, (coordinate, 0, end))
        values.append(sp.simplify(integral))
    assert values == [sp.Rational(2, 5), sp.Rational(1, 5)]
    k = sp.Matrix([2, -3, 4])
    rotation = sp.Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    assert transverse_tilt_tensor(rotation*k) == rotation*transverse_tilt_tensor(k)*rotation.T
    assert transverse_tilt_tensor([0, 0, 1]) == sp.diag(*values[:1]*2, values[1])


def test_action_keeps_gradient_mass_and_truncates_products():
    result = example()
    k, omega = sp.symbols("k omega", real=True)
    action = transverse_optical_action(result, [0, 0, k])
    for axis, factor in enumerate((sp.Rational(2, 5), sp.Rational(2, 5), sp.Rational(1, 5))):
        mass = action.mass[axis, axis]
        stiffness = action.stiffness[axis, axis]
        assert mass == result.J0+factor*k**2*result.J2/2
        assert sp.series(stiffness/mass, k, 0, 3).removeO().expand() == 1+factor*k**2
        assert sp.Poly(stiffness, k).degree() == 2
        assert sp.simplify(sp.diff(stiffness-omega**2*mass, omega, 2)+2*mass) == 0
    product = action.mass*action.frequency_squared
    assert sp.simplify(product-action.stiffness) != sp.zeros(3)  # O(k^4) is not retained.
    assert all(sp.series(x, k, 0, 3).removeO() == 0 for x in product-action.stiffness)


def test_positive_representative_matches_canonical_symbol_and_trace_sector():
    ct, cl, k = sp.symbols("ct cl k", positive=True)
    c = positive_optical_curvature(ct, cl)
    stiffness = micropolar_fourier_stiffness([0, 0, k], c)[3:, 3:]
    assert sp.simplify(stiffness-sp.diag(ct*k**2, ct*k**2, cl*k**2)) == sp.zeros(3)
    assert sp.simplify(3*c.trace_curvature+c.symmetric_curvature
                       -9*cl**2/(2*(4*ct+3*cl))) == 0
    assert c.symmetric_curvature.is_positive and c.skew_curvature.is_positive
    # Negative c_tr alone is not a negative trace-sector energy.
    example_coefficients = positive_optical_curvature(10, 1)
    assert example_coefficients.trace_curvature < 0
    assert example_coefficients.trace_curvature+example_coefficients.symmetric_curvature/3 > 0


@pytest.mark.parametrize("bad", [0, -1, sp.I, sp.oo, sp.nan, sp.zoo])
def test_invalid_mode_domain(bad):
    with pytest.raises(ValueError):
        OpticalModeJet(bad, 2, 0, 1)
    with pytest.raises(ValueError):
        OpticalModeJet(1, 2, 0, bad)


@pytest.mark.parametrize("slopes", [(0, 1), (1, 1), (1, -1)])
def test_degenerate_slope_system(slopes):
    modes = tuple(OpticalModeJet(1, slope, 0, 1) for slope in slopes)
    with pytest.raises(ValueError):
        correlated_optical_preparation(modes, (sp.Rational(1, 2),)*2, 1)


def test_other_invalid_inputs_and_symbolic_hypotheses():
    modes = (OpticalModeJet(1, 2, 0, 1), OpticalModeJet(1, 1, 0, 1))
    for weights, target in (((1, 1), 1), ((0, 1), 1), ((-1, 2), 1),
                            ((sp.Rational(1, 2),)*2, 0)):
        with pytest.raises(ValueError):
            correlated_optical_preparation(modes, weights, target)
    with pytest.raises(ValueError):
        correlated_optical_preparation((modes[0], OpticalModeJet(2, 1, 0, 1)),
                                       (sp.Rational(1, 2),)*2, 1)
    with pytest.raises(ValueError):
        transverse_tilt_tensor([1, 2])
    with pytest.raises(ValueError):
        transverse_tilt_tensor([1, sp.I, 0])
    with pytest.raises(ValueError):
        positive_optical_curvature(-1, 2)
    unknown = sp.Symbol("unknown", real=True)
    result = correlated_optical_preparation(
        (OpticalModeJet(1, 2, 0, unknown), modes[1]), (sp.Rational(1, 2),)*2, 1)
    assert result.domain_condition == (unknown > 0)
