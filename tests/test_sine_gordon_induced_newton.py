import pytest
import sympy as sp

from substrate_framework.dimensional_sine_gordon import (
    DimensionalSineGordonCoefficients,
    dimensional_sine_gordon_coefficients,
    dimensional_sine_gordon_scales,
)
from substrate_framework.induced_gravity import induced_inverse_newton_ledger
from substrate_framework.sine_gordon_induced_newton import (
    CONFORMAL_COUPLING,
    induced_newton_coefficient,
    seeley_dewitt_curvature_coefficient,
    sine_gordon_induced_newton,
)


def test_seeley_dewitt_curvature_coefficient_is_one_sixth_minus_xi() -> None:
    minimal = seeley_dewitt_curvature_coefficient(0)
    assert minimal.curvature_coefficient == sp.Rational(1, 6)
    assert minimal.vanishes_at_conformal_coupling is False
    conformal = seeley_dewitt_curvature_coefficient(CONFORMAL_COUPLING)
    assert conformal.curvature_coefficient == 0
    assert conformal.vanishes_at_conformal_coupling is True
    # the sine-Gordon vacuum mass enters a_2 only through -m**2, never the R coeff
    massive = seeley_dewitt_curvature_coefficient(sp.Rational(1, 12), sp.Symbol("m2", positive=True))
    assert massive.curvature_coefficient == sp.Rational(1, 12)
    assert massive.mass_coefficient == -sp.Symbol("m2", positive=True)


def test_proper_time_coefficient_is_exact_and_linear_in_field_count() -> None:
    one = induced_newton_coefficient(1, 0)
    assert one.coefficient == 1 / (12 * sp.pi)
    assert one.coefficient_per_field == 1 / (12 * sp.pi)
    # exact (1-6 xi)/(12 pi) at a generic rational coupling
    tilted = induced_newton_coefficient(1, sp.Rational(1, 12))
    assert sp.simplify(tilted.coefficient - (1 - 6 * sp.Rational(1, 12)) / (12 * sp.pi)) == 0
    # linear in N
    six = induced_newton_coefficient(6, 0)
    assert sp.simplify(six.coefficient - 6 * one.coefficient) == 0


def test_sign_gate_flips_at_conformal_coupling() -> None:
    assert induced_newton_coefficient(1, 0).sign == 1  # minimal: attractive
    assert induced_newton_coefficient(1, 0).attractive is True
    assert induced_newton_coefficient(1, CONFORMAL_COUPLING).sign == 0  # induced EH vanishes
    assert induced_newton_coefficient(1, CONFORMAL_COUPLING).coefficient == 0
    wrong = induced_newton_coefficient(1, sp.Rational(1, 2))  # xi > 1/6
    assert wrong.sign == -1
    assert wrong.attractive is False


def test_unknown_regulator_is_a_declared_premise_not_a_default() -> None:
    with pytest.raises(ValueError):
        induced_newton_coefficient(1, 0, regulator="momentum_cutoff")


def test_symbolic_coupling_without_decidable_sign_raises() -> None:
    with pytest.raises(ValueError):
        induced_newton_coefficient(1, sp.Symbol("xi", real=True))


def test_field_count_must_be_a_positive_integer() -> None:
    for bad in (0, -1, sp.Rational(3, 2)):
        with pytest.raises(ValueError):
            induced_newton_coefficient(bad, 0)


def test_composition_gives_induced_coupling_in_upstream_coefficients() -> None:
    lam, tension, mu = sp.symbols("lambda_ T mu", positive=True)
    coefficients = DimensionalSineGordonCoefficients(lam, tension, mu)
    result = sine_gordon_induced_newton(coefficients, 1, 0)
    s = result.coefficient.coefficient
    # Delta(1/G) = s*hbar/(ell^2 c^3) = s*mu*lambda^2/T^2, exactly
    assert sp.simplify(result.induced_inverse_newton - s * mu * lam**2 / tension**2) == 0
    assert sp.simplify(
        result.induced_inverse_newton - result.induced_inverse_newton_in_coefficients
    ) == 0
    # cross-check against the accepted induced-gravity ledger
    scales = dimensional_sine_gordon_scales(coefficients)
    ledger = induced_inverse_newton_ledger(
        cutoff_length=scales.length,
        induced_coefficient=s,
        signal_speed=scales.signal_speed,
        action_scale=scales.action,
    )
    assert sp.simplify(ledger.induced_inverse_newton - result.induced_inverse_newton) == 0


def test_induced_coupling_is_target_blind_no_newton_constant_installed() -> None:
    lam, tension, mu = sp.symbols("lambda_ T mu", positive=True)
    coefficients = DimensionalSineGordonCoefficients(lam, tension, mu)
    result = sine_gordon_induced_newton(coefficients, 1, 0)
    # the induced coupling is built only from the model coefficients (and pi),
    # never from an installed Newton constant or observed scale
    assert result.induced_inverse_newton.free_symbols <= {lam, tension, mu}
    assert not any(str(symbol) == "G" for symbol in result.induced_inverse_newton.free_symbols)


def test_induced_planck_length_is_ell_over_sqrt_s_and_none_off_sign() -> None:
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    scales = dimensional_sine_gordon_scales(coefficients)
    attractive = sine_gordon_induced_newton(coefficients, 1, 0)
    expected = sp.simplify(scales.length / sp.sqrt(attractive.coefficient.coefficient))
    assert sp.simplify(attractive.induced_planck_length - expected) == 0
    # conformal: induced EH vanishes -> no real induced Planck length
    conformal = sine_gordon_induced_newton(coefficients, 1, CONFORMAL_COUPLING)
    assert conformal.induced_planck_length is None
    assert conformal.induced_inverse_newton == 0
    # wrong sign -> no real induced Planck length
    wrong = sine_gordon_induced_newton(coefficients, 1, sp.Rational(1, 2))
    assert wrong.induced_planck_length is None
