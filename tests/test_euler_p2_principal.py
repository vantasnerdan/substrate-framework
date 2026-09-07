import sympy as sp

from substrate_framework.euler_p2_principal import (
    charged_column_slow_block,
    hf_detuned_return_generator,
    hf_general_first_harmonic,
    hf_growth_bracket,
    hf_hill_potential,
    hf_resonant_slow_matrix,
)


def test_hf_trace_free_elimination_and_smooth_center_coefficient():
    a, a_prime, bc = sp.symbols("a a_prime bc", real=True)
    assert hf_hill_potential(a, a_prime, bc) == a_prime - a**2 - bc

    omega = sp.symbols("omega", positive=True)
    first = hf_general_first_harmonic(
        2 * omega, omega, sp.Rational(1, 16), sp.Rational(1, 4)
    )
    bracket = hf_growth_bracket(sp.Rational(1, 16), sp.Rational(1, 4))
    assert first == sp.Rational(15, 128) * omega
    assert bracket == sp.Rational(15, 256)
    assert first / (4 * (omega / 2)) == bracket
    assert hf_resonant_slow_matrix(first, omega / 2).eigenvals() == {
        bracket: 1,
        -bracket: 1,
    }


def test_detuning_exposes_paired_hyperbolic_tongue():
    lam, detuning, coupling = sp.symbols("lam detuning coupling", real=True)
    generator = hf_detuned_return_generator(detuning, coupling)
    assert sp.simplify(
        (lam * sp.eye(2) - generator).det()
        - (lam**2 + detuning**2 - coupling**2)
    ) == 0


def test_charged_column_block_has_inertial_frequency_square():
    x, z_eff, r_eff = sp.symbols("x z_eff r_eff", real=True)
    block = charged_column_slow_block(x, z_eff, r_eff)
    assert sp.simplify(block**2 + x**2 * z_eff * r_eff * sp.eye(2)) == sp.zeros(2)
