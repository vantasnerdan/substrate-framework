"""Machine checks for src/substrate_framework/chromomagnetic_sectors.py (P229, #47/#48)."""

import numpy as np
import sympy as sp
import pytest

from substrate_framework import chromomagnetic_sectors as cs
from substrate_framework.chromomagnetic_two_loop import b, g
from substrate_framework.su3 import fundamental_generators


# ---------------------------------------------------------------------------
# Sector coefficients and the transverse dominance of the log
# ---------------------------------------------------------------------------


def test_sector_coefficients_sum_to_11_over_3():
    coeffs = cs.sector_heat_kernel_coefficients()
    assert sum(coeffs.values()) == sp.Rational(11, 3)


def test_each_transverse_sector_carries_11_over_6():
    coeffs = cs.sector_heat_kernel_coefficients()
    assert coeffs["transverse_spin_up"] == sp.Rational(11, 6)
    assert coeffs["transverse_spin_down"] == sp.Rational(11, 6)


def test_longitudinal_ghost_cancellation():
    """Individually the spin-0 pair and the ghost carry ∓1/3 of a (gB)^2 s
    coefficient; they cancel exactly, leaving the log to the transverse spins."""
    coeffs = cs.sector_heat_kernel_coefficients()
    assert coeffs["longitudinal_pair"] == sp.Rational(-1, 3)
    assert coeffs["ghost"] == sp.Rational(1, 3)
    assert coeffs["longitudinal_pair"] + coeffs["ghost"] == 0


def test_transverse_share_is_one():
    assert cs.transverse_share_of_log() == sp.Integer(1)


def test_mutation_ghost_weight_breaks_cancellation():
    """Ghost weight -1 instead of -2 must change the total coefficient."""
    s, gB = sp.symbols("s gB", positive=True)
    series = sp.series(
        gB / sp.sinh(gB * s) * (2 * sp.cosh(2 * gB * s) + 2 - 1), gB, 0, 4
    )
    total = sp.Rational(sp.expand(series.removeO()).coeff(gB, 2) / s)
    assert total != sp.Rational(11, 3)


# ---------------------------------------------------------------------------
# SU(3) color decomposition
# ---------------------------------------------------------------------------


def test_su3_roots_three_tachyonic_sectors():
    assert cs.su3_charged_roots() == (sp.Integer(1), sp.Rational(1, 2), sp.Rational(1, 2))


def test_su3_potential_log_coefficient():
    """SU(3) lambda_3 background: the (gH)^2 ln(gH) coefficient is 3C/2 —
    one full-charge root plus two half-charge roots."""
    gH, mu2 = sp.symbols("gH mu2", positive=True)
    v = cs.su3_one_loop_potential(gH, mu2)
    from substrate_framework.chromomagnetic_background import one_loop_log_coefficient

    c = one_loop_log_coefficient()
    expanded = sp.expand_log(v.subs(mu2, 1), force=True).expand()
    target = expanded.coeff(sp.log(gH))
    assert sp.simplify(target / gH**2 - 3 * c / 2) == 0


# ---------------------------------------------------------------------------
# Stability
# ---------------------------------------------------------------------------


def test_im_v_two_loop_exact():
    expected = -3 * g**2 * b**2 * sp.log(2) / (256 * sp.pi**3)
    assert sp.simplify(cs.im_v_two_loop() - expected) == 0


def test_im_v_one_loop_value():
    gB = sp.symbols("gB", positive=True)
    assert sp.simplify(cs.im_v_one_loop(gB) - gB**2 / (8 * sp.pi)) == 0


def test_verdict_structure():
    v = cs.stability_verdict()
    assert v["verdict"].startswith("unstable")
    assert v["fluctuation_tachyon"].startswith("E^2")


def test_su3_root_charges_from_cartan_matrices():
    """The root charges must follow from the Cartan generator itself: with
    T3 = diag(1,-1,0)/2 the E_ij charge is t_ii - t_jj -> (1, 1/2, 1/2) in
    magnitude (lambda_3 background); T8 = diag(1,1,-2)/(2 sqrt(3)) would give
    (0, sqrt(3)/2, sqrt(3)/2) — guards the generator label against drift."""
    generators = fundamental_generators()
    t3 = generators[2]
    t8 = generators[7]
    pairs = [(0, 1), (0, 2), (1, 2)]
    q3 = sorted(abs(t3[i, i] - t3[j, j]) for i, j in pairs)
    q8 = sorted(abs(sp.simplify(t8[i, i] - t8[j, j])) for i, j in pairs)
    assert q3 == sorted(cs.su3_charged_roots())
    assert q8 != q3  # mutation guard: the lambda_8 pattern is NOT ours
