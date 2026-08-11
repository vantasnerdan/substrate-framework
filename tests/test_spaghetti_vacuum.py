"""Machine checks for src/substrate_framework/spaghetti_vacuum.py (P229, issue #49)."""

import numpy as np
import sympy as sp
import pytest

from substrate_framework import spaghetti_vacuum as sv


# ---------------------------------------------------------------------------
# Abrikosov parameter from the defining condensate density
# ---------------------------------------------------------------------------


def test_condensate_density_periodicity():
    tau = 1j * 1.0
    for z in (0.3 + 0.2j, 0.7 + 0.55j):
        assert sv.condensate_density(z, tau) == pytest.approx(
            sv.condensate_density(z + 1, tau), rel=1e-10
        )
        assert sv.condensate_density(z, tau) == pytest.approx(
            sv.condensate_density(z + tau, tau), rel=1e-10
        )


def test_beta_square_matches_literature():
    assert sv.beta_square(ngrid=120) == pytest.approx(1.180340, rel=1e-5)


def test_beta_triangular_matches_literature():
    assert sv.beta_triangular(ngrid=120) == pytest.approx(1.159595, rel=1e-5)


def test_triangular_lattice_is_optimal():
    opt = sv.optimal_lattice()
    assert opt["optimal"] == "triangular"
    assert opt["beta"] == pytest.approx(1.159595, rel=1e-5)


def test_mutation_oblique_lattice_has_larger_beta():
    """A non-optimal lattice (e.g. tau = 2i) must have beta above triangular."""
    beta_oblique = sv.abrikosov_beta(2j, ngrid=120)
    assert beta_oblique > sv.beta_triangular(ngrid=120)


# ---------------------------------------------------------------------------
# Classical energy and total potential
# ---------------------------------------------------------------------------


def test_classical_energy_factor_below_unity():
    beta = 1.159595
    f = sv.classical_energy_factor(beta)
    assert 0 < f < 1
    assert f == pytest.approx(1 - 1 / beta)


def test_spaghetti_minimum_deeper_than_savvidy():
    """The lattice minimum is deeper; depth ratio exp(24 pi^2/(11 g^2 beta))."""
    g, mu2 = sp.symbols("g mu2", positive=True)
    beta = sp.Rational(1159595, 1000000)
    ln_ratio, v_min = sv.spaghetti_minimum(g, mu2, beta)
    from substrate_framework.chromomagnetic_background import (
        one_loop_log_coefficient,
        savvidy_minimum,
    )

    c = one_loop_log_coefficient()
    f = 1 - 1 / beta
    assert sp.simplify(ln_ratio + f * 24 * sp.pi**2 / (11 * g**2)) == 0
    _, v_savvidy = savvidy_minimum(g, mu2)
    # depth ratio (exponent doubles through (gH_min)^2)
    ratio = sp.simplify(v_min / v_savvidy)
    assert sp.simplify(sp.log(ratio) - 48 * sp.pi**2 / (11 * g**2 * beta)) == 0


def test_spaghetti_stationarity_numeric():
    gval, mu2val, beta = 0.5, 1.0, 1.159595
    H = sp.symbols("H", positive=True)
    v = sv.spaghetti_potential(H, gval, mu2val, beta)
    ln_ratio, _ = sv.spaghetti_minimum(gval, mu2val, beta)
    h_min = np.exp(float(ln_ratio)) / gval
    dv = sp.diff(v, H)
    assert abs(float(dv.subs(H, h_min))) < 1e-12
