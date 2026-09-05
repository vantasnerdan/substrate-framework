"""Actual pressure/preparation and physical-current tests for C-CST-013."""

import pytest
import sympy as s

from substrate_framework import euler_fourier as ef
from substrate_framework.euler_acoustic import (
    observed_acoustic_cell_rows, prepared_acoustic_cell_rows,
)


def one_wave():
    return ef.trig(2), ef.trig(2, kind="sin"), {}


def test_full_prepared_one_wave_matches_independent_oblique_covariance():
    u = one_wave()
    kappa = s.Matrix([s.Rational(3, 5), 0, s.Rational(4, 5)])
    t = s.Symbol("t", real=True)
    for d, speed_squared in ((s.Matrix([0, 1, 0]), s.Rational(9, 50)),
                             (s.Matrix([s.Rational(4, 5), 0, -s.Rational(3, 5)]),
                              s.Rational(351, 1250))):
        v = 2*d
        rows = prepared_acoustic_cell_rows(u, {}, kappa, d, v)
        # In a one-coordinate field the full microscopic pressure removes z.
        cartesian_u = s.Matrix([s.cos(s.Symbol("z")), s.sin(s.Symbol("z")), 0])
        pn = s.diag(1, 1, 0)
        independent_initial = -pn*((kappa.dot(cartesian_u))*d+kappa*cartesian_u.dot(d))
        z = s.Symbol("z")
        actual_initial = s.Matrix([sum(value*s.exp(s.I*wave[2]*z)
                                      for wave, value in field.items()) for field in rows.initial_rate])
        assert all(s.simplify(s.expand_complex(q)) == 0 for q in
                   (actual_initial-independent_initial).subs(z, s.Symbol("zr", real=True)))
        assert not any(rows.forcing_rate)
        cell = tuple(ef.add(ef.scale(rows.initial_rate[i], t),
                            ef.scale(rows.forcing_constant[i], t**2/2)) for i in range(3))
        rate = tuple(ef.add(rows.initial_rate[i], ef.scale(rows.forcing_constant[i], t))
                     for i in range(3))
        observation = observed_acoustic_cell_rows(u, {}, kappa, d+t*v, cell, rate)
        assert observation.acceleration == -speed_squared*(d+t*v)


def test_nonconstant_pressure_and_eulerian_current_against_cartesian_integrals():
    y, z, t = s.symbols("y z t", real=True)
    u = (ef.add(ef.trig(2), ef.scale(ef.trig(1, kind="sin"), 2)),
         ef.trig(2, kind="sin"), ef.scale(ef.trig(1), 2))
    pressure = ef.scale(ef.add(*(ef.mul(v, v) for v in u)), -s.Rational(1, 2))
    dp = tuple(ef.derivative(pressure, i) for i in range(3))
    cell = (ef.add(dp[1], ef.scale(ef.trig(2, kind="sin"), t)), dp[2], ef.scale(dp[1], -1))
    rate = (ef.trig(2, kind="sin"), {}, {})
    result = observed_acoustic_cell_rows(u, pressure, [0, 1, 0], [0, 0, 0], cell, rate)

    uc = s.Matrix([s.cos(z)+2*s.sin(y), s.sin(z), 2*s.cos(y)])
    pc = -uc.dot(uc)/2
    grad = s.Matrix([0, s.diff(pc, y), s.diff(pc, z)])
    chic = s.Matrix([grad[1]+t*s.sin(z), grad[2], -grad[1]])
    ratec = chic.diff(t)

    def mean(expr):
        return s.simplify(s.integrate(s.integrate(s.expand_trig(expr), (y, 0, 2*s.pi)),
                                     (z, 0, 2*s.pi))/(4*s.pi**2))

    pressure_one = (grad[1]*chic).applyfunc(mean)
    pressure_two = (grad*chic[1]).applyfunc(mean)
    raw = (uc[1]*ratec+uc*ratec[1]).applyfunc(mean)+pressure_one+pressure_two
    correction = (-uc[1]*chic+uc*chic[1]).applyfunc(mean)
    slow = s.diag(1, 0, 1)
    assert result.acceleration == slow*raw == s.Matrix([s.Rational(3, 2), 0, 0])
    assert result.current_correction == slow*correction == s.Matrix([-t/2, 0, 0])
    # Omitting either pressure return is exposed, even where they cancel.
    assert slow*pressure_one != s.zeros(3, 1)
    assert slow*pressure_two != s.zeros(3, 1)
    assert result.acceleration != slow*(2*uc[1]*ratec).applyfunc(mean)+slow*(pressure_one+pressure_two)


def test_preparation_keeps_displacement_and_velocity_sources_distinct():
    u = one_wave()
    d = [0, 0, 1]
    rows = prepared_acoustic_cell_rows(u, {}, [1, 0, 0], d, [0, 0, 0])
    common = prepared_acoustic_cell_rows(u, {}, [1, 0, 0], [0, 0, 0], d)
    assert not any(rows.initial_rate)  # Microscopic pressure kills this column.
    assert not any(common.initial_rate)
    assert not any(rows.forcing_constant) and not any(common.forcing_constant)
    response = observed_acoustic_cell_rows(u, {}, [1, 0, 0], d, ({}, {}, {}), ({}, {}, {}))
    # Actual negative SV stiffness; this API never substitutes a positive target.
    assert response.acceleration == s.Matrix([0, 0, s.Rational(1, 2)])
    sh = prepared_acoustic_cell_rows(u, {}, [1, 0, 0], [0, 1, 0], [0, 0, 0])
    sh_velocity = prepared_acoustic_cell_rows(u, {}, [1, 0, 0], [0, 0, 0], [0, 1, 0])
    assert any(sh.initial_rate) and any(sh_velocity.forcing_constant)
    assert sh.initial_rate != sh_velocity.forcing_constant


def test_response_rejects_incompatible_domains_and_stationary_pressure():
    u = one_wave()
    with pytest.raises(ValueError, match="unit"):
        prepared_acoustic_cell_rows(u, {}, [0, 0, 2], [1, 0, 0], [0, 0, 0])
    with pytest.raises(ValueError, match="transverse"):
        prepared_acoustic_cell_rows(u, {}, [1, 0, 0], [1, 0, 0], [0, 0, 0])
    with pytest.raises(ValueError, match="stationary Euler"):
        prepared_acoustic_cell_rows(u, ef.trig(2), [1, 0, 0], [0, 1, 0], [0, 0, 0])
    with pytest.raises(ValueError, match="mean-zero"):
        observed_acoustic_cell_rows(u, {}, [1, 0, 0], [0, 1, 0],
                                    ({ef.ZERO: 1}, {}, {}), ({}, {}, {}))
    with pytest.raises(ValueError, match="real, commensurate"):
        prepared_acoustic_cell_rows(({(0, 0, 1): 1}, {}, {}), {},
                                     [1, 0, 0], [0, 1, 0], [0, 0, 0])
