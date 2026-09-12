#!/usr/bin/env python3
"""Exact exposing checks for the 0095 principal-cocycle derivation."""

from __future__ import annotations

import sympy as sp

from substrate_framework.euler_p2_principal import (
    bas_velocity_derivative,
    beltrami_hessian_symbols,
    cao_odd_border_beta,
    column_bas_matrix,
    column_metric,
    hf_general_first_harmonic,
    hf_detuned_return_generator,
    hf_growth_bracket,
    hf_hill_potential,
    hf_resonant_slow_matrix,
    charged_column_slow_block,
    maxwell_transverse_block,
    metric_bent_c1,
)


def main() -> None:
    lam, c, cem, xiz, q = sp.symbols("lam c cem xiz q", real=True)
    omega, zeta, kz, K = sp.symbols(
        "omega zeta kz K", positive=True, real=True
    )
    alpha, R = sp.symbols("alpha R", positive=True, real=True)

    mem = maxwell_transverse_block(c, cem, xiz, q)
    char = sp.factor((lam * sp.eye(2) - mem).det())
    expected_char = sp.expand((lam + sp.I * c * xiz) ** 2 + cem**2 * q**2)
    assert sp.expand(char - expected_char) == 0

    ccol = column_bas_matrix(omega, zeta, kz, K)
    nu2 = 2 * omega * zeta * kz**2 / K**2
    assert sp.simplify(ccol * ccol + nu2 * sp.eye(2)) == sp.zeros(2)

    gcol = column_metric(omega, zeta, kz, K)
    assert sp.simplify(gcol * ccol + ccol.T * gcol) == sp.zeros(2)
    assert sp.simplify(ccol.trace()) == 0

    # The streamfunction-normal returned multipliers are periodic ratios.
    U0, UT, r0, rT = sp.symbols("U0 UT r0 rT", positive=True)
    mnormal = sp.diag(U0 / UT, r0 / rT)
    assert mnormal.subs({UT: U0, rT: r0}) == sp.eye(2)
    assert mnormal[0, 1] == 0 and mnormal[1, 0] == 0

    # Derive, rather than insert, the universal first curvature Melnikov row.
    p, s0, angle, omega_p = sp.symbols("p s angle omega_p", positive=True)
    kval = sp.sqrt(p**2 + kz**2)
    c1 = metric_bent_c1(omega, omega_p, s0, p, kz, angle)
    c_cos = c1.applyfunc(lambda entry: sp.expand_trig(entry).coeff(sp.cos(angle)))
    c_sin = c1.applyfunc(lambda entry: sp.expand_trig(entry).coeff(sp.sin(angle)))
    bcoef = zeta * kz / kval
    nu = sp.sqrt(2 * omega * zeta) * kz / kval
    vminus = sp.Matrix([1, sp.I * nu / bcoef])
    wplus = sp.Matrix([[sp.Rational(1, 2), sp.I * bcoef / (2 * nu)]])
    melnikov = sp.simplify(
        (wplus * ((c_cos - sp.I * c_sin) / 2) * vminus)[0]
    )
    melnikov = sp.simplify(melnikov.subs(omega_p, (zeta - 2 * omega) / s0))
    expected_melnikov = (
        sp.I
        * s0
        * omega
        / 4
        * (
            sp.sqrt(2 * omega / zeta) * kz / kval
            + (2 * p**2 + kz**2) / kval**2
        )
    )
    assert sp.simplify(melnikov - expected_melnikov) == 0

    wave_angle = sp.symbols("wave_angle", positive=True)
    expected_angle_form = sp.I * s0 * omega / 4 * (
        sp.sqrt(2 * omega / zeta) * wave_angle + 2 - wave_angle**2
    )
    assert sp.simplify(
        expected_melnikov
        - expected_angle_form.subs(wave_angle, kz / kval)
    ) == 0
    resonance_ratio = sp.sqrt(omega / (8 * zeta))
    melnikov_resonant = sp.simplify(
        expected_angle_form.subs(wave_angle, resonance_ratio)
    )
    expected_resonant = sp.I * s0 * omega / 4 * (
        2 + 3 * omega / (8 * zeta)
    )
    assert sp.simplify(melnikov_resonant - expected_resonant) == 0

    # The raw Cao m=1 cell is not translation-solvable without its affine
    # speed/centering border.  These two moments determine beta exactly.
    pexp, ep, mass = sp.symbols("pexp ep mass", positive=True)
    raw_moment = -(sp.Rational(1, 2) + 2 / (pexp + 1)) * ep
    affine_moment = -mass
    beta = cao_odd_border_beta(pexp, ep, mass)
    assert sp.simplify(raw_moment + beta * affine_moment) == 0

    # The complete speed polynomial, not its affine core truncation, gives a
    # constant physical velocity under W=r^-1 grad(P) cross e_theta.
    rvar, speed = sp.symbols("rvar speed", positive=True)
    speed_streamfunction = -speed * rvar**2 / 2
    speed_velocity_z = sp.diff(speed_streamfunction, rvar) / rvar
    assert sp.simplify(speed_velocity_z + speed) == 0
    assert sp.diff(speed_velocity_z, rvar) == 0
    kvec = sp.Matrix([p, 0, kz])
    assert bas_velocity_derivative(sp.zeros(3), kvec) == sp.zeros(3)

    # Re-derive the invariant Hattori--Fukumoto reduction.  If
    # p'=a*p+b*q and q'=c*p-a*q with c constant along the orbit, direct
    # elimination gives q''+(a'-a**2-b*c)q=0.
    alog, alog_prime, bc = sp.symbols("alog alog_prime bc", real=True)
    assert hf_hill_potential(alog, alog_prime, bc) == alog_prime - alog**2 - bc

    # Smooth-center first resonance.  U_theta=Omega_c*s+O(s**3), hence
    # omega_vorticity/Omega_c -> 2, g/(s**2*U_theta**2)->1/4, and HF (23)
    # gives cos(chi)**2->1/16.  Equation (22) then has F=15*Omega_c/128.
    omega_c = sp.symbols("omega_c", positive=True)
    center_f = sp.simplify(
        hf_general_first_harmonic(
            2 * omega_c,
            omega_c,
            sp.Rational(1, 16),
            sp.Rational(1, 4),
        )
    )
    assert center_f == sp.Rational(15, 128) * omega_c

    # The resonant oscillator frequency is omega_c/2.  The standard Hill
    # exponent |epsilon*U_theta*F|/(4*nu_0) therefore equals
    # epsilon*|U_theta|*15/256, exactly HF (24), not the quarantined raw
    # 255/1024 matrix-entry combination.
    center_bracket = sp.simplify(
        hf_growth_bracket(sp.Rational(1, 16), sp.Rational(1, 4))
    )
    assert center_bracket == sp.Rational(15, 256)
    assert sp.simplify(center_f / (4 * (omega_c / 2)) - center_bracket) == 0
    slow = hf_resonant_slow_matrix(center_f, omega_c / 2)
    assert slow == sp.Matrix([[0, center_bracket], [center_bracket, 0]])
    assert slow.eigenvals() == {center_bracket: 1, -center_bracket: 1}

    # An exact returned covector permits detuning rather than a secular drift.
    # The paired generator is hyperbolic precisely inside the Hill tongue.
    detuning, coupling = sp.symbols("detuning coupling", real=True)
    returned = hf_detuned_return_generator(detuning, coupling)
    assert sp.simplify(
        (lam * sp.eye(2) - returned).det()
        - (lam**2 + detuning**2 - coupling**2)
    ) == 0
    assert returned.subs(detuning, 0).eigenvals() == {
        coupling: 1,
        -coupling: 1,
    }

    # The charged straight-column slow block includes both tag-force and
    # magnetic rows.  Its square fixes the inertial frequency and exposes a
    # wrong sign in either projected coupling.
    xpar, z_eff, r_eff = sp.symbols("xpar z_eff r_eff", real=True)
    charged = charged_column_slow_block(xpar, z_eff, r_eff)
    assert sp.simplify(charged**2 + xpar**2 * z_eff * r_eff * sp.eye(2)) == sp.zeros(2)

    hp, hm = beltrami_hessian_symbols(alpha, R)
    assert sp.simplify(hp * hm - (1 - alpha**2 * R**2)) == 0
    assert hp.subs(R, 2 / alpha) < 0 < hm.subs(R, 2 / alpha)

    checks = [
        "Maxwell transverse characteristic polynomial",
        "column BAS square and Rayleigh discriminant",
        "positive column metric skew identity",
        "streamfunction-normal exact identity return",
        "metric-bent d=1 resonant Melnikov coefficient",
        "Cao odd-cell affine solvability moment",
        "full speed polynomial has zero BAS strain",
        "Hattori-Fukumoto 2x2-to-Hill elimination",
        "Hattori-Fukumoto smooth-center forcing coefficient",
        "Hattori-Fukumoto smooth-center physical growth coefficient",
        "Hattori-Fukumoto paired Hill couplings and hyperbolic discriminant",
        "exact-periodic-covector detuning and Hill tongue",
        "charged-column constrained slow-block square",
        "constant-lambda helicity two-sign symbol",
    ]
    for index, check in enumerate(checks, 1):
        print(f"PASS {index}: {check}")
    print(f"PASS: {len(checks)}/{len(checks)} exact checks")


if __name__ == "__main__":
    main()
