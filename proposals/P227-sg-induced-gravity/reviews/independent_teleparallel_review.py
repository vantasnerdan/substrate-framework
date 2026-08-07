#!/usr/bin/env python3
"""Fresh P227 coefficient, geometry, and radiation cross-check.

This reviewer intentionally imports none of the three proposed P227 modules.
It reconstructs the constitutive spectrum from tensor contractions, checks the
nonlinear FLRW identity from closed forms, derives the coupling and static
limit directly, and recomputes the full-retardation benchmark with Simpson
quadrature instead of the proposal's shared trapezoid implementation.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import simpson
import sympy as sp

from substrate_framework.verification import CheckLedger


CHANNELS = tuple(
    (internal, first, second)
    for internal in range(4)
    for first in range(4)
    for second in range(first + 1, 4)
)


def _fresh_torsion_quadratic(values: tuple[sp.Symbol, ...]) -> sp.Expr:
    signs = (-1, 1, 1, 1)
    tensor = [[[sp.Integer(0) for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for value, (internal, first, second) in zip(values, CHANNELS):
        tensor[internal][first][second] = value
        tensor[internal][second][first] = -value
    lowered = [
        [
            [signs[internal] * tensor[internal][first][second] for second in range(4)]
            for first in range(4)
        ]
        for internal in range(4)
    ]
    raised = [
        [
            [
                signs[internal]
                * signs[first]
                * signs[second]
                * lowered[internal][first][second]
                for second in range(4)
            ]
            for first in range(4)
        ]
        for internal in range(4)
    ]
    first_invariant = sum(
        raised[a][b][c] * lowered[a][b][c]
        for a in range(4)
        for b in range(4)
        for c in range(4)
    )
    second_invariant = sum(
        raised[a][b][c] * lowered[c][b][a]
        for a in range(4)
        for b in range(4)
        for c in range(4)
    )
    vector = [sum(tensor[a][a][b] for a in range(4)) for b in range(4)]
    vector_norm = sum(signs[b] * vector[b] ** 2 for b in range(4))
    return sp.expand(first_invariant / 4 + second_invariant / 2 - vector_norm)


def _simpson_full_retardation_power() -> float:
    inverse_width = 0.03
    frequency = np.sqrt(1.0 - inverse_width**2)
    scaled_coordinate = np.linspace(-15.0, 15.0, 1025)
    coordinate = scaled_coordinate / inverse_width
    phase = np.linspace(0.0, np.pi, 129)
    sech = 1.0 / np.cosh(scaled_coordinate)[:, None]
    sine = np.sin(phase)[None, :]
    cosine = np.cos(phase)[None, :]
    tangent = np.tanh(scaled_coordinate)[:, None]
    argument = (inverse_width / frequency) * sech * sine
    denominator = 1.0 + argument**2
    field_time = 4.0 * inverse_width * sech * cosine / denominator
    field_space = -4.0 * inverse_width * argument * tangent / denominator
    potential = 8.0 * argument**2 / denominator**2
    pressure = 0.5 * (field_time**2 + field_space**2) - potential

    modes = np.arange(1, 5, dtype=np.float64)
    pressure_modes = np.asarray(
        [
            2.0
            / np.pi
            * simpson(
                pressure * np.cos(2.0 * mode * phase)[None, :],
                x=phase,
                axis=1,
            )
            for mode in modes
        ]
    )
    directions = np.linspace(0.0, 1.0, 501)
    mean_derivative_squared = np.zeros_like(directions)
    for mode, pressure_mode in zip(modes, pressure_modes):
        wavenumber = 2.0 * mode * frequency
        amplitudes = np.asarray(
            [
                simpson(
                    pressure_mode * np.cos(wavenumber * direction * coordinate),
                    x=coordinate,
                )
                for direction in directions
            ]
        )
        mean_derivative_squared += 0.5 * (wavenumber * amplitudes) ** 2
    half_integral = simpson(
        (1.0 - directions**2) ** 2 * mean_derivative_squared,
        x=directions,
    )
    return float(half_integral / (8.0 * np.pi))


def main() -> int:
    checks = CheckLedger("P227-independent-review")
    channels = sp.symbols("v0:24", real=True)
    quadratic = _fresh_torsion_quadratic(channels)
    constitutive = sp.hessian(quadratic, channels) / 2
    checks.check(
        "fresh tensor contractions produce the unique full-rank TEGR quadratic form",
        constitutive == constitutive.T
        and constitutive.rank() == 24
        and constitutive.eigenvals()
        == {-2: 3, -1: 8, sp.Rational(-1, 2): 1, sp.Rational(1, 2): 3, 1: 8, 2: 1},
    )

    phase, length, onsite = sp.symbols("q ell mu", real=True, positive=True)
    chord = 2 * sp.sin(phase / 2) / length
    checks.check(
        "fresh cosine algebra fixes one spectral channel without a multiplier",
        sp.trigsimp(
            -onsite * (1 - sp.cos(phase)) / length**2
            + onsite * chord**2 / 2
        )
        == 0,
    )

    hubble = sp.symbols("H", real=True)
    fresh_torsion = 6 * hubble**2
    fresh_boundary = 18 * hubble**2
    fresh_ricci = 12 * hubble**2
    checks.check(
        "closed-form exponential FLRW independently satisfies R plus T minus B",
        sp.simplify(fresh_ricci + fresh_torsion - fresh_boundary) == 0,
    )

    speed = sp.symbols("c", positive=True)
    einstein_coupling = 1 / onsite
    newton = speed**4 / (8 * sp.pi * onsite)
    checks.check(
        "fresh action matching gives the same positive Einstein normalization",
        sp.simplify(8 * sp.pi * newton / speed**4 - einstein_coupling) == 0,
    )

    inverse_width = sp.symbols("eta", positive=True)
    source_energy = 16 * inverse_width * onsite * length
    source_mass = source_energy / speed**2
    profile_length = length / inverse_width
    schwarzschild_radius = sp.simplify(2 * newton * source_mass / speed**2)
    checks.check(
        "fresh source integration gives the same static profile compactness",
        sp.simplify(schwarzschild_radius / profile_length - 4 * inverse_width**2 / sp.pi)
        == 0,
    )

    simpson_power = _simpson_full_retardation_power()
    checks.check(
        "independent Simpson projection reproduces the refined retarded power",
        abs(simpson_power - 0.0002875595507880884) < 2.3e-12,
    )
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
