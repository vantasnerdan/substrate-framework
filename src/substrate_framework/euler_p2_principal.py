"""Exact principal-symbol helpers for localized Euler P2 carrier tests.

The functions here expose algebraic identities only.  They do not assert a
generator domain, a finite-curvature Cao expansion, nonlinear persistence, or
a coercive Hill-vortex Hessian.  Their current scientific evidence is the
active P253/0095 attempt.
"""

from __future__ import annotations

import sympy as sp


def maxwell_transverse_block(c_frame, c_em, xi_z, xi_abs):
    """Return one polarization's comoving Maxwell 2x2 symbol."""
    return sp.Matrix(
        [
            [-sp.I * c_frame * xi_z, -sp.I * c_em**2 * xi_abs],
            [-sp.I * xi_abs, -sp.I * c_frame * xi_z],
        ]
    )


def column_bas_matrix(omega, zeta_col, k_z, k_abs):
    """Return the fixed-column BAS matrix in the ``(e_alpha,e_2)`` frame."""
    return sp.Matrix(
        [[0, -zeta_col * k_z / k_abs], [2 * omega * k_z / k_abs, 0]]
    )


def column_metric(omega, zeta_col, k_z_abs, k_abs):
    """Return the positive column metric for positive input coefficients."""
    return sp.diag(
        2 * omega * k_z_abs / k_abs,
        zeta_col * k_z_abs / k_abs,
    )


def metric_bent_c1(omega, omega_prime, radius, p, q, alpha):
    """Return the universal metric-only first-curvature BAS coefficient."""
    k_abs = sp.sqrt(p**2 + q**2)
    return sp.Matrix(
        [
            [
                -omega * radius * sp.sin(alpha),
                q
                * radius
                * (
                    5 * omega * p**2
                    + 3 * omega * q**2
                    + 2 * omega_prime * radius * p**2
                    + omega_prime * radius * q**2
                )
                * sp.cos(alpha)
                / k_abs**3,
            ],
            [
                -2
                * omega
                * q
                * radius
                * (2 * p**2 + q**2)
                * sp.cos(alpha)
                / k_abs**3,
                omega * p**2 * radius * sp.sin(alpha) / k_abs**2,
            ],
        ]
    )


def beltrami_hessian_symbols(alpha, xi_abs):
    """Return energy-helicity symbols on the two curl helicities."""
    return 1 - alpha * xi_abs, 1 + alpha * xi_abs


def cao_odd_border_beta(p, lane_emden_energy, lane_emden_mass):
    """Return the affine coefficient cancelling the odd translation moment."""
    return -(
        (p + 5) * lane_emden_energy
        / (2 * (p + 1) * lane_emden_mass)
    )


def bas_velocity_derivative(velocity_gradient, covector):
    """Differentiate the full-pressure BAS matrix with the covector fixed."""
    k2 = (covector.T * covector)[0]
    return -velocity_gradient + 2 * covector * covector.T * velocity_gradient / k2


def hf_hill_potential(a_log, a_log_prime, b_times_c):
    """Return Hattori--Fukumoto's Hill potential from the trace-free system."""
    return a_log_prime - a_log**2 - b_times_c


def hf_general_first_harmonic(
    omega_vorticity, omega_angular, wave_cosine_sq, circulation_ratio
):
    """Return the first harmonic in Hattori--Fukumoto equation (22)."""
    c2 = wave_cosine_sq
    s2 = 1 - c2
    return omega_angular + 2 * c2 * (
        (1 - 4 * s2) * omega_vorticity
        - omega_angular
        - circulation_ratio
        * ((2 * s2 - 1) * omega_vorticity + omega_angular / 2)
    )


def hf_growth_bracket(wave_cosine_sq, circulation_ratio):
    """Return the dimensionless first-resonance bracket in HF equation (24)."""
    s2 = 1 - wave_cosine_sq
    return (
        s2 / 2
        - sp.Rational(3, 8)
        - (sp.Rational(3, 8) - s2 / 4) * circulation_ratio
    )


def hf_resonant_slow_matrix(first_harmonic, oscillator_frequency):
    """Return the paired rotating-wave matrix for the resonant Hill equation."""
    coupling = first_harmonic / (4 * oscillator_frequency)
    return sp.Matrix([[0, coupling], [coupling, 0]])


def hf_detuned_return_generator(detuning, coupling):
    """Return the first Floquet generator after removing the resonant ``-I``."""
    return sp.Matrix([[sp.I * detuning, coupling], [coupling, -sp.I * detuning]])


def charged_column_slow_block(x_parallel, z_effective, r_effective):
    """Return the constrained frequency-last Euler--Lorentz column block."""
    return x_parallel * sp.Matrix([[0, -z_effective], [r_effective, 0]])
