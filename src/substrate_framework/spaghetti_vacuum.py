"""Flux-tube-lattice (spaghetti) vacuum vs the constant Savvidy background.

Conditional, unpromoted primitive (P229, issue #49). Implements the
Ambjorn-Olesen construction: the Nielsen-Olesen tachyon condenses into a
two-dimensional flux-tube lattice (lowest-Landau-level theta condensate),
the classical configuration closer to the quantum vacuum than the constant
chromomagnetic field.

Classical energy density after condensing the lowest tachyon mode W:

    E_cl/V = H^2/2 - gH |Phi|^2 + (g^2/2) beta |Phi|^4

with beta = <rho^2>/<rho>^2 the Abrikosov overlap parameter of the
condensate density rho = |W|^2. Minimizing over |Phi|^2 = H/(g beta):

    E_cl/V = (H^2/2)(1 - 1/beta),

strictly below the constant-field H^2/2 because beta > 1 for every lattice.
The triangular lattice minimizes beta.

The condensate density is evaluated convention-free: for lattice shape tau
(cell (1, tau)),

    rho(z) = exp(-2 pi (Im z)^2 / Im tau) * |theta_3(pi z | tau)|^2,

which is exactly periodic on (1, tau) (magnetic translations); beta(tau) is
the cell average ratio, computed numerically from this definition and
cross-checked against the literature values beta_square = 1.180340,
beta_triangular = 1.159595 (Abrikosov; review: Bordag, Symmetry 15 (2023)
1137, eqs. (146)-(152)).

Quantization around this configuration (the approximation used by the
literature): add the homogeneous-background one-loop potential of
chromomagnetic_background. The total effective potential

    V(H) = (1 - 1/beta) H^2/2 + C (gH)^2 [ln(gH/mu^2) - 1/2],  C = 11/(48 pi^2)

has ln(gH_min/mu^2) = -(1 - 1/beta) 24 pi^2/(11 g^2): the lattice vacuum is
deeper than Savvidy's, with the depth ratio enhanced exponentially at weak
coupling.

Convention note (documented, not resolved): the review's eq. (152) quotes
theta3(0, q0) ~= 1.2713 at q0 = e^{-pi}; under the standard nome convention
theta3(0, e^{-pi}) = 1.086434, while 1.2713 equals theta3(0, e^{-2}), so the
review's nome convention differs and its exact parameter elimination cannot
be reconstructed from the open text. The numbers here are derived from the
defining structures, not transcribed.

This module is conditional infrastructure; presence here implies no claim
promotion (AGENTS_START_HERE.md section 7).
"""

from __future__ import annotations

import numpy as np
import sympy as sp

# ---------------------------------------------------------------------------
# Abrikosov overlap parameter of the theta condensate
# ---------------------------------------------------------------------------


def condensate_density(z: complex, tau: complex) -> float:
    """rho(z) = exp(-2 pi (Im z)^2 / Im tau) |theta_3(pi z | tau)|^2.

    Exactly periodic on the lattice (1, tau); see module docstring.
    """
    import mpmath as mp

    q = mp.e ** (1j * mp.pi * tau)
    th = mp.jtheta(3, mp.pi * z, q)
    return float(mp.e ** (-2 * mp.pi * z.imag**2 / tau.imag) * abs(th) ** 2)


def abrikosov_beta(tau: complex, ngrid: int = 200) -> float:
    """beta(tau) = <rho^2>/<rho>^2 over the cell (1, tau), direct quadrature."""
    xs = np.linspace(0, 1, ngrid, endpoint=False)
    ys = np.linspace(0, abs(tau.imag), ngrid, endpoint=False)
    vals = np.empty((ngrid, ngrid))
    for i, x in enumerate(xs):
        for j, y in enumerate(ys):
            z = x + 1j * y
            # z in cell coordinates: point = x*1 + (y/Im tau)*tau
            zz = x + (y / tau.imag) * tau
            vals[i, j] = condensate_density(zz, tau)
    return float((vals**2).mean() / vals.mean() ** 2)


def beta_square(ngrid: int = 200) -> float:
    return abrikosov_beta(1j, ngrid)


def beta_triangular(ngrid: int = 200) -> float:
    return abrikosov_beta(np.exp(1j * np.pi / 3), ngrid)


def optimal_lattice() -> dict:
    """Compare the two classical candidates; triangular must win."""
    sq, tr = beta_square(), beta_triangular()
    return {
        "square": sq,
        "triangular": tr,
        "optimal": "triangular" if tr < sq else "square",
        "beta": min(sq, tr),
    }


# ---------------------------------------------------------------------------
# Classical energy and the total effective potential
# ---------------------------------------------------------------------------


def classical_energy_factor(beta: float) -> float:
    """E_cl/(H^2/2) = 1 - 1/beta after minimizing the condensate amplitude."""
    if beta <= 1:
        raise ValueError("beta > 1 required for condensation gain")
    return 1.0 - 1.0 / beta


def condensate_amplitude2(H: float, g: float, beta: float) -> float:
    """Minimizing condensate |Phi|^2 = H/(g beta)."""
    return H / (g * beta)


def spaghetti_potential(H, g, mu2, beta):
    """V(H) = (1 - 1/beta) H^2/2 + C (gH)^2 [ln(gH/mu^2) - 1/2]."""
    from substrate_framework.chromomagnetic_background import one_loop_log_coefficient

    c = one_loop_log_coefficient()
    return (1 - 1 / beta) * H**2 / 2 + c * (g * H) ** 2 * (
        sp.log(g * H / mu2) - sp.Rational(1, 2)
    )


def spaghetti_minimum(g, mu2, beta):
    """Minimum of the spaghetti potential.

    ln(g H_min/mu^2) = -(1 - 1/beta) 24 pi^2/(11 g^2);
    V_min = -C (gH_min)^2 / 2 < 0.
    """
    from substrate_framework.chromomagnetic_background import one_loop_log_coefficient

    c = one_loop_log_coefficient()
    f = 1 - 1 / beta
    ln_ratio = sp.simplify(-f / (2 * c * g**2))
    gH_min = mu2 * sp.exp(ln_ratio)
    v_min = sp.simplify(-c * gH_min**2 / 2)
    return ln_ratio, v_min


def depth_ratio_to_savvidy(beta) -> sp.Expr:
    """V_min(spaghetti)/V_min(Savvidy) at fixed mu, g:
    exp[(1/beta) 48 pi^2/(11 g^2)] (the exponent doubles through (gH_min)^2);
    exponentially deeper at weak coupling."""
    g_ = sp.Symbol("g", positive=True)
    return sp.exp((48 * sp.pi**2 / (11 * g_**2)) / beta)
