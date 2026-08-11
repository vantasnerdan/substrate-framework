"""Charged fields in a covariantly constant chromomagnetic background.

Conditional, unpromoted primitive (P229, issue #45). Encodes the Yang-Mills
fluctuation spectrum around a constant chromomagnetic field (SU(2), background
in color direction a=3, magnitude B; the background enters only through the
covariant combination gB), in background Feynman gauge, and evaluates the
one-loop effective potential by two routes sharing one declared regulator:

- route A: Schwinger proper-time determinant in closed hyperbolic form;
- route B: the same proper-time integral written as an explicit Landau-level
  sum, validating degeneracy, spin, and ghost bookkeeping level by level.

Frozen conventions (proposal P229): mostly-plus signature; SU(2) with the
background in a single Cartan direction, so exactly one charged complex vector
(colors a=1,2) and the complex Faddeev-Popov ghost see the background; the
neutral (a=3) sector is B-independent and drops out of V(B) - V(0).

Spectrum (derived from the discretized operator in the test suite, not
asserted): charged vector levels E^2 = p_z^2 + gB(2n+1) - 2 gB s_3 with
s_3 in {+1, 0, 0, -1}; charged ghost levels E^2 = p_z^2 + gB(2n+1). The
n=0, s_3=+1 level E^2 = p_z^2 - gB is the Nielsen-Olesen unstable mode.

This module is conditional infrastructure; presence here implies no claim
promotion (AGENTS_START_HERE.md section 7).
"""

from __future__ import annotations

import numpy as np
import sympy as sp
from scipy.integrate import quad

# ---------------------------------------------------------------------------
# Spectrum encoders (algebraic claims; the operator discretization is the oracle)
# ---------------------------------------------------------------------------


def charged_vector_mass2(gB: float, n: int, s3: int, pz2: float = 0.0) -> float:
    """Squared eigenvalue of a charged vector mode: pz2 + gB(2n+1) - 2 gB s3."""
    if n < 0:
        raise ValueError("Landau level n must be >= 0")
    if s3 not in (+1, 0, -1):
        raise ValueError("s3 must be +1, 0, or -1")
    return pz2 + gB * (2 * n + 1) - 2.0 * gB * s3


def charged_ghost_mass2(gB: float, n: int, pz2: float = 0.0) -> float:
    """Squared eigenvalue of a charged ghost (complex scalar) mode."""
    if n < 0:
        raise ValueError("Landau level n must be >= 0")
    return pz2 + gB * (2 * n + 1)


def spin_generator_12() -> sp.Matrix:
    """Vector Lorentz generator (S_12)^mu_nu = -i(d^mu_1 d^2_nu - d^mu_2 d^1_nu)."""
    return sp.Matrix(
        4,
        4,
        lambda i, j: -sp.I * ((1 if (i, j) == (1, 2) else 0) - (1 if (i, j) == (2, 1) else 0)),
    )


def landau_orbital_operator(gB: float, L: float, npoints: int) -> np.ndarray:
    """Discretized -D_perp^2 in the p_y = 0 sector (Landau gauge) on [-L/2, L/2].

    Oracle only: the algebraic spectrum is the claim; this operator verifies it.
    """
    x = np.linspace(-L / 2, L / 2, npoints)
    dx = x[1] - x[0]
    d2 = (
        np.diag(-2.0 * np.ones(npoints))
        + np.diag(np.ones(npoints - 1), 1)
        + np.diag(np.ones(npoints - 1), -1)
    ) / dx**2
    return -d2 + np.diag((gB * x) ** 2)


def vector_fluctuation_operator(gB: float, L: float, npoints: int) -> np.ndarray:
    """Discretized charged-vector fluctuation operator (orbital + 2 gB S_12)."""
    horb = landau_orbital_operator(gB, L, npoints)
    s12 = np.array(spin_generator_12().tolist(), dtype=complex)
    return np.kron(horb, np.eye(4)) + np.kron(np.eye(npoints), 2.0 * gB * s12)


# ---------------------------------------------------------------------------
# Proper-time integrand: exact structure
# ---------------------------------------------------------------------------


def heat_kernel_series_coefficients(order: int = 4) -> dict:
    """Exact small-s coefficients of gB/sinh(gB s) * 2 cosh(2 gB s).

    Returns {2: b2, 4: b4, ...} with the (gB)^k s^{k-1} coefficients as
    Rationals.  b2 = 11/3 is the Seeley-DeWitt coefficient behind the
    one-loop logarithm (the SU(2) Cartan-background beta coefficient).
    """
    s, gB = sp.symbols("s gB", positive=True)
    series = sp.series(gB / sp.sinh(gB * s) * 2 * sp.cosh(2 * gB * s), gB, 0, order + 2)
    series = sp.expand(series.removeO())
    out = {}
    for k in range(2, order + 1, 2):
        coef = sp.simplify(series.coeff(gB, k) / s ** (k - 1))
        out[k] = sp.Rational(coef)
    return out


def net_propertime_integrand(s, gB):
    """One-loop proper-time integrand after the exact ghost cancellation.

    Per unit 4-volume, per charge sector the vector spin trace is
    2 cosh(2 gB s) + 2 and the complex ghost weighs -2; the neutral (+2)
    spin-0 vector components cancel the ghosts exactly.  V1(B) is minus the
    s-integral of this against the longitudinal heat kernel (4 pi s)^-1.
    """
    orbital_transverse = gB / (4 * sp.pi * sp.sinh(gB * s))
    spin_plus_ghost = 2 * sp.cosh(2 * gB * s)  # (2 cosh + 2) - 2
    return orbital_transverse * spin_plus_ghost / (4 * sp.pi * s)


# ---------------------------------------------------------------------------
# One-loop potential: numeric evaluation with declared UV cutoff
# ---------------------------------------------------------------------------

_B2 = float(heat_kernel_series_coefficients(4)[2])  # 11/3
_B4 = float(heat_kernel_series_coefficients(4)[4])  # 127/180


def _integrand_minus_vacuum(sv: float, gB: float) -> float:
    """K(B) - K(0) with K = transverse * spin / (16 pi^2 s), stable for all s."""
    x = gB * sv
    if x < 0.02:
        return _B2 * gB**2 / (16 * np.pi**2) + _B4 * gB**4 * sv**2 / (16 * np.pi**2)
    ex = np.exp(-2 * x)
    # gB/sinh(x) * 2 cosh(2x) = 2 gB e^{x} (1 + e^{-4x}) / (1 - e^{-2x})
    kb = 2 * gB * np.exp(x) * (1 + np.exp(-4 * x)) / (1 - ex) / (16 * np.pi**2 * sv)
    k0 = 2.0 / (16 * np.pi**2 * sv * sv)
    return kb - k0


def one_loop_potential_cutoff(gB: float, m2: float, eps: float) -> float:
    """V1(B; eps, m2) = -Integral_eps^inf ds/s e^{-s m2} [K(B) - K(0)].

    Declared regulators: proper-time UV cutoff eps, IR mass m2 > gB (required
    so the Nielsen-Olesen growth e^{+gB s} stays damped; the m2 -> gB limit is
    where the imaginary part opens up).
    """
    if m2 <= gB:
        raise ValueError("IR mass must exceed gB so the unstable mode stays damped")

    def f(sv: float) -> float:
        x = gB * sv
        if x < 0.02:
            diff = _B2 * gB**2 / (16 * np.pi**2) + _B4 * gB**4 * sv**2 / (16 * np.pi**2)
            return -np.exp(-sv * m2) * diff / sv
        # fold the IR damping into the growing piece to avoid 0 * inf
        exx = np.exp(-sv * (m2 - gB))
        kb = 2 * gB * exx * (1 + np.exp(-4 * x)) / (1 - np.exp(-2 * x)) / (16 * np.pi**2 * sv)
        k0 = 2.0 * np.exp(-sv * m2) / (16 * np.pi**2 * sv * sv)
        return -(kb - k0) / sv

    val, _ = quad(f, eps, np.inf, limit=2000, epsabs=1e-13, epsrel=1e-11)
    return val


def one_loop_potential_modesum_cutoff(gB: float, m2: float, eps: float, nmax: int = 400) -> float:
    """Route B: the same cutoff integral as an explicit Landau-level sum.

    V1 = -Int_eps ds/s e^{-s m2} (gB/2pi)/(4 pi s) Sum_n
         [e^{-s gB(2n-1)} + e^{-s gB(2n+3)} + 2 e^{-s gB(2n+1)} - 2 e^{-s gB(2n+1)}]
    with the B=0 subtraction K0 = 2/(16 pi^2 s^2).  The ghost cancellation is
    written out, not pre-simplified, so a bookkeeping slip changes the answer.
    """
    if m2 <= gB:
        raise ValueError("IR mass must exceed gB so the unstable mode stays damped")
    # the sum must reach n ~ 1/(gB eps); scale the truncation to the UV cutoff
    nmax = max(nmax, int(30.0 / (gB * eps)) + 50)
    n = np.arange(nmax)

    def f(sv: float) -> float:
        # fold the IR damping into every level exponential (tachyon grows)
        e_spin_up = np.exp(-sv * (gB * (2 * n + 1) - 2 * gB + m2))
        e_spin_down = np.exp(-sv * (gB * (2 * n + 1) + 2 * gB + m2))
        e_zero = np.exp(-sv * (gB * (2 * n + 1) + m2))
        vec = e_spin_up + e_spin_down + 2.0 * e_zero
        ghost = -2.0 * e_zero
        ssum = (gB / (2 * np.pi)) * float(np.sum(vec + ghost)) / (4 * np.pi * sv)
        k0 = 2.0 * np.exp(-sv * m2) / (16 * np.pi**2 * sv * sv)
        return -(ssum - k0) / sv

    val, _ = quad(f, eps, np.inf, limit=2000, epsabs=1e-13, epsrel=1e-11)
    return val


# ---------------------------------------------------------------------------
# Analytic one-loop potential (MS-bar) and the Savvidy minimum
# ---------------------------------------------------------------------------


def one_loop_log_coefficient() -> sp.Rational:
    """Exact coefficient C in V1(B) - V1(0) = C (gB)^2 [ln(gB/mu^2) - 1/2].

    C = b2 / (4 pi)^2 with b2 = 11/3: V = -Int ds/s (4 pi s)^-1 [K-K(0)] and
    the proper-time identity Int_eps ds/s e^{-s gB} = -gamma_E - ln(eps gB)
    map the small-s coefficient b2 onto +(b2/(4 pi)^2) (gB)^2 ln(gB).  Both
    the factor and the sign are pinned end-to-end by the cutoff-route oracle
    in the test suite.
    """
    return heat_kernel_series_coefficients(2)[2] / (4 * sp.pi) ** 2


def one_loop_potential_msbar(gB, mu2):
    """V1(B) - V1(0) in MS-bar: C (gB)^2 [ln(gB/mu2) - 1/2]."""
    c = one_loop_log_coefficient()
    return c * gB**2 * (sp.log(gB / mu2) - sp.Rational(1, 2))


def total_potential(B, g, mu2):
    """V(B) = B^2/2 + V1(gB)."""
    return B**2 / 2 + one_loop_potential_msbar(g * B, mu2)


def savvidy_minimum(g, mu2):
    """Stationary point of the total potential.

    dV/dB = B + 2 C g^2 B ln(gB/mu2) = 0  =>  ln(gB_min/mu2) = -1/(2 C g^2).
    Returns (ln_ratio, V_min) with V_min = -C (gB_min)^2 / 2 < 0.
    """
    c = one_loop_log_coefficient()
    ln_ratio = sp.simplify(-1 / (2 * c * g**2))
    gB_min = mu2 * sp.exp(ln_ratio)
    v_min = sp.simplify(-c * gB_min**2 / 2)
    return ln_ratio, v_min


def minimum_at_lambda_scale(g) -> sp.Expr:
    """Scheme-independent content: with one-loop running b0 = 22/3 (SU(2)),

    1/g^2(mu) = (b0/8 pi^2) ln(mu/Lambda), the minimum sits at gB_min = Lambda^2.
    Returns the exact exponent ratio ln(gB_min/mu2) / ln(Lambda^2/mu2) = 1.
    """
    mu, lam = sp.symbols("mu Lambda", positive=True)
    b0 = sp.Rational(22, 3)
    g2_inv = b0 / (8 * sp.pi**2) * sp.log(mu / lam)
    ln_ratio = sp.simplify(-1 / (2 * one_loop_log_coefficient()) * g2_inv)
    return sp.simplify(ln_ratio / sp.log(lam**2 / mu**2))


# ---------------------------------------------------------------------------
# Route C: full operator-zeta evaluation (no cutoff, no IR mass)
# ---------------------------------------------------------------------------


def _operator_zeta(s, x):
    """Zeta function of the ghost-cancelled one-loop operator, per unit volume.

    zeta(s) = (x/2pi) x^(1-s) [e^{i pi (1-s)} + 2^(1-s) (zeta(s-1,1/2) + zeta(s-1,3/2))]
              / (4 pi (s-1))
    The e^{i pi (1-s)} term is the n=0, s3=+1 tachyon level; it carries the
    imaginary part.  V1 = -zeta'(0) (complex-field normalization).
    """
    import mpmath as mp

    return (
        (x / (2 * mp.pi))
        * x ** (1 - s)
        * (
            mp.exp(1j * mp.pi * (1 - s))
            + 2 ** (1 - s) * (mp.zeta(s - 1, 0.5) + mp.zeta(s - 1, 1.5))
        )
        / (4 * mp.pi * (s - 1))
    )


def one_loop_potential_zeta(gB: float, dps: int = 40) -> complex:
    """Zeta-regularized one-loop potential V1 = -zeta'(0) at gB (complex).

    Re V1 is the zeta-scheme real potential; Im V1 is the Nielsen-Olesen
    decay term (sign convention: principal branch, positive for e^{+i pi}).
    """
    import mpmath as mp

    with mp.workdps(dps):
        return complex(-mp.diff(lambda s: _operator_zeta(s, gB), 0))


def zeta_scheme_constant() -> float:
    """The zeta-scheme constant c1 in Re V1 = C x^2 [ln x + c1].

    Measured from one_loop_potential_zeta; the MS-bar value -1/2 is the
    literature convention (Savvidy restatements recorded in P229).  The
    difference is the scheme conversion, not physics.
    """
    xs = np.array([0.4, 0.6, 0.8, 1.0, 1.3])
    re = np.array([one_loop_potential_zeta(float(x)).real for x in xs])
    A = np.column_stack([xs**2 * np.log(xs), xs**2, xs**4 * np.log(xs), xs**4])
    coef, *_ = np.linalg.lstsq(A, re, rcond=None)
    return float(coef[1] / coef[0])


# ---------------------------------------------------------------------------
# Unstable mode: imaginary part of the one-loop potential
# ---------------------------------------------------------------------------


def nielsen_olesen_imaginary_part(gB) -> sp.Expr:
    """Imaginary part of the one-loop potential from the n=0, s_3=+1 tachyon.

    Derived from the mode integral, not asserted: per unit volume the unstable
    level contributes (gB/2pi) Int_{pz^2 < gB} dpz/(2pi) (1/2) sqrt(gB-pz^2)
    per real degree of freedom of the charged pair (factor 2), and
    Int_{-a}^{a} dp sqrt(a^2-p^2)/(2pi) = a^2/4, giving
    |Im V1| = (gB/2pi) * (1/2) * gB/4 * 2 = (gB)^2/(8 pi).
    The physical sign is set by the i-epsilon prescription (Im V < 0 marks the
    decay width of the Savvidy configuration).  Matches the value quoted by
    Nielsen-Olesen (1978) and restated by Savvidy (eq. (12) of the 2026 EPJC
    review) — the P229 extraction records both.
    """
    pz = sp.symbols("pz", real=True)
    a = sp.sqrt(gB)
    integral = sp.integrate(sp.sqrt(a**2 - pz**2) / (2 * sp.pi), (pz, -a, a))
    return sp.simplify((gB / (2 * sp.pi)) * sp.Rational(1, 2) * integral * 2)
