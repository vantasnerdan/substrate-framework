"""Sector decomposition and stability accounting for the chromomagnetic background.

Conditional, unpromoted primitive (P229, issues #47 and #48). Decomposes the
one-loop fluctuation determinant around the SU(2) Savvidy background into its
physical sectors.  The one-loop sector accounting is exact; the structured
stability summary also reports conditional P229 two-loop/literature results
that are not part of the narrow C-CMV-001/C-CMV-003 promotion candidates.

Sectors (background Feynman gauge, one SU(2) Cartan background):

- transverse charged vector, s3 = +1: proper-time spin factor e^{+2 gB s};
  contains the n=0 tachyon (E^2 = p_z^2 - gB);
- transverse charged vector, s3 = -1: factor e^{-2 gB s};
- longitudinal/timelike charged vector, s3 = 0 (twofold): factor 2;
- complex Faddeev-Popov ghost: weight -2;
- neutral gluon (color a=3): B-independent, drops out of V(B) - V(0).

The exact cancellation (s3=0 pair) + (ghost) = 0 leaves 2 cosh(2 gB s); the
entire one-loop logarithm is carried by the two transverse spin sectors,
each contributing b2/2 = 11/6 to the heat-kernel coefficient.

SU(3) color decomposition: a Cartan background along lambda_3 has charged
roots with charges (1, 1/2, 1/2) in units of g H, derived from the canonical
``T_3=lambda_3/2`` generator supporting C-LIE-001; each root contributes one
SU(2)-like complex charged vector with its own q g H.

At one loop, every b > 0 has the n=0, s3=+1 tachyon and the exact imaginary
part magnitude b^2/(8 pi).  P229 also records a conditional two-loop result
and literature interpretation, but those require separate review and are not
promoted by exposing this module.

This module is conditional infrastructure; presence here implies no claim
promotion (AGENTS_START_HERE.md section 7).
"""

from __future__ import annotations

import sympy as sp

from .chromomagnetic_background import nielsen_olesen_imaginary_part
from .chromomagnetic_two_loop import B2_tadpole, b, g
from .su3 import fundamental_generators

# ---------------------------------------------------------------------------
# Per-sector proper-time factors
# ---------------------------------------------------------------------------

SECTOR_FACTORS = {
    "transverse_spin_up": ("exp(+2 gB s)", 1),
    "transverse_spin_down": ("exp(-2 gB s)", 1),
    "longitudinal_pair": ("2", 2),
    "ghost": ("-2", -2),
}


def sector_spin_factor(sector: str, s, gB):
    """Proper-time spin factor of one sector (multiplies the orbital kernel)."""
    if sector == "transverse_spin_up":
        return sp.exp(2 * gB * s)
    if sector == "transverse_spin_down":
        return sp.exp(-2 * gB * s)
    if sector == "longitudinal_pair":
        return sp.Integer(2)
    if sector == "ghost":
        return sp.Integer(-2)
    raise ValueError(f"unknown sector {sector}")


def sector_heat_kernel_coefficients(order: int = 4) -> dict:
    """Exact small-s (gB)^2 s coefficient of gB/sinh(gB s) per sector."""
    s, gB = sp.symbols("s gB", positive=True)
    out = {}
    for sector in SECTOR_FACTORS:
        series = sp.series(
            gB / sp.sinh(gB * s) * sector_spin_factor(sector, s, gB), gB, 0, order + 2
        )
        series = sp.expand(series.removeO())
        out[sector] = sp.Rational(sp.simplify(series.coeff(gB, 2) / s))
    return out


def transverse_share_of_log() -> sp.Rational:
    """Fraction of the one-loop log carried by the transverse spin sectors."""
    coeffs = sector_heat_kernel_coefficients()
    total = sum(coeffs.values())
    transverse = coeffs["transverse_spin_up"] + coeffs["transverse_spin_down"]
    return sp.Rational(sp.simplify(transverse / total))


# ---------------------------------------------------------------------------
# SU(3) color decomposition
# ---------------------------------------------------------------------------


def su3_charged_roots() -> tuple:
    """Charges of the off-diagonal SU(3) roots w.r.t. a lambda_3 background.

    The values are computed from the canonical ``T_3=lambda_3/2`` matrix:
    the E_12 root sees the full field, while E_13 and E_23 see half.
    """
    t3 = fundamental_generators()[2]
    pairs = ((0, 1), (0, 2), (1, 2))
    return tuple(
        sp.simplify(abs(t3[first, first] - t3[second, second]))
        for first, second in pairs
    )


def su3_one_loop_potential(gH, mu2):
    """SU(3) one-loop potential in a lambda_3 background: one SU(2)-like
    complex charged vector per root with charge q, V_q = C (q gH)^2
    [ln(q gH / mu2) - 1/2]."""
    from substrate_framework.chromomagnetic_background import one_loop_log_coefficient

    c = one_loop_log_coefficient()
    return sum(
        c * (q * gH) ** 2 * (sp.log(q * gH / mu2) - sp.Rational(1, 2))
        for q in su3_charged_roots()
    )


# ---------------------------------------------------------------------------
# Stability accounting
# ---------------------------------------------------------------------------


def tachyon_condition(gB) -> sp.Expr:
    """The n=0, s3=+1 level: E^2 = p_z^2 - gB; tachyonic for p_z^2 < gB."""
    return sp.Symbol("p_z") ** 2 - gB


def im_v_one_loop(gB):
    """One-loop imaginary part (Nielsen-Olesen), magnitude b^2/(8 pi)."""
    return nielsen_olesen_imaginary_part(gB)


def im_v_two_loop():
    """Two-loop imaginary part from (3 g^2/2) B2^2 with B2 complex."""
    return sp.expand_complex(sp.Rational(3, 2) * g**2 * B2_tadpole() ** 2).as_real_imag()[1]


def stability_verdict() -> dict:
    """Return P229's mixed one-/two-loop stability summary.

    This reporting helper includes conditional two-loop literature results.
    It is not the verifier for the narrower one-loop promotion candidates.
    """
    return {
        "classical_minimum_exists": True,
        "fluctuation_tachyon": "E^2 = p_z^2 - gB < 0 for p_z^2 < gB, every gB > 0",
        "im_v_one_loop": "b^2/(8 pi) != 0",
        "im_v_two_loop": "-3 g^2 b^2 ln 2/(256 pi^3) != 0",
        "ring_resummation": "removes the one-loop Im part, insufficient at two loops (arXiv:2112.01043)",
        "verdict": "unstable against fluctuations; no perturbative stability at two-loop order",
    }
