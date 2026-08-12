"""Two-loop effective potential of SU(2) gluodynamics in a chromomagnetic background.

Conditional, unpromoted primitive (P229, issue #46). Reproduces and audits the
T=0 two-loop result of Bordag-Skalozub, EPJC 82 (2022) 390 (arXiv:2112.01043),
eqs. (56)-(58), from their defining mode sums (their eq. (3)) at a=0, Feynman
gauge xi=1.

Chain of derivation (every step oracle-checked in the test suite):

1. The two-loop contribution at A0=0 is W2 = (3 g^2/2) B2(0,b)^2, from (56)
   with the gauge term absent at xi=1.
2. B2(0,b) = (b/4pi) Sum_{n,sigma=±2} Int d^2k/(2pi)^2 1/(k^2 + b(2n+1+sigma) - i0).
   Zeta-regularized: the per-level 2D tadpole finite part is -ln(M^2/mu^2)/(4pi);
   the pole and the ln b and ln mu terms vanish identically because
   G(1) = 1 + zeta(0,1/2) + zeta(0,3/2) = 0 (tachyon level + spin sectors).
3. The remaining level sum is
   Sum_{n,sigma} ln(2n+1+sigma-i0) = ln 2 - i pi,
   from Sum_{m>=0} ln(2m+a) = ln 2 * zeta(0,a/2) - zeta'(0,a/2) with
   zeta'(0,a) = ln Gamma(a) - ln sqrt(2 pi).
   Hence B2(0,b) = -b (ln 2 - i pi)/(16 pi^2).

Findings recorded against the printed paper (each pinned by a test):
- F1: their (58) prints the exponent correction 3 ln^2(2)/(98 pi^2) in the
  exponential and 3 ln^2(2)/(88 pi^2) in the expansion; minimization of their
  own printed (57) gives 88 -- the 98 is a typo.
- F2: their (57) prints the one-loop imaginary part as -i b^2/(8 pi^2); the
  operator-zeta and the mode integral both give b^2/(8 pi) (Nielsen-Olesen),
  see chromomagnetic_background.nielsen_olesen_imaginary_part.
- F3: their (57) prints the two-loop term g^2 ln^2(2) b^2/(128 pi^4).  The
  derivation from their own defining sums gives (3 g^2/2) B2^2, whose real
  part is 3 g^2 b^2 (ln^2 2 - pi^2)/(512 pi^4) -- a factor 4/3 in the ln^2 2
  coefficient and an additional -pi^2 structure from the complex square.
  Independently reproduced by the P229 literature extraction.

This module is conditional infrastructure; presence here implies no claim
promotion (AGENTS_START_HERE.md section 7).
"""

from __future__ import annotations

import sympy as sp

b, g, mu2 = sp.symbols("b g mu2", positive=True)


# ---------------------------------------------------------------------------
# Level-sum identities (exact)
# ---------------------------------------------------------------------------


def zeta_zero_sum() -> sp.Rational:
    """zeta(0,1/2) + zeta(0,3/2) = -1: the spin-sector divergence cancellation."""
    return sp.zeta(0, sp.Rational(1, 2)) + sp.zeta(0, sp.Rational(3, 2))


def log_level_sum() -> sp.Expr:
    """Sum_{n,sigma=±2} ln(2n+1+sigma-i0) = ln 2 - i pi.

    Sigma=+2 sector: tachyon ln(-1-i0) = -i pi plus Sum_{m>=0} ln(2m+1);
    sigma=-2 sector: Sum_{m>=0} ln(2m+3).  Each Sum_{m>=0} ln(2m+a) =
    ln 2 * zeta(0,a/2) - zeta'(0,a/2).
    """

    def s_log(a):
        # zeta'(0, a) = ln Gamma(a) - ln sqrt(2 pi)  (Lerch identity; verified
        # numerically against mpmath differentiation in the test suite)
        zeta_prime = sp.log(sp.gamma(a / 2)) - sp.Rational(1, 2) * sp.log(2 * sp.pi)
        return sp.log(2) * sp.zeta(0, a / 2) - zeta_prime

    spin_up = -sp.I * sp.pi + s_log(sp.Integer(1))
    spin_down = s_log(sp.Integer(3))
    return sp.simplify(spin_up + spin_down)


def B2_tadpole():
    """The magnetic tadpole B2(0,b) = -b (ln 2 - i pi)/(16 pi^2) (exact)."""
    return -b * log_level_sum() / (16 * sp.pi**2)


def B2_tadpole_mpmath() -> complex:
    """Independent mpmath evaluation of G'(1)/(16 pi^2) per unit b.

    G(s) = e^{-i pi (1-s)} + 2^(1-s) (zeta(s-1,1/2)+zeta(s-1,3/2));
    B2/b = G'(1)/(16 pi^2).  Tachyon branch e^{-i pi(1-s)} here; the imaginary
    sign is the -i0 prescription convention and is tracked separately.
    """
    import mpmath as mp

    with mp.workdps(40):
        def G(s):
            return mp.exp(-1j * mp.pi * (1 - s)) + 2 ** (1 - s) * (
                mp.zeta(s - 1, 0.5) + mp.zeta(s - 1, 1.5)
            )

        return complex(mp.diff(G, 1) / (16 * mp.pi**2))


# ---------------------------------------------------------------------------
# Full T=0 potential and its minimum
# ---------------------------------------------------------------------------


def W_two_loop_derived():
    """Full T=0 effective potential with the derived two-loop term.

    W = b^2/(2 g^2) + 11 b^2/(48 pi^2) (ln(b/mu^2) - 1/2) + (3 g^2/2) B2^2.
    The one-loop imaginary part b^2/(8 pi) is documented in
    chromomagnetic_background and added here for completeness.
    """
    w1 = b**2 / (2 * g**2) + sp.Rational(11) * b**2 / (48 * sp.pi**2) * (
        sp.log(b / mu2) - sp.Rational(1, 2)
    )
    return w1 + sp.Rational(3, 2) * g**2 * B2_tadpole() ** 2 - sp.I * b**2 / (8 * sp.pi)


def W_two_loop_printed():
    """The potential as printed in arXiv:2112.01043 eq. (57) (for comparison)."""
    return (
        b**2 / (2 * g**2)
        + sp.Rational(11) * b**2 / (48 * sp.pi**2) * (sp.log(b / mu2) - sp.Rational(1, 2))
        - sp.I * b**2 / (8 * sp.pi**2)
        + g**2 * sp.log(2) ** 2 * b**2 / (128 * sp.pi**4)
    )


def re_W2_derived():
    """Re W2 (derived) = 3 g^2 b^2 (ln^2 2 - pi^2)/(512 pi^4)."""
    return sp.expand(sp.re(sp.Rational(3, 2) * g**2 * B2_tadpole() ** 2))


def re_W2_printed():
    """Re W2 (printed) = g^2 ln^2(2) b^2/(128 pi^4)."""
    return g**2 * sp.log(2) ** 2 * b**2 / (128 * sp.pi**4)


def minimum_exponent(w_total) -> sp.Expr:
    """The exponent of b_min = mu^2 e^X, solved perturbatively to O(g^2).

    dW/db = b/g^2 [1 + (11 g^2/(24 pi^2)) X + (g^4/b^2) * dW2/db-terms] = 0
    gives X = -24 pi^2/(11 g^2) + g^2 X1; the function returns g^2 X1 as the
    exact correction coefficient.  Derived by SymPy from the input potential.
    """
    x_sym = sp.symbols("X")
    w_real = sp.re(sp.expand_complex(w_total))
    w2 = w_real - (
        b**2 / (2 * g**2)
        + sp.Rational(11) * b**2 / (48 * sp.pi**2) * (sp.log(b / mu2) - sp.Rational(1, 2))
    )
    dw = sp.diff(w_real, b).subs(b, mu2 * sp.exp(x_sym))
    # factor b/g^2 and identify the bracket
    bracket = sp.simplify(dw / (mu2 * sp.exp(x_sym) / g**2))
    # bracket = 1 + (11 g^2/(24 pi^2)) X + g^4 * R; solve X = X0 + g^2 X1
    x0 = -24 * sp.pi**2 / (11 * g**2)
    resid = sp.series(bracket.subs(x_sym, x0), g, 0, 5).removeO().expand()
    # resid must be cancelled by (11 g^2/(24 pi^2)) g^2 X1 -> X1 = -resid*24 pi^2/(11 g^4)
    x1 = sp.simplify(-resid * 24 * sp.pi**2 / (11 * g**4))
    return sp.simplify(g**2 * x1)


# ---------------------------------------------------------------------------
# The potential as a function of b/Lambda^2 (RG-improved, one-loop running)
# ---------------------------------------------------------------------------


def rg_improved_potential(two_loop: str = "derived"):
    """Formal W/Lambda^4 expression with one-loop running substituted.

    One-loop running for SU(2), b0 = 22/3: 1/g^2(mu) = (b0/8 pi^2) ln(mu/Lambda).
    Setting mu^2 = b (RG improvement) turns the tree term into
    b^2/2 * 1/g^2 = (11/(48 pi^2)) b^2 ln(b/Lambda^2), so with x = b/Lambda^2:

        W/Lambda^4 = (11/(48 pi^2)) x^2 [ln x - 1/2] + (two-loop) x^2 / ln x

    with the two-loop coefficient 3 ln^2 2 /(176 pi^2) (printed) or
    9 (ln^2 2 - pi^2)/(704 pi^2) (derived) — see tests.  For a nonzero
    two-loop coefficient this expression is defined only for ``x>0, x!=1``
    and has opposite-sign infinite one-sided limits at the one-loop minimum
    ``x=1``.  It therefore supplies neither a global two-loop potential nor a
    corrected minimum.  The ``two_loop='none'`` expression retains the exact
    one-loop minimum at ``x=1`` in the declared Lambda scheme.
    """
    x = sp.symbols("x", positive=True)
    one_loop = sp.Rational(11, 1) / (48 * sp.pi**2) * x**2 * (sp.log(x) - sp.Rational(1, 2))
    if two_loop == "derived":
        coef = 9 * (sp.log(2) ** 2 - sp.pi**2) / (704 * sp.pi**2)
    elif two_loop == "printed":
        coef = 3 * sp.log(2) ** 2 / (176 * sp.pi**2)
    elif two_loop == "none":
        coef = sp.Integer(0)
    else:
        raise ValueError(two_loop)
    return one_loop + coef * x**2 / sp.log(x)
