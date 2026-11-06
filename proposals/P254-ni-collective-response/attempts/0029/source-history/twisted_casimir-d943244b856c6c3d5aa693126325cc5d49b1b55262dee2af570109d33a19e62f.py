"""Conditional twisted-spectrum one-loop vacuum energetics on T^2 x R^2.

Authority status: conditional, unpromoted infrastructure linked to open goal
issue #26 (vantasnerdan/substrate-framework).  No accepted claim backs these
symbols; reuse requires declaring this status.

Geometry and conventions.  Two transverse coordinates are periodically
identified with side ``L`` (square torus T^2); the remaining two directions
(x0, x3) are noncompact.  A real scalar with twist ``alpha in R^2/Z^2`` obeys
psi(x + L e_i) = exp(2 pi i alpha_i) psi(x), so the transverse operator
-nabla_perp^2 has eigenvalues lambda_n(alpha) = (2 pi/L)^2 |n + alpha|^2,
``n in Z^2``.  The one-loop (Tr ln, Euclidean) density per real scalar is

    V1(alpha) = -(1/(8 pi L^2)) sum_n lambda_n [ln(lambda_n/mu^2) - 1].

Exact statements implemented here (proofs sketched in docstrings):

[T-A]  E2(-1; alpha) = 0 for every alpha, where E2(s; alpha) is the shifted
       Epstein zeta of the square lattice.  For alpha not in Z^2 this follows
       from the symmetric functional equation
       pi^-s Gamma(s) E2(s; alpha) = pi^(s-1) Gamma(1-s) Shat(1-s; alpha),
       whose right side is finite at s = -1 while Gamma(-1) diverges.
       For alpha in Z^2, E2(s) = 4 zeta(s) beta(s) and beta(-1) = 0.
       Consequence: the mu-dependence of V1 is proportional to E2(-1; alpha)
       and therefore vanishes identically, so
           V1(alpha) = (pi/(2 L^4)) E2'(-1; alpha).
       This is NOT a claim that an absolute vacuum density is
       scheme-independent: a local, alpha-independent counterterm can still
       shift V1 by an arbitrary C/L^4.  What is insensitive to such
       counterterms is the holonomy/twist DIFFERENCE V1(alpha) - V1(0),
       which is the only quantity this module asserts; the subtraction
       baseline is the periodic vacuum alpha = 0.
[T-B]  E2'(-1; alpha) = -S(alpha)/pi^3 for alpha not in Z^2, where
       S(alpha) = sum_{m != 0} cos(2 pi m . alpha)/|m|^4 (absolutely
       convergent).  At alpha = 0, E2'(-1; 0) = -beta'(-1)/3 with beta the
       Dirichlet beta function; the classical identity beta'(-1) = 2 G/pi
       (G = Catalan) is used as numeric evidence, cross-checked by the
       independent regulated mode-sum route.
"""

from __future__ import annotations

from typing import Any, Sequence

import numpy as np
import sympy as sp
from mpmath import mp, mpf, pi as mp_pi

from .su2_doublets import su2_fundamental_ledger

# ---------------------------------------------------------------------------
# special functions
# ---------------------------------------------------------------------------


def dirichlet_beta(s: Any) -> mpf:
    """Dirichlet beta via the Hurwitz split 4^-s (zeta(s,1/4)-zeta(s,3/4))."""

    s = mpf(s)
    return mp.power(4, -s) * (mp.zeta(s, mpf(1) / 4) - mp.zeta(s, mpf(3) / 4))


def _as_alpha(alpha: Sequence[float]) -> tuple[mpf, mpf]:
    if len(alpha) != 2:
        raise ValueError("alpha must be a 2-vector")
    return mpf(alpha[0]), mpf(alpha[1])


def _is_integral(alpha: Sequence[float]) -> bool:
    a1, a2 = _as_alpha(alpha)
    return mp.almosteq(a1, mp.nint(a1)) and mp.almosteq(a2, mp.nint(a2))


def _lipschitz_line(u: mpf, a: mpf) -> mpf:
    """sum_{m in Z} cos(2 pi m a)/(m^2 + u^2)^2 for u > 0 (closed form).

    Derived from the Lipschitz summation formula
        sum_m e^{2 pi i m a}/(m^2+u^2) = (pi/u) cosh(pi u (1-2 abar))/sinh(pi u),
    abar = frac(a) in [0,1), by applying -(1/(2u)) d/du.
    """

    abar = a - mp.floor(a)
    k = mp_pi * (1 - 2 * abar)
    sh, ch = mp.sinh(mp_pi * u), mp.cosh(mp_pi * u)
    sk, ck = mp.sinh(k * u), mp.cosh(k * u)
    # d/du [pi cosh(k u) / (u sinh(pi u))]
    g_prime = (
        mp_pi
        * (k * sk * u * sh - ck * (sh + mp_pi * u * ch))
        / (u**2 * sh**2)
    )
    return -g_prime / (2 * u)


def dual_s_sum(alpha: Sequence[float], t: Any = 2, *, tol: Any = None) -> mpf:
    """Return Shat(t; alpha) = sum_{m != 0} cos(2 pi m . alpha) / |m|^{2t}.

    Absolutely convergent for Re t > 1.  The m1 = 0 line is the 1D periodic
    zeta (Bernoulli closed form at t = 2); each m1 != 0 line uses the
    Lipschitz closed form at t = 2 (exponential convergence in m1) and a
    direct extrapolated sum otherwise.
    """

    a1, a2 = _as_alpha(alpha)
    t = mpf(t)
    if not t > 1:
        raise ValueError("dual_s_sum requires Re t > 1")
    tol = mpf(tol) if tol is not None else mp.eps * mpf(10) ** 8
    a2bar = a2 - mp.floor(a2)
    if _is_integral((a1, a2)):
        # exact: E2(2t ... ); for the square lattice S(0,0) at general t>1 is
        # E2(t; 0) = 4 zeta(t) beta(t).
        return 4 * mp.zeta(t) * dirichlet_beta(t)

    if t == 2:
        # sum_{m2 != 0} cos(2 pi m2 a2)/m2^4 = -(2 pi)^4 B4(a2bar)/24
        b4 = a2bar**4 - 2 * a2bar**3 + a2bar**2 - mpf(1) / 30
        total = -((2 * mp_pi) ** 4) * b4 / 24

        def line(m1: int) -> mpf:
            return 2 * mp.cos(2 * mp_pi * m1 * a1) * _lipschitz_line(mpf(m1), a2)

    else:

        def axis_term(m2: int, m1: int) -> mpf:
            return mp.cos(2 * mp_pi * m2 * a2bar) / (mpf(m1) ** 2 + m2**2) ** t

        def zline(m1: int) -> mpf:
            # full m2 line; for m1 = 0 the m2 = 0 term is excluded below
            start = 1 if m1 == 0 else 0
            head = mpf(0) if m1 == 0 else axis_term(0, m1)
            return head + 2 * mp.nsum(
                lambda m2: axis_term(m2, m1), [max(start, 1), mp.inf], tol=tol
            )

        total = zline(0)

        def line(m1: int) -> mpf:
            return 2 * mp.cos(2 * mp_pi * m1 * a1) * zline(m1)

    return total + mp.nsum(line, [1, mp.inf], tol=tol)


def epstein_square(s: Any, alpha: Sequence[float], *, max_mode: int = 400) -> mpf:
    """Shifted square-lattice Epstein zeta E2(s; alpha).

    alpha integral: exact factorization 4 zeta(s) beta(s).
    alpha non-integral, Re s < 0: functional equation
        E2(s) = pi^(2s-1) Gamma(1-s) rgamma(s) Shat(1-s; alpha)
    (rgamma = 1/Gamma, entire).  Re s > 1: direct lattice sum.
    The pole strip 0 <= Re s <= 1 is not implemented.
    """

    s = mpf(s)
    a1, a2 = _as_alpha(alpha)
    if _is_integral((a1, a2)):
        return 4 * mp.zeta(s) * dirichlet_beta(s)
    if s < 0:
        return (
            mp.power(mp_pi, 2 * s - 1)
            * mp.gamma(1 - s)
            * mp.rgamma(s)
            * dual_s_sum((a1, a2), 1 - s)
        )
    if s > 1:
        total = mpf(0)
        for n1 in range(-max_mode, max_mode + 1):
            x = n1 + a1
            for n2 in range(-max_mode, max_mode + 1):
                y = n2 + a2
                r2 = x * x + y * y
                if r2 == 0:
                    continue
                total += mp.power(r2, -s)
        return total
    raise ValueError("pole strip 0 <= Re s <= 1 not implemented")


def epstein_square_derivative_minus_one(alpha: Sequence[float]) -> mpf:
    """E2'(-1; alpha): -S(alpha)/pi^3 off the lattice; -beta'(-1)/3 on it."""

    a1, a2 = _as_alpha(alpha)
    if _is_integral((a1, a2)):
        # d/ds [4 zeta beta] at -1 = 4 zeta(-1) beta'(-1) since beta(-1) = 0.
        return 4 * mp.zeta(-1) * mp.diff(dirichlet_beta, -1)
    return -dual_s_sum((a1, a2), 2) / mp_pi**3


def epstein_square_value_minus_one(alpha: Sequence[float]) -> mpf:
    """E2(-1; alpha); [T-A] states this is identically zero."""

    return epstein_square(-1, alpha)


# ---------------------------------------------------------------------------
# one-loop densities
# ---------------------------------------------------------------------------


def one_loop_density_scalar(alpha: Sequence[float], L: Any = 1) -> mpf:
    """Zeta-regularized one-loop density (pi/(2 L^4)) E2'(-1; alpha).

    Exact mu-independence follows from [T-A]: the mu-dependent piece of the
    Tr ln density is proportional to E2(-1; alpha) = 0.
    """

    L = mpf(L)
    if not L > 0:
        raise ValueError("L must be positive")
    return (mp_pi / (2 * L**4)) * epstein_square_derivative_minus_one(alpha)


def vacuum_energy_difference(
    twists: Sequence[Sequence[float]],
    L: Any = 1,
    *,
    polarizations: int = 2,
) -> mpf:
    """Twisted-sector minus periodic one-loop density.

    polarizations x sum_a [V1(alpha_a) - V1(0)].  For a gauge field on a flat
    background the Faddeev-Popov cancellation leaves D-2 = 2 real scalars per
    adjoint component.
    """

    if polarizations < 1:
        raise ValueError("polarizations must be >= 1")
    reference = one_loop_density_scalar((0, 0), L)
    return polarizations * sum(
        one_loop_density_scalar(alpha, L) - reference for alpha in twists
    )


def regulated_mode_sum_difference(
    twists: Sequence[Sequence[float]],
    regulator: float,
) -> float:
    """Route 2: direct Gaussian-regulated spectral-moment difference.

    Returns D(Lambda) = sum_a U(alpha_a) - N_twists U(0) with
    U(alpha) = sum_n lam lam-log, lam = |n+alpha|^2, weight exp(-lam/Lambda).
    Compare against -(2/pi) L^4 DeltaV / polarizations from route 1:
    DeltaV = -(pi/L^4) D for polarizations = 2 (per scalar -pi/(2 L^4) D).
    """

    if not regulator > 0:
        raise ValueError("regulator must be positive")
    extent = int(6 * np.sqrt(regulator)) + 2
    grid = np.arange(-extent, extent + 1, dtype=float)

    def moment(alpha: Sequence[float]) -> float:
        a1, a2 = float(alpha[0]), float(alpha[1])
        x = grid + a1
        y = grid + a2
        r2 = x[:, None] ** 2 + y[None, :] ** 2
        mask = r2 > 0
        vals = np.where(mask, r2 * np.log(np.where(mask, r2, 1.0)), 0.0)
        vals *= np.exp(-r2 / regulator)
        return float(vals.sum())

    reference = moment((0.0, 0.0))
    return float(sum(moment(alpha) for alpha in twists) - len(twists) * reference)


# ---------------------------------------------------------------------------
# preprint-specific constants and symbolic algebra
# ---------------------------------------------------------------------------

#: SU(2)_L adjoint twist set of the toronic background (preprint Eqs. 28-30).
TORON_ADJOINT_TWISTS = ((0.5, 0.0), (0.0, 0.5), (0.5, 0.5))


def lattice_adjoint_transverse_spectrum(n_side: int) -> np.ndarray:
    """Dimensionless adjoint transverse Laplacian spectrum, N x N lattice.

    Returns eigenvalues of -(L/2pi)^2 nabla_perp^2 with unit lattice spacing,
    sorted ascending, for three real components with the constant toronic
    sign flips Ad_P = diag(-1,-1,1) on the x1-seam and Ad_Q = diag(1,-1,-1)
    on the x2-seam.  The flat adjoint bundle exists, so this spectrum must
    reproduce the analytic twisted Epstein spectrum; it serves as the
    validation of the lattice discretization used elsewhere.
    """

    if n_side < 4:
        raise ValueError("n_side must be >= 4")
    n = n_side
    dim = 3
    seam1 = np.diag([-1.0, -1.0, 1.0])  # Ad_P: (T1,T2,T3) -> (-T1,-T2,T3)
    seam2 = np.diag([1.0, -1.0, -1.0])  # Ad_Q: (T1,T2,T3) -> (T1,-T2,-T3)
    size = n * n * dim
    lap = np.zeros((size, size), dtype=float)

    def flat(n1: int, n2: int, comp: int) -> int:
        return (n1 * n + n2) * dim + comp

    for n1 in range(n):
        for n2 in range(n):
            for mu in (0, 1):
                forward = ((n1 + 1) % n, n2) if mu == 0 else (n1, (n2 + 1) % n)
                if mu == 0:
                    link = seam1 if n1 == n - 1 else np.eye(dim)
                else:
                    link = seam2 if n2 == n - 1 else np.eye(dim)
                for a in range(dim):
                    i0 = flat(n1, n2, a)
                    j0 = flat(forward[0], forward[1], a)
                    lap[i0, i0] += 1
                    lap[j0, j0] += 1
                    for b in range(dim):
                        lap[i0, flat(forward[0], forward[1], b)] -= link[a, b]
                        lap[flat(forward[0], forward[1], b), i0] -= link[a, b]
    return np.sort(np.linalg.eigvalsh(lap).real) * (n / (2 * np.pi)) ** 2


def analytic_twist_spectrum(
    alpha: Sequence[float], max_mode: int
) -> np.ndarray:
    """Dimensionless |n + alpha|^2 eigenvalues, sorted ascending."""

    a1, a2 = float(alpha[0]), float(alpha[1])
    values = [
        (n1 + a1) ** 2 + (n2 + a2) ** 2
        for n1 in range(-max_mode, max_mode + 1)
        for n2 in range(-max_mode, max_mode + 1)
    ]
    return np.sort(np.array(values))


# ---------------------------------------------------------------------------
# verified verdict constants
# ---------------------------------------------------------------------------


#: Preprint Eq. (63)/(64): DeltaV = -(5 pi G/8) L_EW^-4 (claimed, negative).
PREPRINT_DELTA_V_COEFFICIENT = -5 * sp.pi * sp.Catalan / 8

#: Corrected coefficient from [T-A]/[T-B]: DeltaV = +(5 G/2) L_EW^-4.
CORRECTED_DELTA_V_COEFFICIENT = 5 * sp.Catalan / 2


def transition_matrices() -> tuple[sp.ImmutableMatrix, sp.ImmutableMatrix]:
    """Toron transition matrices P = i sigma3, Q = i sigma1 (exact)."""

    generators = su2_fundamental_ledger().generators  # T_a = sigma_a / 2
    return (2 * sp.I * generators[2], 2 * sp.I * generators[0])


def wilson_loop_commutator() -> sp.ImmutableMatrix:
    """Return P Q P^-1 Q^-1 exactly; the 't Hooft flux claim is -I."""

    p, q = transition_matrices()
    return sp.simplify(p * q * p.inv() * q.inv())


def matrix_commutant_basis(
    matrices: Sequence[sp.MatrixBase],
) -> tuple[sp.ImmutableMatrix, ...]:
    """Exact basis of {V : [V, M] = 0 for every M in matrices} (2x2).

    Solves the entrywise linear system [V, M] = 0 with
    V = [[a, b], [c, d]] directly; no vec/Kronecker ordering conventions
    are involved, so the result is independent of reshape choices.
    """

    a, b, c, d = sp.symbols("a b c d")
    unknown = sp.Matrix([[a, b], [c, d]])
    equations: list[sp.Expr] = []
    for matrix in matrices:
        m = sp.Matrix(matrix)
        if m.shape != (2, 2):
            raise ValueError("matrix_commutant_basis handles 2x2 matrices")
        for row in (unknown * m - m * unknown).tolist():
            equations.extend(sp.expand(entry) for entry in row)
    system, _ = sp.linear_eq_to_matrix(equations, [a, b, c, d])
    basis = []
    for vector in system.nullspace():
        basis.append(
            sp.ImmutableMatrix(2, 2, lambda i, j: sp.expand(vector[2 * i + j]))
        )
    return tuple(basis)


def adjoint_twists() -> tuple[tuple[sp.Rational, sp.Rational], ...]:
    """Adjoint twist vectors from exact conjugation by P and Q.

    Returns the twists of the T1, T2, T3 adjoint components under
    eta -> P eta P^-1 (direction 1) and eta -> Q eta Q^-1 (direction 2),
    as ((dir1, dir2)) per component.
    """

    p, q = transition_matrices()
    generators = su2_fundamental_ledger().generators

    def sign(u: sp.ImmutableMatrix, v: sp.ImmutableMatrix) -> sp.Rational:
        product = sp.simplify(u * v * u.inv())
        if product == v:
            return sp.Rational(0)
        if product == -v:
            return sp.Rational(1, 2)
        raise ValueError("conjugation is not a sign flip on this generator")

    return tuple(
        (sign(p, t), sign(q, t)) for t in generators
    )


def higgs_quartic_from_gap_ansatz() -> dict[str, sp.Expr]:
    """Exact coefficient matching for the preprint's Eq. (76) gap ansatz.

    DeltaV(v) = -c [mu0^2 + m_W^2(v)]^2 with m_W^2 = g^2 v^2/4 gives
    V_eff(v) = DeltaV(v) - DeltaV(0); matching
    V_eff = -mu2_eff (v^2/2) + lambda_eff (v^2/4 ... ) via H+H = v^2/2
    yields lambda_eff = -c g^4/4 (negative), to be compared with the
    preprint's asserted Eq. (83) lambda_eff = +c g^4/8.
    """

    c, g, v, mu0 = sp.symbols("c g v mu0", positive=True)
    m_w2 = g**2 * v**2 / 4
    delta_v = -c * (mu0**2 + m_w2) ** 2
    v_eff = sp.expand(delta_v - delta_v.subs(v, 0))
    mu2_eff, lambda_eff = sp.symbols("mu2_eff lambda_eff")
    target = -mu2_eff * v**2 / 2 + lambda_eff * v**4 / 4
    solution = sp.solve(
        [
            sp.Eq(v_eff.coeff(v, 2), target.coeff(v, 2)),
            sp.Eq(v_eff.coeff(v, 4), target.coeff(v, 4)),
        ],
        [mu2_eff, lambda_eff],
        dict=True,
    )
    if len(solution) != 1:
        raise ValueError("coefficient matching is not unique")
    match = solution[0]
    return {
        "mu2_eff": sp.simplify(match[mu2_eff]),
        "lambda_eff": sp.simplify(match[lambda_eff]),
        "preprint_lambda_eff": c * g**4 / 8,
    }
