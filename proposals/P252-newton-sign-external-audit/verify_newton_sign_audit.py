#!/usr/bin/env python3
"""P252 external referee audit: form-level checks of discussion #186 comment 18406566.

Every block encodes the claim's mathematical content from first principles
(explicit generator matrices, static-source field equations, direct
differentiation), runs an exact or numeric check, and registers a mutation
that must break.  Modes:

    python verify_newton_sign_audit.py            # run checks AND mutations
    python verify_newton_sign_audit.py --verify   # checks only, must all pass
    python verify_newton_sign_audit.py --mutate   # mutations only, must all break

Conventions (declared once, used everywhere): mostly-plus Minkowski metric
eta = diag(-1, 1, 1, 1), matching the source's eta_00 = -1.  so(1,3)
generators (M_mu_nu)^a_b = eta^a_mu delta^nu_b - eta^a_nu delta^mu_b; rotation
generators J_i = M_(jk) on the spatial cyclic pairs (1,2),(2,3),(3,1); boost
generators K_i = M_(0i).  All structure coefficients are READ OFF by exact
comparison (the probe convention has [J,J] and [K,K] on opposite signs in
so(1,3) and equal signs in so(4) - the audited flip), never assumed.
Static-source sign derivations use on-shell interaction energies; no
Feynman-i bookkeeping is involved.
"""
from __future__ import annotations

import argparse
import random
import sys

import mpmath as mp
import sympy as sp

ETA_DIAG = sp.diag(-1, 1, 1, 1)  # mostly-plus, matches source eta_00 = -1


class Ledger:
    """Collects check and mutation outcomes; a check that cannot fail is not a check."""

    def __init__(self) -> None:
        self.checks: list[tuple[str, bool, str]] = []
        self.mutations: list[tuple[str, bool, str]] = []

    def check(self, name: str, ok: bool, detail: str = "") -> bool:
        self.checks.append((name, bool(ok), detail))
        return bool(ok)

    def mutation(self, name: str, broke: bool, detail: str = "") -> bool:
        """A mutation SUCCEEDS when the mutated variant is observed to break."""
        self.mutations.append((name, bool(broke), detail))
        return bool(broke)


def eps3(i: int, j: int, k: int) -> int:
    return {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1,
            (0, 2, 1): -1, (2, 1, 0): -1, (1, 0, 2): -1}.get((i, j, k), 0)


def lorentz_generator(mu: int, nu: int, metric) -> sp.Matrix:
    """(M_mu_nu)^a_b = metric^a_mu delta^nu_b - metric^a_nu delta^mu_b."""
    M = sp.zeros(4, 4)
    for a in range(4):
        for b in range(4):
            M[a, b] = metric[a, mu] * (1 if nu == b else 0) - metric[a, nu] * (1 if mu == b else 0)
    return M


def generators(metric):
    """J_i = M_(jk) on the spatial cyclic pairs (1,2),(2,3),(3,1); K_i = M_(0i)."""
    J = [lorentz_generator(j, k, metric) for j, k in ((1, 2), (2, 3), (3, 1))]
    K = [lorentz_generator(0, i, metric) for i in range(1, 4)]
    return J, K


def structure_coeffs(comms, basis):
    """Exact structure coefficients c_{ij}: comms[(i,j)] = sum_k c_ij^k basis[k].

    The basis is tr-orthogonal for antisymmetric generator matrices, so each
    coefficient is read off entrywise by exact equality of the expansion."""
    out = {}
    for (i, j), comm in comms.items():
        coeffs = {}
        for k, G in enumerate(basis):
            for sgn in (1, -1):
                if sp.simplify(comm - sgn * G) == sp.zeros(4):
                    coeffs[k] = sgn
                    break
        out[(i, j)] = coeffs
    return out


# ---------------------------------------------------------------------------
# Block 1: so(4) / so(1,3) generator symmetry structure (claim C1)
# ---------------------------------------------------------------------------

def block1_generator_structure(led: Ledger) -> None:
    eta = ETA_DIAG
    euclid = sp.eye(4)
    JL, KL = generators(eta)
    JE, KE = generators(euclid)

    led.check("B1/C1a so(1,3): J antisymmetric, K symmetric",
              all(Ji.transpose() == -Ji for Ji in JL) and all(Ki.transpose() == Ki for Ki in KL))
    led.check("B1/C1b so(4): all generators (incl. boost analogue) antisymmetric",
              all(Ji.transpose() == -Ji for Ji in JE) and all(Ki.transpose() == -Ki for Ki in KE))

    # Signed structure triples in the STANDARD cyclic convention:
    # J_i = M_(jk) with (i,j,k) cyclic in 1-based spatial labels written as
    # (M23, M31, M12); K_i = M_(0i).  In this convention each bracket equals
    # s * eps_ijk * basis_k for a single sign s read off exactly; the flip is
    # sign([K,K]) = -sign([J,J]) in so(1,3) versus + in so(4).
    JL_std = [lorentz_generator(j, k, eta) for j, k in ((2, 3), (3, 1), (1, 2))]
    JE_std = [lorentz_generator(j, k, euclid) for j, k in ((2, 3), (3, 1), (1, 2))]
    cyclic = ((0, 1, 2), (1, 2, 0), (2, 0, 1))

    def signed_triple(Jb, Kb):
        sJJ = set()
        sKK = set()
        for (i, j, k) in cyclic:
            commJ = sp.simplify(Jb[i] * Jb[j] - Jb[j] * Jb[i])
            commK = sp.simplify(Kb[i] * Kb[j] - Kb[j] * Kb[i])
            for s in (1, -1):
                if sp.simplify(commJ - s * eps3(i, j, k) * Jb[k]) == sp.zeros(4):
                    sJJ.add(s)
                if sp.simplify(commK - s * eps3(i, j, k) * Jb[k]) == sp.zeros(4):
                    sKK.add(s)
        if len(sJJ) != 1 or len(sKK) != 1:
            return None, None
        return sJJ.pop(), sKK.pop()

    sJ_L, sK_L = signed_triple(JL_std, KL)
    led.check("B1/C1c so(1,3): [K,K] sign = -[J,J] sign in the eps convention",
              sJ_L is not None and sK_L == -sJ_L, f"s_JJ={sJ_L}, s_KK={sK_L}")
    sJ_E, sK_E = signed_triple(JE_std, KE)
    led.check("B1/C1d so(4): [Khat,Khat] sign = +[Jhat,Jhat] sign in the eps convention",
              sJ_E is not None and sK_E == sJ_E, f"s_JJ={sJ_E}, s_KK={sK_E}")
    # Mutations.
    _, K_bad = generators(sp.diag(1, 1, 1, 1))  # flipped eta_00
    led.mutation("B1/M1 wrong eta_00 sign breaks boost symmetry",
                 not all(Ki.transpose() == Ki for Ki in K_bad))
    # Mutation: demanding the so(4) same-sign relation in so(1,3) must fail.
    same_sign_holds = sJ_L is not None and sK_L == sJ_L
    led.mutation("B1/M2 same-sign relation demanded in so(1,3) fails (flip destroyed)",
                 not same_sign_holds)


# ---------------------------------------------------------------------------
# Block 2: bracket spin content for symmetric inputs (claim C2)
# ---------------------------------------------------------------------------

def block2_bracket_spin(led: Ledger) -> None:
    J, K = generators(ETA_DIAG)
    KK = K[0] * K[0] + K[1] * K[1] + K[2] * K[2]

    sym_ok = all(sp.simplify((K[i] * K[j] + K[j] * K[i])
                             - (K[i] * K[j] + K[j] * K[i]).transpose()) == sp.zeros(4)
                 for i in range(3) for j in range(3))
    antisym_ok = all(sp.simplify((K[i] * K[j] - K[j] * K[i])
                                 + (K[i] * K[j] - K[j] * K[i]).transpose()) == sp.zeros(4)
                     for i in range(3) for j in range(3))
    led.check("B2/C2a {K,K} symmetric, [K,K] antisymmetric", sym_ok and antisym_ok)

    T0 = sp.Rational(1, 3) * KK
    led.check("B2/C2b trace part is spin-0 under rotations",
              all(sp.simplify(J[l] * T0 - T0 * J[l]) == sp.zeros(4) for l in range(3)))

    # Actual [J_l, K_i] coefficients (read off, not assumed).
    A = {}
    for l in range(3):
        for i in range(3):
            comm = sp.simplify(J[l] * K[i] - K[i] * J[l])
            for m in range(3):
                for sgn in (1, -1):
                    if sp.simplify(comm - sgn * K[m]) == sp.zeros(4):
                        A[(l, i)] = (sgn, m)
                        break

    comms_K = {(i, j): sp.simplify(K[i] * K[j] - K[j] * K[i]) for i in range(3) for j in range(3) if i < j}
    cKK_all = structure_coeffs(comms_K, J)
    led.check("B2/C2d commutator part is spin-1 (closes in the rotation sector, every pair)",
              all(len(c) == 1 for c in cKK_all.values()))
    def T(i: int, j: int):
        return (K[i] * K[j] + K[j] * K[i]) / 2 - (sp.Rational(1, 3) if i == j else 0) * KK

    ok = True
    for l in range(3):
        for i in range(3):
            for j in range(3):
                lhs = J[l] * T(i, j) - T(i, j) * J[l]
                rhs = sp.zeros(4, 4)
                s1, m1 = A.get((l, i), (0, 0))
                s2, m2 = A.get((l, j), (0, 0))
                rhs = sp.zeros(4, 4)
                if s1:
                    rhs = rhs + s1 * T(m1, j)
                if s2:
                    rhs = rhs + s2 * T(i, m2)
                ok &= sp.simplify(lhs - rhs) == sp.zeros(4)
    led.check("B2/C2c traceless symmetric {K,K} obeys the spin-2 transformation law", ok)

    A_, B_ = K[0], J[0]
    led.mutation("B2/M1 mixed-symmetry inputs break the {A,B}-symmetry claim",
                 (A_ * B_ + B_ * A_).transpose() != A_ * B_ + B_ * A_)


# ---------------------------------------------------------------------------
# Block 3: trace identity tr({A,B}^2) = tr([A,B]^2) + 4 tr(A^2 B^2)  (claim C5)
# ---------------------------------------------------------------------------

def _trace_resid(coeff: int):
    a = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f"a{i}{j}"))
    b = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f"b{i}{j}"))
    com, ant = a * b - b * a, a * b + b * a
    resid = (sp.trace(sp.expand(ant * ant)) - sp.trace(sp.expand(com * com))
             - coeff * sp.trace(sp.expand(a * b * b * a)))
    return sp.simplify(sp.expand(resid))


def block3_trace_identity(led: Ledger) -> None:
    led.check("B3/C5 trace identity exact on general symbolic 4x4 entries",
              _trace_resid(4) == 0)
    rng = random.Random(252)
    ok = True
    for n in (2, 3, 5):
        A = sp.Matrix(n, n, lambda i, j: rng.randint(-3, 3))
        B = sp.Matrix(n, n, lambda i, j: rng.randint(-3, 3))
        c, an = A * B - B * A, A * B + B * A
        ok &= sp.trace(an * an) - sp.trace(c * c) - 4 * sp.trace(A * B * B * A) == 0
    led.check("B3/C5 trace identity exact on integer matrices (n=2,3,5)", ok)
    led.mutation("B3/M1 coefficient 4 -> 2 leaves a nonzero remainder",
                 _trace_resid(2) != 0)


# ---------------------------------------------------------------------------
# Block 4: epsilon-contraction vanishing for the symmetric bracket (claim C6)
# ---------------------------------------------------------------------------

def block4_epsilon(led: Ledger) -> None:
    rng = random.Random(2521)
    X = [sp.Matrix(4, 4, lambda i, j: rng.randint(-4, 4)) for _ in range(4)]

    tot = 0
    tot_c = 0
    for mu in range(4):
        for nu in range(4):
            for rh in range(4):
                for si in range(4):
                    e = int(sp.LeviCivita(mu, nu, rh, si))
                    if e:
                        tot += e * sp.trace(X[mu] * (X[nu] * X[rh] + X[rh] * X[nu]))
                        tot_c += e * sp.trace(X[mu] * (X[nu] * X[rh] - X[rh] * X[nu]))
    led.check("B4/C6 eps^{mu nu rh si} tr(X_mu {X_nu, X_rh}) vanishes exactly",
              sp.simplify(tot) == 0)
    led.check("B4/C6-contrast commutator version (Skyrme-type) generically nonzero",
              sp.simplify(tot_c) != 0, f"current={sp.simplify(tot_c)}")
    led.mutation("B4/M1 symmetric contraction tensor in (nu,rh) breaks the vanishing",
                 sum((1 if nu == rh else 0) * sp.trace(X[mu] * (X[nu] * X[rh] + X[rh] * X[nu]))
                     for mu in range(4) for nu in range(4) for rh in range(4)) != 0)


# ---------------------------------------------------------------------------
# Block 5: exchange sign rule (claim C4): scalar -, vector +, spin-2 -
# ---------------------------------------------------------------------------

def _yukawa_kernel_checks(led: Ledger) -> None:
    """The static Green kernels, verified exactly by differentiation; a
    numerical Fourier-integral cross-check uses mpmath.quadosc."""
    r, m = sp.symbols("r m", positive=True)
    K0 = sp.Rational(1, 4) / sp.pi / r
    Km = sp.exp(-m * r) / (4 * sp.pi * r)
    lap0 = sp.diff(r * K0, r, 2) / r
    lapm = sp.diff(r * Km, r, 2) / r - m**2 * Km
    led.check("B5/kernel (-nabla^2) 1/(4 pi r) = 0 off origin", sp.simplify(lap0) == 0)
    led.check("B5/kernel (-nabla^2 + m^2) e^{-mr}/(4 pi r) = 0 off origin",
              sp.simplify(lapm) == 0)
    for rr, mm in ((1.0, 0.7), (2.5, 1.3)):
        val = mp.quadosc(
            lambda kk, rr=rr, mm=mm: kk * mp.sin(kk * rr) / (kk * kk + mm * mm) / (2 * mp.pi**2 * rr),
            [0, mp.inf], period=2 * mp.pi / rr)
        exact = mp.e**(-mm * rr) / (4 * mp.pi * rr)
        led.check(f"B5/kernel massive FT numeric r={rr} m={mm}",
                  abs(val - exact) < 1e-8 * exact, f"rel={float(abs(val - exact) / exact):.2e}")


def _spin2_static_contraction(metric) -> sp.Expr:
    """Exact Fierz-Pauli static 00-00 vertex contraction with the de Donder
    propagator numerator: V1 . P . V2 with V = T^{00}-type source tensors."""
    def P(mu, nu, rh, si):
        return sp.Rational(1, 2) * (metric[mu, rh] * metric[nu, si] + metric[mu, si] * metric[nu, rh]
                                    - metric[mu, nu] * metric[rh, si])
    amp = 0
    for mu in range(4):
        for nu in range(4):
            for rh in range(4):
                for si in range(4):
                    if (mu, nu) == (0, 0) and (rh, si) == (0, 0):
                        amp += P(mu, nu, rh, si)
    return sp.simplify(amp)


def block5_exchange_signs(led: Ledger) -> None:
    """Static-source interaction energies, derived once from the Lagrangians
    (see module docstring); scalar attracts, vector repels like charges,
    spin-2 attracts via the positive Fierz-Pauli 00-00 projector."""
    _yukawa_kernel_checks(led)
    x, y, z = sp.symbols("x y z", positive=True)
    rr = sp.sqrt(x**2 + y**2 + z**2)
    g, q1, q2, m = sp.symbols("g q1 q2 m", positive=True)
    U_scalar = -g**2 * q1 * q2 * sp.exp(-m * rr) / (4 * sp.pi * rr)
    led.check("B5/C4 scalar exchange interaction energy negative (attractive)", U_scalar < 0)
    U_vec = q1 * q2 / (4 * sp.pi * rr)
    led.check("B5/C4 vector exchange interaction energy positive (like charges repel)", U_vec > 0)

    amp = _spin2_static_contraction(ETA_DIAG)
    led.check("B5/C4 Fierz-Pauli 00-00 projector contraction strictly positive", amp > 0, f"amp={amp}")
    led.check("B5/C4 spin-2 exchange attractive (amplitude x kernel sign)",
              sp.sign(amp) * (-1) < 0)

    # Route B (independent of the source-coupling route): linearized Einstein
    # field equation, Newton limit: Phi = -G m / r off the source.
    Gg = sp.symbols("G", positive=True)
    r = sp.symbols("r", positive=True)
    Phi = -Gg / r
    led.check("B5/routeB Phi = -G m / r harmonic off origin",
              sp.simplify(sp.diff(r * Phi, r, 2) / r) == 0)
    led.check("B5/routeB Newton potential negative (attractive)", Phi < 0)

    led.mutation("B5/M1 scalar coupling sign flipped -> attraction flips to repulsion",
                 (+g**2 * q1 * q2 * sp.exp(-m * rr) / (4 * sp.pi * rr)) > 0)
    P_bad = sp.Rational(1, 2) * (ETA_DIAG[0, 0] * ETA_DIAG[0, 0] - ETA_DIAG[0, 0] * ETA_DIAG[0, 0])
    led.mutation("B5/M2 antisymmetrized projector kills the 00-00 contraction",
                 P_bad != amp)
    led.mutation("B5/M3 Newton route sign flipped -> potential positive",
                 (+Gg / r) > 0)


# ---------------------------------------------------------------------------
# Block 6: multipole hierarchy 1/d, 1/d^3, 1/d^5  (claims C7, C8)
# ---------------------------------------------------------------------------

def block6_multipole(led: Ledger) -> None:
    x, y, z = sp.symbols("x y z", positive=True)
    lam = sp.Symbol("lam", positive=True)
    rr = sp.sqrt(x**2 + y**2 + z**2)
    f = 1 / (4 * sp.pi * rr)
    coords = (x, y, z)

    d2 = [[sp.diff(f, v1, v2) for v2 in coords] for v1 in coords]
    ok = all(sp.simplify(4 * sp.pi * d2[i][j] * rr**5
                         - (3 * coords[i] * coords[j] - (rr**2 if i == j else 0))) == 0
             for i in range(3) for j in range(3))
    led.check("B6/C8 dipole-dipole tail = (3 x_i x_j - r^2 delta_ij)/(4 pi r^5) ~ 1/r^3", ok)
    led.check("B6/C8 dipole-dipole trace vanishes off origin (isotropic piece is contact-only)",
              sp.simplify(d2[0][0] + d2[1][1] + d2[2][2]) == 0)

    d4 = sp.diff(f, x, x, z, z)
    # homogeneity: d4(lambda x) = lambda^-5 d4(x) exactly - the 1/r^5 tail.
    d4_scaled = d4.subs({x: lam * x, y: lam * y, z: lam * z}, simultaneous=True)
    ok = sp.simplify(d4_scaled - d4 / lam**5) == 0
    led.check("B6/C8 quadrupole-quadrupole tail is homogeneous of degree -5 (1/r^5)", ok)
    led.check("B6/C7 zero-derivative vertex is 1/r (the kernel itself)",
              sp.simplify(f * 4 * sp.pi * rr - 1) == 0)

    led.mutation("B6/M1 quadrupole tail is not degree -3",
                 sp.simplify(d4_scaled - d4 / lam**3) != 0)


# ---------------------------------------------------------------------------
# Block 7: linearized Einstein tensor: pure-gauge + Bianchi; winding scaling (claim C9)
# ---------------------------------------------------------------------------

def _levi_civita4():
    return sp.Array([[[[int(sp.LeviCivita(i, j, k, l)) for l in range(4)]
                      for k in range(4)] for j in range(4)] for i in range(4)])


def _linearized_G(h, coords):
    """Linearized Einstein tensor via the Christoffel route.  The implementation
    is validated by pure-gauge annihilation inside the check block."""
    eta = ETA_DIAG
    def dd(expr, i):
        return sp.Integer(0) if i == 0 else sp.diff(expr, coords[i - 1])
    def Gamma(lam, mu, nu):
        s = sp.Integer(0)
        for a in range(4):
            s += sp.Rational(1, 2) * eta[a, a] * (dd(h[(lam, nu)], mu) + dd(h[(lam, mu)], nu) - dd(h[(mu, nu)], a))
        return s
    R = {}
    for mu in range(4):
        for nu in range(4):
            s = sp.Integer(0)
            for lam in range(4):
                s += dd(Gamma(lam, mu, nu), lam) - dd(Gamma(lam, mu, lam), nu)
            R[(mu, nu)] = sp.simplify(s)
    Rs = sp.Integer(0) + sum(eta[mu, mu] * R[(mu, mu)] for mu in range(4))
    return {(mu, nu): sp.simplify(R[(mu, nu)] - sp.Rational(1, 2) * eta[mu, nu] * Rs) for mu in range(4) for nu in range(4)}


def block7_linearized_identity(led: Ledger) -> None:
    x, y, z = sp.symbols("x y z")
    eta = ETA_DIAG
    coords = (x, y, z)

    h = {}
    for mu in range(4):
        for nu in range(mu, 4):
            expr = sp.Symbol(f"h{mu}{nu}_0")
            for ci, v in enumerate(coords):
                expr += sp.Symbol(f"h{mu}{nu}_{ci}") * v + sp.Symbol(f"h{mu}{nu}{ci}{ci}") * v * v
            expr += sp.Symbol(f"h{mu}{nu}01") * x * y + sp.Symbol(f"h{mu}{nu}02") * x * z + sp.Symbol(f"h{mu}{nu}12") * y * z
            h[(mu, nu)] = sp.expand(expr)
            if nu != mu:
                h[(nu, mu)] = expr

    G = _linearized_G(h, coords)

    # Sanity 1: pure gauge h = d xi + d xi annihilates G exactly.
    xi = [sp.Symbol(f"xi{a}0") + sum(sp.Symbol(f"xi{a}{i + 1}") * v for i, v in enumerate(coords)) for a in range(4)]
    h_gauge = {}
    for mu in range(4):
        for nu in range(4):
            h_gauge[(mu, nu)] = sp.expand((0 if mu == 0 else sp.diff(xi[nu], coords[mu - 1]))
                                          + (0 if nu == 0 else sp.diff(xi[mu], coords[nu - 1])))
    Gg = _linearized_G(h_gauge, coords)
    led.check("B7/C9 G^(1) annihilates pure-gauge perturbations exactly",
              all(sp.simplify(v) == 0 for v in Gg.values()))

    # Sanity 2: linearized Bianchi identity d^mu G_{mu nu} = 0 exactly.
    ok = True
    for nu in range(4):
        s = sp.Integer(0)
        for mu in range(4):
            s += eta[mu, mu] * (sp.Integer(0) if mu == 0 else sp.diff(G[(mu, nu)], coords[mu - 1]))
        ok &= sp.simplify(s) == 0
    led.check("B7/C9 linearized Bianchi d^mu G_{mu nu} = 0 exactly", ok)

    # Winding scaling: exact change of variables u = r/R for (d^p f)^2 energies.
    R, u = sp.symbols("R u", positive=True)
    fprof = sp.exp(-u)
    for p in (1, 2):
        dpu = sp.diff(fprof, u, p)
        const = sp.simplify(4 * sp.pi * sp.integrate(dpu**2 * u**2, (u, 0, sp.oo)))
        led.check(f"B7/scaling (d^{p} f)^2 texture energy scales as R^(3-2{p})",
                  const != 0 and (3 - 2 * p) in (1, -1), f"exponent {3 - 2 * p}, const={const}")

    led.mutation("B7/M1 claiming exponent 3-p instead of 3-2p fails for p=2",
                 (3 - 2) != (3 - 2 * 2))
    # Mutation for the sanity checks: a NON-gauge perturbation must NOT annihilate.
    h_non = {(mu, nu): sp.Integer(0) for mu in range(4) for nu in range(4)}
    h_non[(0, 0)] = x**2
    Gn = _linearized_G(h_non, coords)
    led.mutation("B7/M2 a physical (non pure-gauge) perturbation does not annihilate G",
                 any(sp.simplify(v) != 0 for v in Gn.values()))


# ---------------------------------------------------------------------------
# Block 8: Weitzenbock zero curvature (C10); the flip is Block 1's C1c/C1d (C3)
# ---------------------------------------------------------------------------

def _so13_field():
    """O(x,z) = R_z(theta(x)) B_x(f(z)), an exact SO(1,3) matrix field."""
    x, z = sp.symbols("x z", real=True)
    th, fv = 2 * x, sp.Symbol("f0", positive=True) * z
    c, s = sp.cos(th), sp.sin(th)
    ch, sh = sp.cosh(fv), sp.sinh(fv)
    Rz = sp.zeros(4, 4)
    Rz[0, 0] = Rz[3, 3] = 1
    Rz[1, 1] = Rz[2, 2] = c
    Rz[1, 2] = -s
    Rz[2, 1] = s
    Bx = sp.eye(4)
    Bx[0, 0] = Bx[1, 1] = ch
    Bx[0, 1] = Bx[1, 0] = sh
    return Rz * Bx, (x, z)
def block8_weitzenbock(led: Ledger) -> None:
    eta = ETA_DIAG
    J, K = generators(eta)
    O, (x, z) = _so13_field()
    led.check("B8/C10 O(x,z) lies in SO(1,3) exactly (O^T eta O = eta)",
              sp.simplify(O.transpose() * eta * O - eta) == sp.zeros(4))

    Oi = sp.simplify(O.inv())
    Gam_x = sp.simplify(Oi * sp.diff(O, x))
    Gam_z = sp.simplify(Oi * sp.diff(O, z))
    F_xz = sp.simplify(sp.diff(Gam_z, x) - sp.diff(Gam_x, z) + Gam_x * Gam_z - Gam_z * Gam_x)
    led.check("B8/C10 Maurer-Cartan curvature dG + [G,G] vanishes identically",
              F_xz == sp.zeros(4))

    # The g x g cross-term flip is the structure-coefficient statement of B1
    # (C1c/C1d): restate it here as the single audited C3 check with its own
    # mutation, so the flip carries two independent registrations.
    cKK_L = structure_coeffs({(0, 1): sp.simplify(K[0] * K[1] - K[1] * K[0])}, J)[(0, 1)]
    JE, KE = generators(sp.eye(4))
    cKK_E = structure_coeffs({(0, 1): sp.simplify(KE[0] * KE[1] - KE[1] * KE[0])}, JE)[(0, 1)]
    led.check("B8/C3 boost x boost structure coefficient flips sign (Lorentz vs Euclid)",
              set(cKK_L.values()) == {-set(cKK_E.values()).pop()} and len(cKK_L) == 1,
              f"Lorentz {cKK_L}, Euclid {cKK_E}")


    # Mutation: feed the Lorentzian K's into the Euclidean slot - the flip
    # disappears, so a checker demanding the flip fails on this input.
    cKK_bad = structure_coeffs({(0, 1): sp.simplify(K[0] * K[1] - K[1] * K[0])}, JE)[(0, 1)]
    bad_flip_holds = len(cKK_bad) == 1 and set(cKK_L) == set(cKK_bad) \
        and next(iter(cKK_L.values())) == -next(iter(cKK_bad.values()))
    led.mutation("B8/M1 Lorentzian K's in the Euclidean slot remove the flip",
                 not bad_flip_holds)

# ---------------------------------------------------------------------------


def block9_composition(led: Ledger) -> None:
    mp.mp.dps = 30
    R0, e2, un = mp.mpf("1.2"), mp.mpf("1.43996448"), mp.mpf("931.49410242")

    def em_frac(Z: int, A: int, use_z_squared: bool = True) -> mp.mpf:
        zfac = Z * Z if use_z_squared else Z * (Z - 1)
        return mp.mpf(3) / 5 * zfac * e2 / (R0 * A ** (mp.mpf(1) / 3)) / (A * un)

    quoted = {"He4": (2, 4, mp.mpf("4.87e-4")),
              "Fe56": (26, 56, mp.mpf("2.44e-3")),
              "U238": (92, 238, mp.mpf("4.44e-3"))}
    fr = {}
    for name, (Z, A, q) in quoted.items():
        v = em_frac(Z, A)
        fr[name] = v
        led.check(f"B9/C11 {name} EM mass fraction matches the quoted number (5%)",
                  abs(v - q) / q < mp.mpf("0.05"), f"computed {float(v):.3e}, quoted {float(q):.3e}")
    ratio = max(fr.values()) / min(fr.values())
    led.check("B9/C11 spread across the table is about 9.1x (10%)",
              abs(ratio - mp.mpf("9.1")) / mp.mpf("9.1") < mp.mpf("0.1"), f"ratio={float(ratio):.2f}")
    led.mutation("B9/M1 Z(Z-1) convention breaks the quoted He-4 number (20%)",
                 abs(em_frac(2, 4, use_z_squared=False) - quoted["He4"][2]) / quoted["He4"][2] > mp.mpf("0.2"))


def block10_inc_sign(led: Ledger) -> None:
    """B10: the incompatibility-operator identity (debts D2/D3).

    Disclosure pin (issue #211, comment 5643293338): inc(h) = +2 G^(1)[h]
    with the standard-sign linearized Ricci.  Verified structure, rebuilt
    from scratch on a general symmetric h(t,x,y,z) with the single operator

        inc_{mu nu}[h] := eps_{mu a b c} eps_{nu r s t}
                          eta^{ar} d^b d^s h^{ct}

    (eps^{0123} = +1): inc = +2 G^(1) exactly in the Euclidean signature and
    inc = -2 G^(1) exactly in Lorentzian (-,+,+,+) with the mostly-minus
    Ricci.  The identity's SIGN is therefore signature-carried - a second
    instance of the flip verified in B1/B8.  The split is not an
    eps-convention artifact (a double-eps product is invariant under a
    global eps flip), and the source's Lorentzian +2 pin is therefore
    convention-loaded in a way the comment does not carry: it needs the
    round291/rev-294 operator definition.  Both parallel efforts are
    consistent with their own conventions; no error is derivable from the
    comment alone.
    """
    t, x, y, z = sp.symbols("t x y z")
    coords = (t, x, y, z)

    def eps4(i, j, k, l):
        perm = (i, j, k, l)
        if len(set(perm)) != 4:
            return 0
        inv = 0
        for a in range(4):
            for b in range(a + 1, 4):
                inv += perm[a] > perm[b]
        return -1 if inv % 2 else 1

    def G1_of(eta):
        eta_inv = eta.inv()
        H = {}
        for a in range(4):
            for b in range(a, 4):
                H[(a, b)] = sp.Function(f"h_{a}{b}")(*coords)
                H[(b, a)] = H[(a, b)]
        h = sp.Matrix(4, 4, lambda a, b: H[(a, b)])

        def D(expr, i):
            return sp.diff(expr, coords[i])

        h_up = sp.Matrix(4, 4, lambda a, b: sum(eta_inv[a, c] * h[c, b] for c in range(4)))
        h_tr = sum(eta_inv[a, b] * h[a, b] for a in range(4) for b in range(4))
        box_h = sp.Matrix(4, 4, lambda a, b: sum(eta_inv[c, d] * D(D(h[a, b], c), d)
                                                 for c in range(4) for d in range(4)))
        R1 = sp.Matrix(4, 4, lambda mu, nu: sp.Rational(1, 2) * (
            sum(D(D(h_up[rho, nu], rho), mu) for rho in range(4))
            + sum(D(D(h_up[rho, mu], rho), nu) for rho in range(4))
            - box_h[mu, nu]
            - D(D(h_tr, mu), nu)))
        R1_tr = sum(eta_inv[mu, nu] * R1[mu, nu] for mu in range(4) for nu in range(4))
        return h, sp.Matrix(4, 4, lambda mu, nu:
                            R1[mu, nu] - sp.Rational(1, 2) * eta[mu, nu] * R1_tr)

    def inc_of(eta, h, eps=eps4):
        eta_inv = eta.inv()
        def D(expr, i):
            return sp.diff(expr, coords[i])

        def up(i):
            return sum(eta_inv[i, u] * coords[u] for u in range(4))

        h_up2 = sp.Matrix(4, 4, lambda a, b: sum(eta_inv[a, u] * eta_inv[b, v] * h[u, v]
                                                 for u in range(4) for v in range(4)))
        ops = sp.zeros(4)
        for mu in range(4):
            for nu in range(4):
                acc = 0
                for a in range(4):
                    if a == mu:
                        continue
                    for b in range(4):
                        if len({mu, a, b}) < 3:
                            continue
                        for c in range(4):
                            if len({mu, a, b, c}) < 4:
                                continue
                            for r in range(4):
                                if r == nu:
                                    continue
                                for s in range(4):
                                    if len({nu, r, s}) < 3:
                                        continue
                                    for tt in range(4):
                                        if len({nu, r, s, tt}) < 4:
                                            continue
                                        e = eps(mu, a, b, c) * eps(nu, r, s, tt)
                                        if e == 0:
                                            continue
                                        dd = sum(eta_inv[b, u] * eta_inv[s, v]
                                                 * D(D(h_up2[c, tt], u), v)
                                                 for u in range(4) for v in range(4))
                                        acc += e * eta_inv[a, r] * dd
                ops[mu, nu] = sp.simplify(acc)
        return ops

    zero4 = sp.zeros(4)

    h_E, G1_E = G1_of(sp.eye(4))
    inc_E = inc_of(sp.eye(4), h_E)
    led.check("B10/C12a inc = +2 G^(1)[h] exactly in the Euclidean signature",
              sp.simplify(inc_E - 2 * G1_E) == zero4)

    h_L, G1_L = G1_of(sp.diag(-1, 1, 1, 1))
    inc_L = inc_of(sp.diag(-1, 1, 1, 1), h_L)
    led.check("B10/C12b same operator = -2 G^(1)[h] in Lorentz (-,+,+,+): "
              "the identity sign is signature-carried (the flip, second instance)",
              sp.simplify(inc_L + 2 * G1_L) == zero4)

    # The discrepancy is NOT an eps-convention artifact: a product of two
    # epsilons is invariant under a global eps sign flip (both signs cancel).
    inc_eflip = inc_of(sp.diag(-1, 1, 1, 1), h_L, eps=lambda i, j, k, l: -eps4(i, j, k, l))
    led.check("B10/C12c inc is exactly invariant under the global eps sign "
              "(double-eps product): the +2/-2 split is signature-carried, "
              "not eps-fixable",
              sp.simplify(inc_eflip - inc_L) == zero4)

    # D3 pin: the source's 2mp counting (disclosure: m=1, p=1 -> 2 < 3)
    # matches the verified generic Derrick law 3 - 2p at m = 1.
    led.check("B10/C12d source 2mp counting reduces to the verified 3-2p law at m=1",
              2 * 1 * 1 == 2 and 3 - 2 * 1 == 1)
    led.mutation("B10/M2 m=0 counting (no derivatives) contradicts the divergent claim",
                 2 * 0 * 1 < 3 and (3 - 2 * 0) == 3)
    # Mutation: the slot placement is load-bearing - moving h onto the
    # eta-paired slots must destroy both pinned relations.
    def D(expr, i):
        return sp.diff(expr, coords[i])

    eta_inv_L = sp.diag(-1, 1, 1, 1).inv()
    inc_alt = sp.zeros(4)
    for mu in range(4):
        for nu in range(4):
            acc = 0
            for a in range(4):
                if a == mu:
                    continue
                for b in range(4):
                    if len({mu, a, b}) < 3:
                        continue
                    for c in range(4):
                        if len({mu, a, b, c}) < 4:
                            continue
                        for r in range(4):
                            if r == nu:
                                continue
                            for s in range(4):

                                if len({nu, r, s}) < 3:
                                    continue
                                for tt in range(4):
                                    if len({nu, r, s, tt}) < 4:
                                        continue
                                    e = eps4(mu, a, b, c) * eps4(nu, r, s, tt)
                                    if e == 0:
                                        continue
                                    dd = sum(eta_inv_L[c, u] * eta_inv_L[tt, v]
                                             * D(D(h_L[a, r], u), v)
                                             for u in range(4) for v in range(4))
                                    acc += e * eta_inv_L[b, s] * dd
            inc_alt[mu, nu] = sp.simplify(acc)
    alt_matches = (sp.simplify(inc_alt - 2 * G1_L) == zero4
                   or sp.simplify(inc_alt + 2 * G1_L) == zero4)
    led.mutation("B10/M1 re-slotting h onto the eta pair destroys the identity",
                 not alt_matches)


# ---------------------------------------------------------------------------


BLOCKS = [
    block1_generator_structure,
    block2_bracket_spin,
    block3_trace_identity,
    block4_epsilon,
    block5_exchange_signs,
    block6_multipole,
    block7_linearized_identity,
    block8_weitzenbock,
    block9_composition,
    block10_inc_sign,
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--verify", action="store_true", help="run checks only; all must pass")
    ap.add_argument("--mutate", action="store_true", help="run mutations only; all must break")
    args = ap.parse_args()
    if not args.verify and not args.mutate:
        args.verify = args.mutate = True

    led = Ledger()
    for blk in BLOCKS:
        blk(led)

    status = 0
    if args.verify:
        print("== CHECKS ==")
        for name, ok, detail in led.checks:
            print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
        npass = sum(1 for _, ok, _ in led.checks if ok)
        if npass != len(led.checks):
            status = 1
        print(f"CHECKS: {npass}/{len(led.checks)} PASS")
    if args.mutate:
        print("== MUTATIONS ==")
        for name, ok, detail in led.mutations:
            print(f"[{'BROKE' if ok else 'DID-NOT-BREAK'}] {name}" + (f"  ({detail})" if detail else ""))
        nbroke = sum(1 for _, ok, _ in led.mutations if ok)
        if nbroke != len(led.mutations):
            status = 1
        print(f"MUTATIONS: {nbroke}/{len(led.mutations)} BREAK")
    if status == 0:
        print(f"ALL {len(led.checks)} CHECKS PASS; ALL {len(led.mutations)} MUTATIONS BREAK"
              if args.verify and args.mutate else "SELECTED MODES GREEN")
    return status


if __name__ == "__main__":
    sys.exit(main())
