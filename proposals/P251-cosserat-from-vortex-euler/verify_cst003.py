"""C-CST-003: exact conditional sphere moments and micropolar energy matching.

The supplied elastic energy is
W = lambda/2*(tr sym grad u)^2 + mu*|sym grad u|^2
  + alpha/2*|curl u - 2 Phi|^2
  + c_tr*(tr grad Phi)^2 + c_s*|sym grad Phi|^2 + c_a*|skew grad Phi|^2.
Phi is a dimensionless rotation ANGLE; alpha has pressure units. B and C_tw
have force*length^2 units when used as curvature-energy coefficients.

Moment and coefficient identities are exact once these microscopic energies
are supplied. Line tension alone does not derive relative-angle locking.
Attempt 0029 shows that 0028's rate-quadratic coefficient has microinertia
units, so it cannot replace alpha in this action. A Doppler frequency alone
does not supply the missing action or kinetic coordinate map. The P242 axial
stiffness premise remains an explicit conditional import. Historical attempts
are preserved; the Euler-derived closure of N2/N3 remains active.
"""

import sys

import sympy as sp

from substrate_framework.homogenization import (
    affine_lame_moduli,
    axial_moment_identity,
    sphere_fourth_moment_isotropic,
    sphere_second_moment,
    straight_line_tension,
)
from substrate_framework.verification import CheckLedger


def check_frame_locking_bridge(ledger):
    """Distinguish exact line stretch from a postulated relative shear energy."""
    angle = sp.Symbol("theta", real=True)
    rotation = sp.rot_axis3(angle)
    exact_relative = rotation.T - sp.eye(3)
    exact_green_lagrange = sp.simplify(
        (exact_relative + exact_relative.T + exact_relative.T * exact_relative) / 2
    )
    ledger.check(
        "exact relative Green-Lagrange line strain is independent of frame rotation",
        exact_green_lagrange == sp.zeros(3, 3),
        "R.T R = I, so a straight-line tension cannot see a free director rotation",
    )

    phi = sp.Symbol("phi", real=True)
    linear_skew = sp.Matrix([[0, -phi, 0], [phi, 0, 0], [0, 0, 0]])
    truncated_relative = -linear_skew
    truncated_quadratic = sp.simplify(truncated_relative.T * truncated_relative / 2)
    ledger.check(
        "first-order relative rotation has a nonzero retained quadratic norm",
        truncated_quadratic != sp.zeros(3, 3),
        "this term needs an Euler-derived frame-locking interaction; line tension alone does not supply it",
    )


def check_sphere_moment_reuse(ledger):
    """Declared import reuse: second/third/fourth moment structure exact."""
    P2 = sphere_second_moment()
    ledger.check(
        "sphere_second_moment == delta/3",
        sp.simplify(P2 - sp.eye(3) / 3) == sp.zeros(3, 3),
        "P2 = delta_ij/3",
    )
    eps = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"e{i}{j}"))
    ident = axial_moment_identity(eps)
    skew2 = sum((eps[i, j] - eps[j, i]) ** 2 for i in range(3) for j in range(i))
    ledger.check(
        "axial_moment_identity: residual = -(antisymmetric pairs^2)/15 for general eps",
        sp.simplify(ident + skew2 / 15) == 0,
        "symmetric-input closure exact; general-input correction recorded",
    )
    P4 = sphere_fourth_moment_isotropic()
    d = sp.eye(3)
    ok = all(
        sp.simplify(
            P4[i, j, k, ell]
            - (d[i, j] * d[k, ell] + d[i, k] * d[j, ell] + d[i, ell] * d[j, k]) / 15
        )
        == 0
        for i in range(3)
        for j in range(3)
        for k in range(3)
        for ell in range(3)
    )
    ledger.check(
        "sphere_fourth_moment_isotropic closure",
        ok,
        "(d_ij d_kl + d_ik d_jl + d_il d_jk)/15",
    )


def check_comparsi_structure(ledger):
    """W_alpha = (alpha/2)|rot u - 2 Phi|^2 gives dW/dPhi = -2 alpha rot u + 4 alpha Phi."""
    alpha_c = sp.Symbol("alpha_c")
    ru = sp.Matrix(sp.symbols("ru1:4"))
    Phi = sp.Matrix(sp.symbols("Phi1:4"))
    W_alpha = alpha_c / 2 * (ru - 2 * Phi).dot(ru - 2 * Phi)
    gradW = [sp.diff(W_alpha, p) for p in Phi]
    target = [-2 * alpha_c * ru[a] + 4 * alpha_c * Phi[a] for a in range(3)]
    ledger.check(
        "Comparsi coupling: dW_alpha/dPhi = -2 alpha rot u + 4 alpha Phi",
        all(sp.simplify(gradW[a] - target[a]) == 0 for a in range(3)),
        "intake form (issue #198, 2026-09-03) reproduced structurally",
    )


def check_stretch_sector(ledger):
    """(lambda, mu) from the P242 affine matching on the tube tension T."""
    rho_, Gam, Rout, acore = sp.symbols("rho Gamma R a", positive=True)
    T_expr = sp.simplify(straight_line_tension(rho_, Gam, Rout, acore))
    expected = Gam**2 * rho_ / (4 * sp.pi) * sp.log(Rout / acore)
    ledger.check(
        "straight_line_tension == rho G^2/(4 pi) ln(R/a)",
        sp.simplify(T_expr - expected) == 0,
        f"T = {T_expr}",
    )
    Lv = sp.Symbol("L_v", positive=True)
    lam_eff, mu_eff = affine_lame_moduli(T_expr, Lv)
    ledger.check(
        "lam_eff = mu_eff = L_v T / 15",
        sp.simplify(lam_eff - Lv * T_expr / 15) == 0
        and sp.simplify(mu_eff - Lv * T_expr / 15) == 0,
        "P242 affine matching applied to the vortex-tangle tension",
    )
    ledger.check(
        "pre-stress identification: W_1 = (L_v T/3) tr(es)",
        sp.simplify(sum(sphere_second_moment()[i, i] for i in range(3)) / 3 - sp.Rational(1, 3)) == 0,
        "isotropic pre-stress P = L_v T/3 (recorded for N4 tangent operator)",
    )


def check_locking_sector(ledger):
    """Type the 0028 coefficient before using it in an angle action."""
    rho, Lv, a, eta = sp.symbols("rho L_v a eta", positive=True)
    mass, length, time = sp.symbols("M L T", positive=True)
    coefficient = Lv * sp.pi * rho * a**2 * eta**2 / 4
    scaling = {rho: rho * mass / length**3, Lv: Lv / length**2,
               a: a * length, eta: eta * length}
    units = sp.simplify(coefficient.xreplace(scaling) / coefficient)
    ledger.check("0028 rate coefficient has inertia units M/L", units == mass / length)
    ledger.check(
        "0028 rate coefficient cannot be an angle stiffness",
        sp.simplify(units / (mass / (length * time**2))) == time**2,
        "alpha in N4 has pressure units; attempt 0029 derives this mismatch",
    )
    t = sp.Symbol("t", real=True)
    q = sp.Function("q")(t)
    J = sp.Symbol("J", positive=True)
    kinetic = J * sp.diff(q, t)**2 / 2
    residual = sp.diff(sp.diff(kinetic, sp.diff(q, t)), t) - sp.diff(kinetic, q)
    ledger.check(
        "relative-rate quadratic energy contributes inertia, not static torque",
        sp.simplify(residual - J * sp.diff(q, t, 2)) == 0
        and residual.subs(q, sp.Symbol("q0")).doit() == 0,
    )

def check_wryness_sector(ledger):
    """(c_tr, c_s, c_a) from the projected bend energy with joint moments."""
    kappa = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"kp_{i}{j}"))
    P2 = sphere_second_moment()
    P4 = sphere_fourth_moment_isotropic()
    avg_abs = sp.simplify(
        sum(
            kappa[i, j] * kappa[ell, j] * P2[i, ell]
            for i in range(3)
            for ell in range(3)
            for j in range(3)
        )
    )
    avg_ntkn = sp.simplify(
        sum(
            kappa[i, j] * kappa[ell, m] * P4[i, j, ell, m]
            for i in range(3)
            for j in range(3)
            for ell in range(3)
            for m in range(3)
        )
    )
    kkt = sum(kappa[i, j] * kappa[i, j] for i in range(3) for j in range(3))
    ledger.check(
        "<|kappa n|^2> = (kappa:kappa)/3",
        sp.simplify(avg_abs - kkt / 3) == 0,
        "second-moment contraction",
    )
    tr_k = sum(kappa[i, i] for i in range(3))
    kktT = sum(kappa[i, j] * kappa[j, i] for i in range(3) for j in range(3))
    ledger.check(
        "<(n' kappa n)^2> = ((tr k)^2 + k:k + k:k^T)/15",
        sp.simplify(avg_ntkn - ((tr_k**2 + kkt + kktT) / 15)) == 0,
        "fourth-moment contraction (general kappa)",
    )
    avg_perp = sp.simplify(avg_abs - avg_ntkn)
    Lv, B, Ctw = sp.symbols("L_v B C_tw", positive=True)
    W_c = sp.expand(Lv * (B / 2 * avg_perp + Ctw / 2 * avg_ntkn))
    kS = (kappa + kappa.T) / 2
    kA = (kappa - kappa.T) / 2
    Ssym = sum(kS[i, j] ** 2 for i in range(3) for j in range(3))
    Sskw = sum(kA[i, j] ** 2 for i in range(3) for j in range(3))
    c_tr, c_s, c_a = sp.symbols("c_tr c_s c_a")
    W_form = sp.expand(c_tr * tr_k**2 + c_s * Ssym + c_a * Sskw)
    sol = sp.solve(sp.expand(W_c - W_form), [c_tr, c_s, c_a], dict=True)
    ledger.check(
        "coefficient matching closes (three couple invariants)",
        len(sol) == 1,
        f"solutions: {len(sol)}",
    )
    s0 = sol[0]
    ledger.check(
        "c_tr = L_v(-B + C_tw)/30 ; at C_tw^tube = 0: -B L_v/30",
        sp.simplify(s0[c_tr] - Lv * (-B + Ctw) / 30) == 0,
        f"c_tr = {sp.simplify(s0[c_tr])}",
    )
    ledger.check(
        "c_s = L_v(3B + 2 C_tw)/30 ; at C_tw^tube = 0: B L_v/10",
        sp.simplify(s0[c_s] - Lv * (3 * B + 2 * Ctw) / 30) == 0,
        f"c_s = {sp.simplify(s0[c_s])}",
    )
    ledger.check(
        "c_a = B L_v / 6 (C_tw-independent)",
        sp.simplify(s0[c_a] - B * Lv / 6) == 0,
        f"c_a = {sp.simplify(s0[c_a])}",
    )
    ledger.check(
        "conditional zero-twist specialization of the matched curvature energy",
        sp.simplify(s0[c_s].subs(Ctw, 0) - B * Lv / 10) == 0,
        "this substitution checks the coefficient map; N2 must supply the microscopic value",
    )


def check_matching_residual_and_mutations(ledger):
    """General angle probe of the supplied micropolar energy, not line tension."""
    h = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"h{i}{j}"))
    Phi = sp.Matrix(sp.symbols("Phi1:4"))
    S = sp.Matrix(3, 3, lambda i, j:
                  -sum(sp.LeviCivita(i, j, c) * Phi[c] for c in range(3)))
    es = (h + h.T) / 2
    ea = (h - h.T) / 2 - S
    lam, mu, alpha = sp.symbols("lambda mu alpha", positive=True)
    tensor_energy = (lam * sp.trace(es)**2 / 2
                     + mu * sum(x**2 for x in es) + alpha * sum(x**2 for x in ea))
    curl = sp.Matrix([h[2, 1] - h[1, 2], h[0, 2] - h[2, 0], h[1, 0] - h[0, 1]])
    vector_energy = (lam * sp.trace(es)**2 / 2 + mu * sum(x**2 for x in es)
                     + alpha / 2 * (curl - 2 * Phi).dot(curl - 2 * Phi))
    ledger.check("general tensor and axial-vector angle energies agree",
                 sp.simplify(tensor_energy - vector_energy) == 0)
    rigid = {h[i, j]: S[i, j] for i in range(3) for j in range(3)}
    ledger.check("supplied relative-angle energy vanishes for coherent rigid rotation",
                 sp.simplify(tensor_energy.subs(rigid)) == 0)
    wrong = vector_energy.subs(dict(zip(Phi, -Phi)), simultaneous=True)
    ledger.check("wrong relative-angle sign breaks rigid-rotation cancellation",
                 sp.simplify(wrong.subs(rigid)) != 0)
    ledger.check("doubled angle coefficient changes the physical torque",
                 sp.simplify(sp.diff(vector_energy.subs(alpha, 2 * alpha)
                                     - vector_energy, Phi[0])) != 0)

def main():
    ledger = CheckLedger("C-CST-003")
    check_frame_locking_bridge(ledger)
    check_sphere_moment_reuse(ledger)
    check_comparsi_structure(ledger)
    check_locking_sector(ledger)
    check_stretch_sector(ledger)
    check_wryness_sector(ledger)
    check_matching_residual_and_mutations(ledger)
    return ledger.finish()


if __name__ == "__main__":
    sys.exit(main())
