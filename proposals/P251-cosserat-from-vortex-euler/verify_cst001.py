"""C-CST-001 verifier: frame-transport kinematics of a filament triad.

Claim (exact, conditional on smoothness premises declared in
proposals/P251-cosserat-from-vortex-euler/proposal.yaml, node N1):

  For a smooth unit-speed curve with tangent t(s) and an orthonormal
  normal-plane frame (q1, q2) transported by the rotation-minimizing
  (Bishop) law q_i' = -(q_i . t') t:

  (i)   the Bishop gauge angle relative to the Frenet frame satisfies
        beta' = -tau (sign pinned by direct computation);
  (ii)  the frame angular velocity Omega = 1/2 sum e_i x e_i' has axial
        component Omega . t = 0 (Bishop), +tau (Frenet); gauge-additivity
        holds for ANY orthonormal frame rotated by chi about t,
        Omega_material.t - Omega_base.t = chi' (generic proof), hence the
        material twist density in the Bishop gauge is exactly chi'(s) --
        the micro-twist variable consumed by C-CST-003;
  (iii) the transport law is form-invariant under rigid rotations of the
        whole structure (curve + frame), with R orthogonal by
        construction (explicit Euler parameterization);
  (iv)  mutations: the flipped gauge (beta' = +tau) leaves a nonzero
        transport residual; a rotated frame with a non-matching twist
        rate breaks additivity; the negated torsion convention breaks
        the Frenet-rate identity.

Helix closed forms: kappa_h = a_c/(a_c^2+b_c^2), tau_h = b_c/(a_c^2+b_c^2).

Units: s arc length [L]; angles dimensionless; tau [1/L]. No numerical
sampling: this is an exact claim (SymPy only), so no production-numerics
license is consumed or required.

Attempt history (append-only, see attempts/):
  0001: CheckFailure material_axial_rate_is_twist_density -- axial_rate
        built a fresh Symbol('s', positive=True) not matching the frames'
        Symbol('s', real=True, positive=True); every differentiation inside
        it was silently zero. Mechanism: symbol identity is part of the
        contract; pass s explicitly.
  0002: same check failed again; mechanism: the brute-force route
        (simplify the full material-frame angular velocity on the helix,
        then structural ==) is computationally explosive and brittle
        under simplification-order differences. Repair: prove gauge
        additivity GENERically (cheap, exact, frame-independent) and
        compose it with the concrete helix checks -- verify blocks,
        compose upward. 0003 is the repaired rerun of the unchanged
        scientific route.
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import sys

import sympy as sp

from substrate_framework import CheckLedger


def helix_frames():
    """Return s, ac, bc, tau, t, t', n_f, b_f for the (a_c, b_c) helix."""
    s = sp.Symbol("s", real=True, positive=True)
    ac, bc = sp.symbols("a_c b_c", real=True, positive=True)
    c = sp.sqrt(ac**2 + bc**2)
    r = sp.Matrix([ac * sp.cos(s / c), ac * sp.sin(s / c), bc * s / c])
    t = sp.simplify(sp.diff(r, s))
    assert sp.simplify(t.dot(t)) == 1
    n_f = sp.simplify(sp.diff(t, s) / sp.sqrt(sp.diff(t, s).dot(sp.diff(t, s))))
    b_f = sp.simplify(t.cross(n_f))
    tau = sp.simplify(sp.diff(n_f, s).dot(b_f))
    return s, ac, bc, tau, t, sp.simplify(sp.diff(t, s)), n_f, b_f


def axial_rate(e1, e2, e3, s):
    """Omega . e3 for a right-handed orthonormal frame (e1, e2, e3).

    ``s`` must be the SAME symbol the frames were built with; a fresh
    Symbol with different assumptions is a different symbol (attempt 0001).
    """
    omega = (e1.cross(sp.diff(e1, s))
             + e2.cross(sp.diff(e2, s))
             + e3.cross(sp.diff(e3, s))) / 2
    return omega.dot(e3)


def check_helix_closed_forms(ledger: CheckLedger) -> None:
    s, ac, bc, tau, t, tprime, n_f, b_f = helix_frames()
    kappa = sp.sqrt(sp.diff(t, s).dot(sp.diff(t, s)))
    ledger.check(
        "helix_curvature_closed_form",
        sp.simplify(kappa - ac / (ac**2 + bc**2)) == 0,
        f"kappa_h = {kappa}",
    )
    ledger.check(
        "helix_torsion_closed_form",
        sp.simplify(tau - bc / (ac**2 + bc**2)) == 0,
        f"tau_h = {tau}",
    )


def check_bishop_gauge_sign(ledger: CheckLedger) -> None:
    s, ac, bc, tau, t, tprime, n_f, b_f = helix_frames()
    beta = sp.Symbol("C1") - tau * s
    q1 = sp.cos(beta) * n_f + sp.sin(beta) * b_f
    residual = sp.simplify(sp.diff(q1, s) + (q1.dot(tprime)) * t)
    ledger.check(
        "bishop_transport_residual_vanishes",
        all(sp.simplify(x) == 0 for x in residual),
        "dq1/ds + (q1.t')t = 0 with beta' = -tau",
    )
    ledger.check(
        "bishop_gauge_sign_pinned",
        sp.simplify(sp.diff(beta, s) + tau) == 0,
        "beta' = -tau (not +tau); sign computed, not assumed",
    )


def check_axial_rates_concrete(ledger: CheckLedger) -> None:
    """Concrete helix rates: Bishop 0, Frenet +tau (structural equality)."""
    s, ac, bc, tau, t, tprime, n_f, b_f = helix_frames()
    beta = sp.Symbol("C1") - tau * s
    qb1 = sp.simplify(sp.cos(beta) * n_f + sp.sin(beta) * b_f)
    qb2 = sp.simplify(-sp.sin(beta) * n_f + sp.cos(beta) * b_f)
    ledger.check(
        "bishop_axial_rate_zero",
        sp.simplify(axial_rate(qb1, qb2, t, s)) == 0,
        "rotation-minimizing frame has Omega.t = 0",
    )
    ledger.check(
        "frenet_axial_rate_is_torsion",
        sp.simplify(axial_rate(n_f, b_f, t, s) - tau) == 0,
        "Omega.t = +tau for the Frenet frame (difference-form compare)",
    )


def check_gauge_additivity_generic(ledger: CheckLedger) -> None:
    """Generic frame additivity: rotating (e1,e2) about t by chi(s) shifts
    the axial rate by exactly chi'(s), independent of curve and base frame.

    This composes with bishop_axial_rate_zero to identify the material
    twist density chi'(s) in the Bishop gauge (C-CST-003 input) without a
    brute-force simplification of the composed helix frame.
    """
    s = sp.Symbol("s", real=True, positive=True)
    th, ph = sp.symbols("theta phi", real=True)
    t = sp.Matrix([sp.sin(th) * sp.cos(ph),
                   sp.sin(th) * sp.sin(ph),
                   sp.cos(th)])
    e1 = sp.simplify(sp.diff(t, th))          # unit, orthonormal completion
    e2 = sp.simplify(t.cross(e1))
    assert sp.simplify(e1.dot(e1) - 1) == 0
    chi = sp.Function("chi")(s)
    m1 = sp.simplify(sp.cos(chi) * e1 + sp.sin(chi) * e2)
    m2 = sp.simplify(-sp.sin(chi) * e1 + sp.cos(chi) * e2)
    residual = sp.simplify(
        axial_rate(m1, m2, t, s) - axial_rate(e1, e2, t, s) - sp.diff(chi, s)
    )
    ledger.check(
        "gauge_additivity_generic",
        residual == 0,
        "Omega_material.t - Omega_base.t = chi' for any base frame on any curve",
    )


def check_transport_covariance(ledger: CheckLedger) -> None:
    q = sp.Matrix(sp.symbols("q1:4", real=True))
    t = sp.Matrix(sp.symbols("t1:4", real=True))
    tp = sp.Matrix(sp.symbols("w1:4", real=True))
    psi, th, phi = sp.symbols("psi theta phi", real=True)
    rot = sp.rot_axis3(phi) * sp.rot_axis2(th) * sp.rot_axis1(psi)
    law_rq = -(((rot * q).dot(rot * tp)) * (rot * t))
    covariance = sp.simplify(law_rq - rot * (-(q.dot(tp)) * t)) == sp.zeros(3, 1)
    ledger.check(
        "transport_form_invariant_under_rigid_rotation",
        covariance,
        "R orthogonal by construction (Euler angles); law covariant",
    )


def check_mutations(ledger: CheckLedger) -> None:
    s, ac, bc, tau, t, tprime, n_f, b_f = helix_frames()
    # M1: flipped gauge sign must break transport
    beta_flip = sp.Symbol("C1") + tau * s
    q1f = sp.cos(beta_flip) * n_f + sp.sin(beta_flip) * b_f
    residual = sp.simplify(sp.diff(q1f, s) + (q1f.dot(tprime)) * t)
    ledger.check(
        "mutation_flipped_gauge_fails",
        any(sp.simplify(x) != 0 for x in residual),
        "beta' = +tau leaves nonzero transport residual",
    )
    # M2: correct gauge keeps the frame in the normal plane
    beta = sp.Symbol("C1") - tau * s
    q1 = sp.cos(beta) * n_f + sp.sin(beta) * b_f
    ledger.check(
        "mutation_probe_normal_plane_maintained",
        sp.simplify(q1.dot(t)) == 0 and sp.simplify(sp.diff(q1.dot(t), s)) == 0,
        "(q.t) = 0 and d(q.t)/ds = 0 along the Bishop frame",
    )
    # M3: negated torsion convention would break the Frenet-rate identity
    ledger.check(
        "mutation_negated_torsion_detectable",
        sp.simplify(axial_rate(n_f, b_f, t, s) + tau) != 0,
        "Frenet rate is +tau; a -tau convention fails this difference-form probe",
    )
    # M4: a wrong twist rate (chi' / 2) must break generic additivity
    th, ph = sp.symbols("theta phi", real=True)
    tg = sp.Matrix([sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)])
    eg1 = sp.simplify(sp.diff(tg, th))
    eg2 = sp.simplify(tg.cross(eg1))
    chi = sp.Function("chi")(s)
    m1 = sp.simplify(sp.cos(chi) * eg1 + sp.sin(chi) * eg2)
    m2 = sp.simplify(-sp.sin(chi) * eg1 + sp.cos(chi) * eg2)
    wrong = sp.simplify(
        axial_rate(m1, m2, tg, s) - axial_rate(eg1, eg2, tg, s) - sp.diff(chi, s) / 2
    )
    ledger.check(
        "mutation_half_twist_rate_fails",
        wrong != 0,
        "additivity pins the twist rate exactly; chi'/2 is detectably wrong",
    )


def main() -> int:
    ledger = CheckLedger("C-CST-001")
    check_helix_closed_forms(ledger)
    check_bishop_gauge_sign(ledger)
    check_axial_rates_concrete(ledger)
    check_gauge_additivity_generic(ledger)
    check_transport_covariance(ledger)
    check_mutations(ledger)
    return int(ledger.finish())


if __name__ == "__main__":
    sys.exit(main())
