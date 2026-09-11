#!/usr/bin/env python3
"""T2 back-reaction compute (beacon, 0157/tilt, chartered post-freeze).
Frozen design: tilt/t2-design.md (committed BEFORE this compute).
SymPy: prestress tensor + PK structure (F1 guard: F_glide IDENTICALLY
0); numeric: fiducial evaluation (F2 core-shift vs B1 7e-4 bars;
F3 Peierls barrier sanity). Sketch coefficients in the freeze are
superseded by receipted exact forms where they differ (frozen rule).
Fiducial (DECLARED): cell ell^3/ring, Lam=2pi.R/ell^3, R/ell=0.1,
ell/a=100, th=0.1. Verdict = parametric inequality + fiducial number.
Exit nonzero on F1 breach or arithmetic contradiction.
"""

import sympy as sp

R_OVER_L = 0.1
L_OVER_A = 100.0
TH = 0.1
B1_BAR = 7e-4


def main() -> None:
    J0, ell, th, b = sp.symbols("J0 ell th b", positive=True)
    # mu enters only via parametric beta below (banked mu_aff form cited
    # in freeze); no separate symbol needed.
    # Prestress from w_c = (J0/ell^3)(4.eta_in-3.eps_zz)th^2/2:
    s = J0 * th ** 2 / (2 * ell ** 3)
    sig = sp.diag(2 * s, 2 * s, -3 * s)
    bv = sp.Matrix([b, 0, 0])
    t = sp.Matrix([0, 0, 1])
    F = (sig * bv).cross(t)
    print(f"PK force (edge b=xhat): F = {list(F)}")
    assert sp.simplify(F[0]) == 0, "F1 FIRES: glide nonzero (bug)"
    assert sp.simplify(F[2]) == 0
    Fy = sp.simplify(F[1])
    print(f"F1 guard: F_glide IDENTICALLY 0 at O(th^2); "
          f"F_climb = {Fy} (receipted exact form)")
    # Relative core shift: s_fac.sigma_xx.b^2 / (mu.b^2), s_fac=1 fiducial
    beta = 10 * sp.pi / sp.log(L_OVER_A) * R_OVER_L ** 3
    beta_n = float(beta.evalf())
    shift = beta_n * TH ** 2  # s_fac = 1
    print(f"beta = {beta_n:.3e} (parametric 10pi/ln(ell/a).(R/ell)^3); "
          f"rel core shift (s=1) = {shift:.3e} vs B1 bar {B1_BAR:.1e}")
    assert shift > 0
    if shift > B1_BAR:
        raise SystemExit("F2 FIRES FOR EDGE: tilt unbinds core "
                         f"({shift:.2e} > {B1_BAR:.1e})")
    print(f"F2 verdict: BOUND (headroom {B1_BAR / shift:.1f}x at s=1; "
          f"unbind iff beta.th^2.s > {B1_BAR:.1e})")
    # Peierls (tilt wave, long-wave limit, M=1/2 fiducial core moment):
    bar = 0.5 * beta_n * TH ** 2
    print(f"F3: Peierls barrier/core ~ {bar:.3e} (0<bar<<1 expected; "
          f"long-wave limit, averaging suppression unpriced)")
    assert 0 < bar < 1e-2, "F3 anomaly: barrier out of estimate band"
    print("T2 verdict: glide unaffected O(th^2); climb activated (mass "
          "transport unpriced); core BOUND; barrier negligible -> mobile.")


if __name__ == "__main__":
    main()
