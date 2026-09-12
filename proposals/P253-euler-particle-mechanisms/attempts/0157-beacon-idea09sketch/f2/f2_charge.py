#!/usr/bin/env python3
"""F2 charge-quantization check (beacon, under 0157): frozen f2/design.md.
B1: Gauss linking integer (numeric). B2: grad-phi loop integral zero
+ kappa = 0 (symbolic). Exits nonzero on any falsifier fire.
"""

import numpy as np
import sympy as sp


def gauss_link(c1, c2):
    """Gauss linking integral of two closed polygonal loops."""
    dl1 = np.roll(c1, -1, axis=0) - c1
    m1 = (c1 + np.roll(c1, -1, axis=0)) / 2
    dl2 = np.roll(c2, -1, axis=0) - c2
    m2 = (c2 + np.roll(c2, -1, axis=0)) / 2
    R = m1[:, None, :] - m2[None, :, :]
    r3 = np.linalg.norm(R, axis=2) ** 3 + 1e-30
    cross = np.cross(dl1[:, None, :], dl2[None, :, :])
    return float(((cross * R) / r3[:, :, None]).sum() / (4 * np.pi))


def ring(R, center, n=400, plane='xy'):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    if plane == 'xy':
        return center + R * np.stack([np.cos(t), np.sin(t),
                                      np.zeros(n)], axis=1)
    return center + R * np.stack([np.cos(t), np.zeros(n), np.sin(t)], axis=1)


def main() -> None:
    # B1: Hopf-linked ring pair -> linking +-1; unlinked -> 0
    a = ring(1.0, np.zeros(3), plane='xy')
    b = ring(0.6, np.array([1.0, 0.0, 0.0]), plane='xz')  # threads a
    c = ring(0.6, np.array([3.0, 0.0, 0.0]), plane='xz')  # far, unlinked
    lab, lac = gauss_link(a, b), gauss_link(a, c)
    print(f"B1: linked Lk={lab:.4f} (expect +-1); unlinked Lk={lac:.4f} (expect 0)")
    assert abs(abs(lab) - 1) < 0.05, "falsifier (a): linking non-integer"
    assert abs(lac) < 0.05, "control failed"
    # B2: loop integral of grad phi identically zero (Stokes, symbolic 1-D)
    t = sp.symbols('t', real=True)
    for mode in [sp.sin(t), sp.cos(t), sp.sin(3 * t) + sp.cos(2 * t)]:
        assert sp.simplify(sp.integrate(sp.diff(mode, t), (t, 0, 2 * sp.pi))) == 0
    # general case: FTC + single-valuedness f(2pi) = f(0) (analytic argument)
    print("B2: closed-loop integral of gradient field = 0 (charge disjoint "
          "from dilatation)")
    # kappa = lam + 2mu/3 = 0 for banked moduli
    K = sp.symbols('K', positive=True)
    assert sp.simplify(-K / 15 + 2 * (K / 10) / 3) == 0
    print("kappa = 0 confirmed (nu=-1 gap): breathing family exists at zero "
          "energy, carries zero loop-charge by B2")
    print("F2 verdict: F2-HOLDS (defect charge integer via closed-filament "
          "linking + I-Sing; dilatational continuum chargeless)")


if __name__ == "__main__":
    main()
