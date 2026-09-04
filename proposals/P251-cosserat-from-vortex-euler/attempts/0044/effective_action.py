"""Derive the local action before orientational averaging or boundary rewrites."""

from __future__ import annotations

import sympy as sp

from bloch_sector import matrices, taylor_matrix
from substrate_framework.homogenization import sphere_fourth_moment_isotropic
from substrate_framework.verification import CheckLedger


def main() -> int:
    ledger = CheckLedger("P251-0044-isotropic-action-repair")
    kappa = sp.symbols("kappa", real=True)
    curvatures = []
    inertia_terms = []
    i0 = None
    k0 = None
    for axis in range(3):
        _, _, hessian, kks = matrices(axis, kappa)
        h0 = taylor_matrix(hessian, kappa, 0).applyfunc(sp.simplify)
        h2 = taylor_matrix(hessian, kappa, 2).applyfunc(sp.simplify)
        o0 = taylor_matrix(kks, kappa, 0).applyfunc(sp.simplify)
        o2 = taylor_matrix(kks, kappa, 2).applyfunc(sp.simplify)
        ledger.check(
            f"axis {axis}: no chiral first energy derivative",
            taylor_matrix(hessian, kappa, 1).applyfunc(sp.simplify) == sp.zeros(2),
        )
        ledger.check(
            f"axis {axis}: no chiral first KKS derivative",
            taylor_matrix(kks, kappa, 1).applyfunc(sp.simplify) == sp.zeros(2),
        )
        inertia = sp.simplify(o0[0, 1] ** 2 / h0[1, 1])
        i2 = sp.simplify(inertia * (2 * o2[0, 1] / o0[0, 1] - h2[1, 1] / h0[1, 1]))
        effective_c = sp.simplify(h2[0, 0] - h0[0, 0] * i2 / inertia)
        if i0 is not None:
            ledger.check(f"axis {axis}: common zero-order inertia", inertia == i0)
        i0, k0 = inertia, h0[0, 0]
        inertia_terms.append(i2)
        curvatures.append(effective_c)
        print(
            f"AXIS {'xyz'[axis]} raw I2={i2}; raw K2={h2[0, 0]}; normalized C={effective_c}"
        )
    ledger.check(
        "single-cell axial curvature is genuinely negative", curvatures[2].is_negative
    )
    ct = sp.factor((curvatures[0] + curvatures[1]) / 2)
    cz = curvatures[2]
    acoef = sp.factor((4 * ct + cz) / 15)
    bcoef = sp.factor((cz - ct) / 15)
    gradient = sp.Matrix(3, 3, sp.symbols("G0:9", real=True))
    moment = sphere_fourth_moment_isotropic()
    direct_average = sum(
        gradient[i, j]
        * gradient[row2, m]
        * (
            ct * sp.KroneckerDelta(i, row2) * sp.KroneckerDelta(j, m) / 3
            + (cz - ct) * moment[i, row2, j, m]
        )
        for i in range(3)
        for j in range(3)
        for row2 in range(3)
        for m in range(3)
    )
    norm2 = sum(value**2 for value in gradient)
    closed_average = acoef * norm2 + bcoef * (
        sp.trace(gradient) ** 2 + sp.trace(gradient * gradient)
    )
    ledger.check(
        "full action tensor average equals invariant form",
        sp.simplify(direct_average - closed_average) == 0,
    )
    ledger.check(
        "isotropic transverse gradient is strictly positive", acoef.is_positive
    )
    ledger.check(
        "isotropic longitudinal gradient is strictly positive",
        (acoef + 2 * bcoef).is_positive,
    )
    curl2 = sum(
        (gradient[(j + 2) % 3, (j + 1) % 3] - gradient[(j + 1) % 3, (j + 2) % 3]) ** 2
        for j in range(3)
    )
    null_term = sp.trace(gradient * gradient) - sp.trace(gradient) ** 2
    positive_density = acoef * curl2 + (acoef + 2 * bcoef) * sp.trace(gradient) ** 2
    ledger.check(
        "positive curl-div action differs by declared null Lagrangian",
        sp.simplify(closed_average - positive_density - (acoef + bcoef) * null_term)
        == 0,
    )

    coordinates = sp.symbols("x y z", real=True)
    field = [sp.Function(f"Q{j}")(*coordinates) for j in range(3)]
    jacobian = sp.Matrix(3, 3, lambda i, j: sp.diff(field[i], coordinates[j]))
    flux = [
        sum(
            field[i] * sp.diff(field[j], coordinates[i])
            - field[j] * sp.diff(field[i], coordinates[i])
            for i in range(3)
        )
        for j in range(3)
    ]
    divergence_flux = sum(sp.diff(flux[j], coordinates[j]) for j in range(3))
    ledger.check(
        "null Lagrangian is explicit boundary divergence",
        sp.simplify(
            divergence_flux - sp.trace(jacobian * jacobian) + sp.trace(jacobian) ** 2
        )
        == 0,
    )
    ledger.check(
        "not pointwise-positive before boundary rewrite",
        sp.simplify(
            closed_average.subs(dict(zip(list(gradient), list(sp.eye(3)))))
        ).is_negative,
    )
    print(f"I0={i0}, K0={k0}, I2={inertia_terms}, C={curvatures}")
    print(
        f"ISOTROPIC C_perp={acoef}; C_parallel={sp.factor(acoef + 2 * bcoef)}; boundary={sp.factor(acoef + bcoef)}"
    )
    print(f"ISOTROPIC physical mass coefficient = rho ell^2 * {i0 / 3}")
    print(
        f"ISOTROPIC omega_squared transverse k2 = b^2 * {sp.factor(acoef / (i0 / 3))}"
    )
    print(
        f"ISOTROPIC omega_squared longitudinal k2 = b^2 * {sp.factor((acoef + 2 * bcoef) / (i0 / 3))}"
    )
    print(
        "Scope: periodic or compact-support macrofields, explicit isotropic affine angle sector, before the parent's common-angle/body lift."
    )
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
