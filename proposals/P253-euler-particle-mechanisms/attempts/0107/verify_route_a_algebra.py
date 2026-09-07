#!/usr/bin/env python3
"""Exact exposing checks for P253/0107 Route-A algebra.

This verifier checks only local vector identities, normalizations, and the
finite-dimensional determinant.  It does not prove cutoff convergence,
witness nondegeneracy, nonlinear IFT hypotheses, or a P2 statement.
"""

from __future__ import annotations

import sympy as sp


def curl(v: sp.Matrix, xyz: tuple[sp.Symbol, ...]) -> sp.Matrix:
    x, y, z = xyz
    return sp.Matrix(
        [
            sp.diff(v[2], y) - sp.diff(v[1], z),
            sp.diff(v[0], z) - sp.diff(v[2], x),
            sp.diff(v[1], x) - sp.diff(v[0], y),
        ]
    )


def divergence(v: sp.Matrix, xyz: tuple[sp.Symbol, ...]) -> sp.Expr:
    return sum(sp.diff(v[i], xyz[i]) for i in range(3))


def divergence_tensor(tensor: sp.Matrix, xyz: tuple[sp.Symbol, ...]) -> sp.Matrix:
    return sp.Matrix(
        [sum(sp.diff(tensor[i, j], xyz[j]) for j in range(3)) for i in range(3)]
    )


def check(label: str, predicate: bool) -> None:
    if not predicate:
        raise AssertionError(label)
    print(f"PASS {label}")


def maxwell_residual(
    electric: sp.Matrix,
    magnetic: sp.Matrix,
    xyz: tuple[sp.Symbol, ...],
    time: sp.Symbol,
    epsilon: sp.Symbol,
    permeability: sp.Symbol,
) -> sp.Matrix:
    electric_t = electric.diff(time)
    magnetic_t = -curl(electric, xyz)
    charge = epsilon * divergence(electric, xyz)
    current = curl(magnetic, xyz) / permeability - epsilon * electric_t
    identity = sp.eye(3)
    sigma = epsilon * (electric * electric.T - electric.dot(electric) * identity / 2)
    sigma += (
        magnetic * magnetic.T - magnetic.dot(magnetic) * identity / 2
    ) / permeability
    lhs = epsilon * (electric_t.cross(magnetic) + electric.cross(magnetic_t))
    force = charge * electric + current.cross(magnetic)
    return (lhs - (divergence_tensor(sigma, xyz) - force)).applyfunc(sp.simplify)


def main() -> None:
    x, y, z, t = sp.symbols("x y z t", real=True)
    xyz = (x, y, z)
    eps, mu = sp.symbols("epsilon_EM mu_EM", positive=True)
    e_components = [sp.Function(name)(x, y, z, t) for name in ("E1", "E2", "E3")]
    b_components = [sp.Function(name)(x, y, z, t) for name in ("B1", "B2", "B3")]
    electric = sp.Matrix(e_components)
    magnetic = sp.Matrix(b_components)
    residual = maxwell_residual(electric, magnetic, xyz, t, eps, mu)
    monopole_defect = (residual + magnetic * divergence(magnetic, xyz) / mu).applyfunc(
        sp.simplify
    )
    check(
        "unconstrained Maxwell identity exposes -(B/mu) div B",
        monopole_defect == sp.zeros(3, 1),
    )

    potential = sp.Matrix(
        [sp.Function(name)(x, y, z, t) for name in ("A1", "A2", "A3")]
    )
    magnetic_curl = curl(potential, xyz)
    check("div curl A enforces magnetic Gauss row", sp.simplify(divergence(magnetic_curl, xyz)) == 0)
    constrained_residual = maxwell_residual(
        electric, magnetic_curl, xyz, t, eps, mu
    )
    check(
        "constrained Maxwell momentum sign: dt(epsilon E cross B)=div sigma-f",
        constrained_residual == sp.zeros(3, 1),
    )

    # Integrated identity (1/2) int x cross curl A = int A.  Integration by
    # parts contributes -delta_(j l); this contraction must be +delta_(i m).
    coeff = sp.MutableDenseMatrix.zeros(3, 3)
    for i in range(3):
        for m in range(3):
            coeff[i, m] = sp.Rational(-1, 2) * sum(
                sp.LeviCivita(i, j, k) * sp.LeviCivita(k, j, m)
                for j in range(3)
                for k in range(3)
            )
    check("impulse curl factor and sign", coeff == sp.eye(3))

    rho = sp.symbols("rho_m", positive=True)
    a, b, c = sp.symbols("a b c")
    moment = sp.Matrix([[0, a, b], [-a, 0, c], [-b, -c, 0]])
    impulse = sp.Matrix(
        [
            rho
            * sp.Rational(1, 2)
            * sum(
                sp.LeviCivita(l, m, n) * moment[m, n]
                for m in range(3)
                for n in range(3)
            )
            for l in range(3)
        ]
    )
    reconstructed = sp.Matrix(
        3,
        3,
        lambda m, n: sum(sp.LeviCivita(l, m, n) * impulse[l] for l in range(3)) / rho,
    )
    check("moment tensor reconstructed from hydrodynamic impulse", reconstructed == moment)

    for axis in range(3):
        basis = sp.eye(3)[:, axis]
        position = sp.Matrix([x, y, z])
        potential = basis.cross(position) / 2
        check(f"compact-center core potential axis {axis}", curl(potential, xyz) == basis)

    h, c_row, g, ki, kc = sp.symbols("H C G K_I K_C")
    lower = sp.Matrix([[h, 0, 0], [0, c_row, 0], [ki, kc, g]])
    check("lower-triangular row determinant", sp.factor(lower.det()) == h * c_row * g)

    omega = sp.symbols("Omega", positive=True)
    zeta = 2 * omega
    k_center = sp.simplify(8 * zeta / omega)
    check("uncharged smooth-center K equals 16", k_center == 16)

    print("ALL 10 ROUTE-A ALGEBRA CHECKS PASSED")


if __name__ == "__main__":
    main()
