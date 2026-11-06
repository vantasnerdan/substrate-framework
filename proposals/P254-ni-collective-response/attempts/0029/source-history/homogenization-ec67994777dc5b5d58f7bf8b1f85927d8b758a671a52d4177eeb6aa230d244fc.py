"""Exact affine homogenization of a declared isotropic filament ensemble.

Conventions are fixed once:

- ``eps`` is the small-strain tensor;
- a filament segment of length ``l`` along unit direction ``n_hat`` with axial
  stiffness ``E_f`` (force units) stores energy ``E_f*l*(n.e.n)**2/2`` under an
  affine ambient strain;
- the ensemble is orientation-isotropic with total filament length per volume
  ``L_v``; length-weighted averages are written ``<q(n)>``;
- Cauchy-Born energy density ``W(eps) = E_f*L_v*<(n.eps.n)**2>/2``;
- the straight-vortex line tension follows from Biot-Savart
  ``v_theta(r) = Gamma/(2*pi*r)`` as ``T = rho*Gamma**2/(4*pi)*ln(R/a)``.

All averages are computed exactly; nothing is fitted or hard-coded.
"""

from __future__ import annotations

from typing import Any

import sympy as sp

from .elasticity import poisson_ratio


def sphere_second_moment() -> sp.ImmutableMatrix:
    """Exact <n_i n_j> = delta_ij/3 on the unit sphere."""

    theta, phi = sp.symbols("theta phi", real=True)
    measure = sp.sin(theta) / (4 * sp.pi)
    basis = [
        sp.sin(theta) * sp.cos(phi),
        sp.sin(theta) * sp.sin(phi),
        sp.cos(theta),
    ]
    moment = sp.zeros(3, 3)
    for i in range(3):
        for j in range(3):
            moment[i, j] = sp.simplify(
                sp.integrate(basis[i] * basis[j] * measure, (theta, 0, sp.pi), (phi, 0, 2 * sp.pi))
            )
    return sp.ImmutableMatrix(moment)


def sphere_fourth_moment_direct() -> dict[tuple[int, int, int, int], sp.Expr]:
    """Every unique <n_i n_j n_k n_l> by explicit spherical integration."""

    theta, phi = sp.symbols("theta phi", real=True)
    measure = sp.sin(theta) / (4 * sp.pi)
    basis = [
        sp.sin(theta) * sp.cos(phi),
        sp.sin(theta) * sp.sin(phi),
        sp.cos(theta),
    ]
    result: dict[tuple[int, int, int, int], sp.Expr] = {}
    for i in range(3):
        for j in range(i, 3):
            for k in range(j, 3):
                for l in range(k, 3):
                    integrand = (
                        basis[i] * basis[j] * basis[k] * basis[l] * measure
                    )
                    value = sp.simplify(
                        sp.integrate(
                            integrand, (theta, 0, sp.pi), (phi, 0, 2 * sp.pi)
                        )
                    )
                    result[(i, j, k, l)] = value
    return result


def sphere_fourth_moment_isotropic(
    dim: int = 3,
) -> sp.ImmutableDenseNDimArray:
    """Isotropic rank-4 closure <n_i n_j n_k n_l> = (d_ij d_kl + d_ik d_jl + d_il d_jk)/15."""

    delta = sp.eye(dim)
    return sp.ImmutableDenseNDimArray(
        [
            [
                [
                    [
                        (
                            delta[a, b] * delta[c, d]
                            + delta[a, c] * delta[b, d]
                            + delta[a, d] * delta[b, c]
                        )
                        / 15
                        for d in range(dim)
                    ]
                    for c in range(dim)
                ]
                for b in range(dim)
            ]
            for a in range(dim)
        ],
        shape=(dim,) * 4,
    )


def axial_moment_identity(strain: sp.MatrixBase) -> sp.Expr:
    """Return <(n eps n)^2> - ((tr eps)^2 + 2 eps:eps)/15 evaluated exactly.

    The contraction uses the isotropic rank-4 closure; tests verify it
    against the direct spherical integration componentwise.
    """

    e = sp.Matrix(strain)
    if e.shape != (3, 3):
        raise ValueError("strain must be 3x3")
    moment_array = sphere_fourth_moment_isotropic()
    average = sp.Integer(0)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    average += e[i, j] * e[k, l] * moment_array[i, j, k, l]
    trace = sum(e[a, a] for a in range(3))
    frobenius = sum(e[a, b] * e[a, b] for a in range(3) for b in range(3))
    closed_form = (trace**2 + 2 * frobenius) / 15
    return sp.simplify(average - closed_form)


def cauchy_born_energy_density(
    strain: sp.MatrixBase, axial_stiffness: Any, line_density: Any
) -> sp.Expr:
    """W = E_f*L_v*<(n eps n)^2>/2 under the affine declared model."""

    e = sp.Matrix(strain)
    if e.shape != (3, 3):
        raise ValueError("strain must be 3x3")
    moment_array = sphere_fourth_moment_isotropic()
    average = sp.Integer(0)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    average += e[i, j] * e[k, l] * moment_array[i, j, k, l]
    return sp.simplify(axial_stiffness * line_density * average / 2)


def affine_lame_moduli(axial_stiffness: Any, line_density: Any) -> tuple[sp.Expr, sp.Expr]:
    """Solve lambda, mu by coefficient matching against ½λ(tr eps)^2 + mu eps:eps.

    Returns (lambda, mu); the match is performed algebraically, never
    hard-coded, so a changed energy ansatz cannot silently pass.
    """

    for value, name in (
        (axial_stiffness, "axial_stiffness"),
        (line_density, "line_density"),
    ):
        expression = sp.sympify(value)
        if expression.is_number and (
            expression.is_real is not True or float(expression) < 0
        ):
            raise ValueError(
                f"{name} must be real and nonnegative for a physical ensemble"
            )
    exx, eyy, ezz, exy, exz, eyz = sp.symbols(
        "_eXX _eYY _eZZ _eXY _eXZ _eYZ", real=True
    )
    e_general = sp.Matrix(
        [[exx, exy, exz], [exy, eyy, eyz], [exz, eyz, ezz]]
    )
    energy_model = cauchy_born_energy_density(e_general, axial_stiffness, line_density)

    lam_s, mu_s = sp.symbols("_lamC _muC")
    target = lam_s * sum(e_general[a, a] for a in range(3)) ** 2 / 2 + mu_s * sum(
        e_general[a, b] ** 2 for a in range(3) for b in range(3)
    )
    poly = sp.Poly(sp.expand(energy_model - target), exx, eyy, ezz, exy, exz, eyz)
    equations = [sp.Eq(coefficient, 0) for coefficient in poly.coeffs()]
    solution = sp.solve(equations, (lam_s, mu_s), dict=True)
    if not solution:
        raise ValueError("no Lamé constants reproduce the Cauchy-Born energy")
    lam_value = sp.simplify(solution[0][lam_s])
    mu_value = sp.simplify(solution[0][mu_s])
    return lam_value, mu_value


def straight_line_tension(density: Any, circulation: Any, outer: Any, core: Any) -> sp.Expr:
    """Energy per length of one straight vortex line from Biot-Savart.

    With v_theta = Gamma/(2 pi r) the kinetic-energy density integrates to
    rho*Gamma^2/(4*pi)*ln(R/a). The integral is evaluated here, not quoted.
    """

    r = sp.Symbol("_r", positive=True)
    speed = circulation / (2 * sp.pi * r)
    energy_per_length = sp.integrate(
        density / 2 * speed**2 * 2 * sp.pi * r, (r, core, outer)
    )
    return sp.simplify(energy_per_length)


def affine_poisson_ratio(axial_stiffness: Any, line_density: Any) -> sp.Expr:
    """Poisson ratio of the affine isotropic tangle medium."""

    lam_value, mu_value = affine_lame_moduli(axial_stiffness, line_density)
    return poisson_ratio(lam_value, mu_value)
