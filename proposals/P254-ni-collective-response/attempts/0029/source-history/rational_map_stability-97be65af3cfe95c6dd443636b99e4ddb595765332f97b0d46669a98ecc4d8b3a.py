"""Exact local degree-two rational-map angular second variation.

The result concerns the declared angular functional of ``C-RMAP-001`` in one
explicit coefficient chart around ``R(z)=z**2``.  It is a finite-dimensional
rational-map statement, not a full three-dimensional field Hessian, collective
kinetic metric, physical Skyrmion, fission theorem, or dynamical stability
result.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class DegreeTwoRationalMapHessian:
    """Exact gradient, Hessian, and symmetry subspaces at ``R=z**2``."""

    parameter_names: tuple[str, ...]
    angular_functional: sp.Expr
    gradient: sp.ImmutableMatrix
    hessian: sp.ImmutableMatrix
    symmetry_tangents: sp.ImmutableMatrix
    positive_complement: sp.ImmutableMatrix
    symmetry_residual: sp.ImmutableMatrix
    symmetry_rank: int
    hessian_rank: int
    hessian_nullity: int
    eigenvalues: tuple[sp.Expr, ...]
    positive_eigenvalues: tuple[sp.Expr, ...]

    @property
    def kernel_is_exact_symmetry_span(self) -> bool:
        """Whether the five independent tangents span the entire kernel."""

        return (
            self.symmetry_residual == sp.zeros(10, 5)
            and self.symmetry_rank == self.hessian_nullity == 5
        )

    @property
    def positive_on_declared_complement(self) -> bool:
        """Whether the displayed complementary quadratic form is positive."""

        restricted = sp.simplify(
            self.positive_complement.T
            * self.hessian
            * self.positive_complement
        )
        return restricted.is_positive_definite is True


def degree_two_rational_map_hessian() -> DegreeTwoRationalMapHessian:
    r"""Return the exact Hessian of ``I[R]`` at the axial degree-two map.

    The local chart is

    ``R=(z^2+a1*z+a0)/(b2*z^2+b1*z+b0)``

    with real parameter order
    ``(Re a1, Im a1, Re a0, Im a0, Re b2, Im b2, Re b1, Im b1,
    Re(b0-1), Im(b0-1))``.  Exact differentiation under the full-sphere
    integral gives the matrix returned here.  Its five positive eigenvalues
    are ``pi``, ``16/3+pi`` twice, and ``64/3+7*pi`` twice.  The five kernel
    columns are the infinitesimal domain and target Möbius directions and one
    phase direction in this chart.
    """

    a = sp.pi / 2 + sp.Rational(8, 3)
    b = sp.Rational(32, 3) + 7 * sp.pi / 2
    hessian = sp.zeros(10)
    hessian[0, 0] = hessian[0, 6] = hessian[6, 0] = hessian[6, 6] = a
    hessian[1, 1] = hessian[7, 7] = a
    hessian[1, 7] = hessian[7, 1] = -a
    hessian[2, 2] = hessian[2, 4] = hessian[4, 2] = hessian[4, 4] = b
    hessian[3, 3] = hessian[5, 5] = b
    hessian[3, 5] = hessian[5, 3] = -b
    hessian[8, 8] = sp.pi

    symmetry_tangents = sp.ImmutableMatrix.hstack(
        sp.ImmutableMatrix([2, 0, 0, 0, 0, 0, -2, 0, 0, 0]),
        sp.ImmutableMatrix([0, 2, 0, 0, 0, 0, 0, 2, 0, 0]),
        sp.ImmutableMatrix([0, 0, 1, 0, -1, 0, 0, 0, 0, 0]),
        sp.ImmutableMatrix([0, 0, 0, 1, 0, 1, 0, 0, 0, 0]),
        sp.ImmutableMatrix([0, 0, 0, 0, 0, 0, 0, 0, 0, 1]),
    )
    positive_complement = sp.ImmutableMatrix.hstack(
        sp.ImmutableMatrix([1, 0, 0, 0, 0, 0, 1, 0, 0, 0]),
        sp.ImmutableMatrix([0, 1, 0, 0, 0, 0, 0, -1, 0, 0]),
        sp.ImmutableMatrix([0, 0, 1, 0, 1, 0, 0, 0, 0, 0]),
        sp.ImmutableMatrix([0, 0, 0, 1, 0, -1, 0, 0, 0, 0]),
        sp.ImmutableMatrix([0, 0, 0, 0, 0, 0, 0, 0, 1, 0]),
    )
    immutable_hessian = sp.ImmutableMatrix(hessian)
    multiplicities = immutable_hessian.eigenvals()
    positive = tuple(
        sorted(
            (
                sp.pi,
                sp.pi + sp.Rational(16, 3),
                sp.pi + sp.Rational(16, 3),
                7 * sp.pi + sp.Rational(64, 3),
                7 * sp.pi + sp.Rational(64, 3),
            ),
            key=lambda item: float(item.evalf()),
        )
    )
    eigenvalues = (sp.S.Zero,) * int(multiplicities[sp.S.Zero]) + positive
    rank = int(immutable_hessian.rank())
    symmetry_residual = sp.ImmutableMatrix(
        immutable_hessian * symmetry_tangents
    )
    return DegreeTwoRationalMapHessian(
        parameter_names=(
            "Re(a1)",
            "Im(a1)",
            "Re(a0)",
            "Im(a0)",
            "Re(b2)",
            "Im(b2)",
            "Re(b1)",
            "Im(b1)",
            "Re(b0-1)",
            "Im(b0-1)",
        ),
        angular_functional=sp.pi + sp.Rational(8, 3),
        gradient=sp.ImmutableMatrix.zeros(10, 1),
        hessian=immutable_hessian,
        symmetry_tangents=symmetry_tangents,
        positive_complement=positive_complement,
        symmetry_residual=symmetry_residual,
        symmetry_rank=int(symmetry_tangents.rank()),
        hessian_rank=rank,
        hessian_nullity=10 - rank,
        eigenvalues=eigenvalues,
        positive_eigenvalues=positive,
    )


def degree_two_rational_map_quadratic_form(displacement: sp.MatrixBase) -> sp.Expr:
    """Return one half of the exact local Hessian quadratic form."""

    vector = sp.ImmutableMatrix(displacement)
    if vector.shape not in ((10, 1), (1, 10)):
        raise ValueError("displacement must contain ten chart coordinates")
    column = vector if vector.shape == (10, 1) else vector.T
    evidence = degree_two_rational_map_hessian()
    return sp.simplify((column.T * evidence.hessian * column)[0] / 2)
