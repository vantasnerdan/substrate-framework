"""Exact real-ell=2 moment tensors and temporal-rank diagnostics.

The tensor map applies to a declared scalar-density coefficient.  It does not
infer a self-consistent field deformation, a localized Floquet mode, a
source-to-waveform law, or physical gravity.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray
import sympy as sp

from .tt_angular import (
    TTPolarizationBasis,
    frobenius_inner_product,
    tt_basis_reconstruct,
    tt_polarization_basis,
)


@dataclass(frozen=True)
class RealL2TTReadout:
    """TT projection and coordinates of a real-ell=2 triple-STF tensor.

    The normalized coordinates contract against unit-Frobenius plus/cross
    tensors.  Conventional matrix readouts are smaller by ``sqrt(2)``.
    """

    basis: TTPolarizationBasis
    source_tensor: sp.Matrix
    projected_tensor: sp.Matrix
    normalized_plus_coordinate: sp.Expr
    normalized_cross_coordinate: sp.Expr
    conventional_plus_readout: sp.Expr
    conventional_cross_readout: sp.Expr


def real_l2_triple_stf_tensor(
    p20: Any = 0,
    m2_cosine: Any = 0,
    m2_sine: Any = 0,
    m1_cosine: Any = 0,
    m1_sine: Any = 0,
) -> sp.Matrix:
    """Map real unnormalized ell=2 density coefficients to ``3*I_STF``.

    The angular basis is

    ``P2(n_z), n_x**2-n_y**2, 2*n_x*n_y, 2*n_x*n_z, 2*n_y*n_z``.

    Each argument is the corresponding radial coefficient
    ``H=4*pi*integral(r**4*h(r), r)``.  Explicitly retaining this convention
    prevents normalized spherical-harmonic factors or moment scales from
    entering silently.
    """

    h20, h2c, h2s, h1c, h1s = map(
        sp.sympify, (p20, m2_cosine, m2_sine, m1_cosine, m1_sine)
    )
    return sp.Matrix(
        [
            [-h20 / 5 + 2 * h2c / 5, 2 * h2s / 5, 2 * h1c / 5],
            [2 * h2s / 5, -h20 / 5 - 2 * h2c / 5, 2 * h1s / 5],
            [2 * h1c / 5, 2 * h1s / 5, 2 * h20 / 5],
        ]
    )


def real_l2_tt_readout(
    tensor: Any,
    direction: Any,
    reference: Any | None = None,
) -> RealL2TTReadout:
    """Return exact arbitrary-view TT coordinates of a symmetric tensor."""

    source = sp.Matrix(tensor)
    if source.shape != (3, 3):
        raise ValueError("tensor must be 3 by 3")
    if sp.simplify(source - source.T) != sp.zeros(3):
        raise ValueError("tensor must be symmetric")
    basis = tt_polarization_basis(direction, reference)
    plus = sp.simplify(frobenius_inner_product(source, basis.plus))
    cross = sp.simplify(frobenius_inner_product(source, basis.cross))
    return RealL2TTReadout(
        basis=basis,
        source_tensor=source,
        projected_tensor=tt_basis_reconstruct(source, basis),
        normalized_plus_coordinate=plus,
        normalized_cross_coordinate=cross,
        conventional_plus_readout=sp.simplify(plus / sp.sqrt(2)),
        conventional_cross_readout=sp.simplify(cross / sp.sqrt(2)),
    )


def temporal_coefficient_rank(
    coefficients: ArrayLike,
    relative_tolerance: float = 1.0e-12,
) -> int:
    """Return the numerical rank of sampled source-coefficient time traces.

    Rows are time samples and columns are declared angular/STF components.
    Callers must remove physically irrelevant DC offsets before using this as
    an oscillatory-mode diagnostic.  Two nonzero but proportional columns
    have rank one; quadrature-phase cosine/sine columns have rank two.
    """

    values: NDArray[np.float64] = np.asarray(coefficients, dtype=float)
    if values.ndim != 2 or 0 in values.shape:
        raise ValueError("coefficients must be a nonempty two-dimensional array")
    if not np.all(np.isfinite(values)):
        raise ValueError("coefficients must be finite")
    if relative_tolerance < 0 or not np.isfinite(relative_tolerance):
        raise ValueError("relative_tolerance must be finite and nonnegative")
    singular_values = np.linalg.svd(values, compute_uv=False)
    if singular_values[0] == 0:
        return 0
    return int(np.count_nonzero(singular_values > relative_tolerance * singular_values[0]))


def linearized_l_mode_residual(
    mode_tt: Any,
    mode_rr: Any,
    mode_r: Any,
    mode: Any,
    background: Any,
    radius: Any,
    ell: Any,
) -> sp.Expr:
    """Return the exact radial-background linearized sine-Gordon residual."""

    order = _nonnegative_integer(ell, "ell")
    radial = sp.sympify(radius)
    if sp.simplify(radial) == 0:
        raise ValueError("radius must be nonzero for the reduced expression")
    value = sp.sympify(mode)
    return sp.simplify(
        sp.sympify(mode_tt)
        - sp.sympify(mode_rr)
        - 2 * sp.sympify(mode_r) / radial
        + order * (order + 1) * value / radial**2
        + sp.cos(sp.sympify(background)) * value
    )


def averaged_mode_equation_defect(
    background: Any,
    averaged_cosine: Any,
    perturbation: Any,
) -> sp.Expr:
    """Return the omitted term when ``cos(P)`` is replaced by its time average."""

    return sp.simplify(
        (sp.cos(sp.sympify(background)) - sp.sympify(averaged_cosine))
        * sp.sympify(perturbation)
    )


def regular_l_mode_origin_mismatch(
    value: Any,
    radial_derivative: Any,
    radius: Any,
    ell: Any,
) -> sp.Expr:
    """Return ``r*psi_r-ell*psi``, which vanishes at leading regular order."""

    order = _nonnegative_integer(ell, "ell")
    return sp.simplify(
        sp.sympify(radius) * sp.sympify(radial_derivative)
        - order * sp.sympify(value)
    )


def _nonnegative_integer(value: Any, name: str) -> int:
    integer = sp.sympify(value)
    if integer.is_integer is not True or integer.is_nonnegative is not True:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(integer)
