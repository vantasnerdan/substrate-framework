"""Conditional exact finite-Fourier Euler algebra.

The signed operator and quadratic-form identities support C-CST-008/009.
Their accepted scopes do not identify an arbitrary Fourier field with an
EPS tube or promote a numerical Fourier truncation to an exact continuum.

Scalar fields map commensurate rational wavevectors to SymPy coefficients
of exp(i k.x); vectors are triples of scalar fields. Products and ``inner``
are bilinear, without implicit conjugation: ``inner`` extracts the zero
coefficient of the real-space product, i.e. its common-cell volume AVERAGE.
Real fields satisfy f[-k] = conjugate(f[k]). No sampling, truncation, or
simulation runs here. Coefficients supplied as floats retain their floats.

For a real divergence-free Beltrami background and physical pressure
p=-rho|u|²/2, the two quadratic forms obey H-K=rho Gram(A xi), with A
given by ``material_kelvin_operator`` and mean-preserving Leray. A positive
coadjoint H therefore does not by itself prove positive material K.
"""

from collections.abc import Mapping, Sequence
from typing import Literal

import sympy as sp

Wave = tuple[int | sp.Rational, int | sp.Rational, int | sp.Rational]
ScalarField = Mapping[Wave, sp.Expr]
VectorField = Sequence[ScalarField]
ZERO = (0, 0, 0)


def add(*fields: ScalarField) -> dict:
    """Add scalar fields, expand coefficients, and discard exact zeros."""
    result = {}
    for field in fields:
        for wave, coefficient in field.items():
            result[wave] = result.get(wave, 0) + coefficient
    return {k: sp.expand(v) for k, v in result.items() if sp.expand(v) != 0}


def scale(field: ScalarField, coefficient) -> dict:
    return add({k: coefficient * value for k, value in field.items()})


def mul(left: ScalarField, right: ScalarField) -> dict:
    """Exact convolution, retaining every product mode."""
    result = {}
    for k, a in left.items():
        for ell, b in right.items():
            wave = tuple(k[j] + ell[j] for j in range(3))
            result[wave] = result.get(wave, 0) + a * b
    return add(result)


def trig(axis: int, harmonic=1, kind: Literal["sin", "cos"] = "cos") -> dict:
    """A real single-axis trigonometric field, including harmonic zero."""
    if axis not in (0, 1, 2):
        raise ValueError("axis must be 0, 1, or 2")
    harmonic = sp.sympify(harmonic)
    if harmonic.is_Rational is not True:
        raise ValueError("harmonic must be rational (commensurate modes)")
    if kind not in ("sin", "cos"):
        raise ValueError("kind must be 'sin' or 'cos'")
    wave = tuple(harmonic if j == axis else 0 for j in range(3))
    minus = tuple(-entry for entry in wave)
    positive = 1 / (2 * sp.I) if kind == "sin" else sp.Rational(1, 2)
    negative = -positive if kind == "sin" else positive
    return add({wave: positive}, {minus: negative})


def derivative(field: ScalarField, axis: int) -> dict:
    if axis not in (0, 1, 2):
        raise ValueError("axis must be 0, 1, or 2")
    return add({k: sp.I * k[axis] * value for k, value in field.items()})


def _vector(vector: VectorField) -> None:
    if len(vector) != 3:
        raise ValueError("a vector field has exactly three components")


def divergence(vector: VectorField) -> dict:
    _vector(vector)
    return add(*(derivative(vector[j], j) for j in range(3)))


def curl(vector: VectorField) -> tuple:
    _vector(vector)
    return tuple(add(derivative(vector[(j + 2) % 3], (j + 1) % 3),
                     scale(derivative(vector[(j + 1) % 3], (j + 2) % 3), -1))
                 for j in range(3))


def cross(left: VectorField, right: VectorField) -> tuple:
    _vector(left)
    _vector(right)
    return tuple(add(mul(left[(j + 1) % 3], right[(j + 2) % 3]),
                     scale(mul(left[(j + 2) % 3], right[(j + 1) % 3]), -1))
                 for j in range(3))


def inner(left: VectorField, right: VectorField) -> sp.Expr:
    """Common-cell average of left dot right; bilinear, not Hermitian."""
    _vector(left)
    _vector(right)
    return sp.expand(sum(mul(left[j], right[j]).get(ZERO, 0) for j in range(3)))


def leray(vector: VectorField, *, mean_mode: Literal["preserve", "zero"] = "preserve") -> tuple:
    """Orthogonal solenoidal projection; explicitly choose its zero mode.

    ``preserve`` is the usual periodic Leray projector. ``zero`` additionally
    projects out constant vectors, and is not interchangeable in orbit H.
    """
    _vector(vector)
    if mean_mode not in ("preserve", "zero"):
        raise ValueError("mean_mode must be 'preserve' or 'zero'")
    result = ({}, {}, {})
    for wave in set().union(*(component.keys() for component in vector)):
        norm2 = sp.sympify(sum(value**2 for value in wave))
        dot = sum(wave[j] * vector[j].get(wave, 0) for j in range(3))
        for j in range(3):
            value = vector[j].get(wave, 0)
            if norm2 != 0:
                value -= sp.sympify(wave[j]) * dot / norm2
            elif mean_mode == "zero":
                value = 0
            if sp.expand(value) != 0:
                result[j][wave] = sp.expand(value)
    return result


def transport(velocity: VectorField, field: VectorField) -> tuple:
    """The vector field (velocity dot grad) field."""
    _vector(velocity)
    _vector(field)
    return tuple(add(*(mul(velocity[j], derivative(field[i], j))
                       for j in range(3))) for i in range(3))


def material_kelvin_operator(background: VectorField, displacement: VectorField) -> tuple:
    """A xi=-P[(u.grad)xi+(D xi)^T u], using mean-preserving Leray.

    The fixed-Kelvin reconstruction equation is xi_t=A xi. Thus A is the
    NEGATIVE stationary circulation-constraint term, not the induced
    Eulerian velocity v=P(xi cross omega). For solenoidal u,xi,
    A xi=v-curl(xi cross u).
    """
    convective = transport(background, displacement)
    stationary = leray(tuple(add(convective[i], *(mul(background[j], derivative(displacement[j], i))
                                                for j in range(3))) for i in range(3)))
    return tuple(scale(component, -1) for component in stationary)


def _density(density):
    density = sp.sympify(density)
    if (density.is_positive is False or density.is_real is False
            or density.is_finite is False
            or (density.is_number and density.is_positive is not True)):
        raise ValueError("density must be positive, real, and finite")
    return density


def coadjoint_matrices(background: VectorField, generators: Sequence[VectorField], *,
                       beltrami_eigenvalue, density=1) -> tuple:
    """Return induced velocities, H, Omega for curl u=lambda u, lambda!=0.

    H_ij=rho<vi.vj-vi.curl(vj)/lambda>,
    Omega_ij=rho<omega.(xi_i cross xi_j)>; vi=P(xi_i cross omega).
    Background/generator divergence and the Beltrami equation are checked.
    Leray preserves constants. Density with symbolic sign is a positive
    caller hypothesis; an explicitly nonpositive value is rejected.
    """
    rho = _density(density)
    lam = sp.sympify(beltrami_eigenvalue)
    if lam.is_zero is not False:
        raise ValueError("Beltrami eigenvalue must be provably nonzero")
    omega = curl(background)
    if divergence(background) or any(add(omega[i], scale(background[i], -lam)) for i in range(3)):
        raise ValueError("background must satisfy div u=0 and curl u=lambda u")
    if any(divergence(generator) for generator in generators):
        raise ValueError("coadjoint generators must be divergence free")
    velocities = tuple(leray(cross(generator, omega)) for generator in generators)
    hessian = sp.Matrix(len(generators), len(generators), lambda i, j: sp.factor(
        rho * (inner(velocities[i], velocities[j]) - inner(velocities[i], curl(velocities[j])) / lam)))
    kks = sp.Matrix(len(generators), len(generators), lambda i, j: sp.factor(
        rho * inner(omega, cross(generators[i], generators[j]))))
    return velocities, hessian, kks


def material_jacobi_matrix(background: VectorField, displacements: Sequence[VectorField],
                           pressure: ScalarField, *, density=1) -> sp.Matrix:
    """K_ij=<xi_i.Hess(p)xi_j>-rho<(u.grad xi_i).(u.grad xi_j)>.

    Pressure is the PHYSICAL pressure (including density), supplied by the
    caller. This directly evaluates the material Jacobi stiffness; it does
    not substitute the coadjoint H. A stationary Euler background and real
    solenoidal displacements are the scientific application hypotheses.
    """
    rho = _density(density)
    if divergence(background) or any(divergence(xi) for xi in displacements):
        raise ValueError("background and material displacements must be divergence free")
    convection = tuple(transport(background, xi) for xi in displacements)
    hess_p = [[derivative(derivative(pressure, i), j) for j in range(3)] for i in range(3)]

    def entry(a, b):
        pressure_part = sum(mul(mul(displacements[a][i], displacements[b][j]), hess_p[i][j])
                            .get(ZERO, 0) for i in range(3) for j in range(3))
        return sp.factor(pressure_part - rho * inner(convection[a], convection[b]))

    return sp.Matrix(len(displacements), len(displacements), entry)
