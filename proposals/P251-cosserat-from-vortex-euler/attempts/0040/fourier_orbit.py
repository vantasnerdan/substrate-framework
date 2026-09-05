"""Exact finite-Fourier Euler orbit operations for this declared cell model."""

import sympy as s

ZERO = (0, 0, 0)


def add(*fields):
    result = {}
    for field in fields:
        for wave, coefficient in field.items():
            result[wave] = result.get(wave, 0)+coefficient
    return {wave: s.expand(coefficient) for wave, coefficient in result.items()
            if s.expand(coefficient) != 0}


def scale(field, coefficient):
    return add({wave: coefficient*value for wave, value in field.items()})


def mul(left, right):
    result = {}
    for k, a in left.items():
        for ell, b in right.items():
            wave = tuple(k[j]+ell[j] for j in range(3))
            result[wave] = result.get(wave, 0)+a*b
    return add(result)


def trig(axis, harmonic=1, kind="cos"):
    k = tuple(harmonic if j == axis else 0 for j in range(3))
    minus = tuple(-entry for entry in k)
    if kind == "sin":
        return {k: 1/(2*s.I), minus: -1/(2*s.I)}
    return {k: s.Rational(1, 2), minus: s.Rational(1, 2)}


def derivative(field, axis):
    return scale({wave: s.I*wave[axis]*value for wave, value in field.items()}, 1)


def divergence(vector):
    return add(*(derivative(vector[j], j) for j in range(3)))


def curl(vector):
    return tuple(add(derivative(vector[(j+2) % 3], (j+1) % 3),
                     scale(derivative(vector[(j+1) % 3], (j+2) % 3), -1))
                 for j in range(3))


def cross(left, right):
    return tuple(add(mul(left[(j+1) % 3], right[(j+2) % 3]),
                     scale(mul(left[(j+2) % 3], right[(j+1) % 3]), -1))
                 for j in range(3))


def inner(left, right):
    return s.expand(sum(mul(left[j], right[j]).get(ZERO, 0) for j in range(3)))


def leray(vector):
    waves = set().union(*(component.keys() for component in vector))
    result = ({}, {}, {})
    for wave in waves:
        norm2 = sum(value**2 for value in wave)
        dot = sum(wave[j]*vector[j].get(wave, 0) for j in range(3))
        for j in range(3):
            value = vector[j].get(wave, 0)
            if norm2:
                value -= s.Rational(wave[j], norm2)*dot
            if s.expand(value) != 0:
                result[j][wave] = s.expand(value)
    return result


def orbit_matrices(background, generators):
    """rho=1 and curl u0=u0; entries are averages per cell volume."""
    omega = curl(background)
    tangents = [leray(cross(generator, omega)) for generator in generators]
    hessian = s.Matrix(len(generators), len(generators),
                       lambda i, j: s.factor(inner(tangents[i], tangents[j])
                                              -inner(tangents[i], curl(tangents[j]))))
    kks = s.Matrix(len(generators), len(generators),
                   lambda i, j: s.factor(inner(omega, cross(generators[i], generators[j]))))
    return tangents, hessian, kks


def actual_tube(a, b):
    sx, sy = trig(0, kind="sin"), trig(1, kind="sin")
    cx, cy = trig(0), trig(1)
    return (scale(sy, -b), scale(sx, a), add(scale(cx, a), scale(cy, b)))


def core_cage_generators(k):
    sx, sy = trig(0, kind="sin"), trig(1, kind="sin")
    cx, cy = trig(0), trig(1)
    axial = trig(2, k) if k else {ZERO: s.Integer(1)}
    angle = (scale(mul(sy, axial), -1), mul(sx, axial), {})
    shape = (mul(mul(sx, cy), axial), scale(mul(mul(cx, sy), axial), -1), {})
    return angle, shape
