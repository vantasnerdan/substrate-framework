"""Exact divergence-free Bloch continuation of frozen attempt 0040."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import sympy as sp

from substrate_framework.verification import CheckLedger

SOURCE = Path(__file__).resolve().parents[1] / "0040" / "fourier_orbit.py"
SPEC = importlib.util.spec_from_file_location("p251_0040_fourier_orbit", SOURCE)
assert SPEC is not None and SPEC.loader is not None
fo = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fo)


def wave_set(vector):
    return set().union(*(component.keys() for component in vector))


def coulomb_potential(vector):
    rotated = fo.curl(vector)
    return tuple(
        {
            wave: coefficient / sum(entry**2 for entry in wave)
            for wave, coefficient in component.items()
        }
        for component in rotated
    )


def bloch_curl(vector, kappa, axis):
    result = ({}, {}, {})
    for wave in wave_set(vector):
        shifted = sp.Matrix([wave[j] + (kappa if j == axis else 0) for j in range(3)])
        value = sp.I * shifted.cross(
            sp.Matrix([vector[j].get(wave, 0) for j in range(3)])
        )
        for j in range(3):
            if value[j] != 0:
                result[j][wave] = sp.expand(value[j])
    return result


def bloch_divergence(vector, kappa, axis):
    return {
        wave: sp.simplify(
            sum(
                sp.I * (wave[j] + (kappa if j == axis else 0)) * vector[j].get(wave, 0)
                for j in range(3)
            )
        )
        for wave in wave_set(vector)
    }


def bloch_leray(vector, kappa, axis):
    result = ({}, {}, {})
    for wave in wave_set(vector):
        shifted = sp.Matrix([wave[j] + (kappa if j == axis else 0) for j in range(3)])
        coefficient = sp.Matrix([vector[j].get(wave, 0) for j in range(3)])
        value = coefficient - shifted * shifted.dot(coefficient) / shifted.dot(shifted)
        for j in range(3):
            if value[j] != 0:
                result[j][wave] = sp.cancel(value[j])
    return result


def conjugate_vector(vector):
    return tuple(
        {
            tuple(-entry for entry in wave): sp.conjugate(value)
            for wave, value in component.items()
        }
        for component in vector
    )


def hermitian_inner(left, right):
    return sum(
        sp.conjugate(coefficient) * right[j].get(wave, 0)
        for j in range(3)
        for wave, coefficient in left[j].items()
    )


def generators():
    angle, _ = fo.core_cage_generators(1)
    seed = fo.mul(fo.trig(0, kind="sin"), fo.trig(1, kind="sin"))
    shape = (
        fo.mul(fo.derivative(seed, 0), fo.trig(2, kind="sin")),
        fo.mul(fo.derivative(seed, 1), fo.trig(2, kind="sin")),
        fo.scale(fo.mul(seed, fo.trig(2)), -2),
    )
    return angle, shape


def matrices(axis, kappa):
    base = fo.actual_tube(sp.Integer(2), sp.Integer(1))
    omega = fo.curl(base)
    original = generators()
    slow = tuple(
        bloch_curl(coulomb_potential(field), kappa, axis) for field in original
    )
    velocities = tuple(
        bloch_leray(fo.cross(field, omega), kappa, axis) for field in slow
    )
    hessian = sp.Matrix(
        2,
        2,
        lambda i, j: sp.factor(
            hermitian_inner(velocities[i], velocities[j])
            - hermitian_inner(velocities[i], bloch_curl(velocities[j], kappa, axis))
        ),
    )
    kks = sp.Matrix(
        2,
        2,
        lambda i, j: sp.factor(
            fo.inner(omega, fo.cross(conjugate_vector(slow[i]), slow[j]))
        ),
    )
    return slow, velocities, hessian, kks


def taylor_matrix(matrix, kappa, order):
    return sp.Matrix(
        matrix.rows,
        matrix.cols,
        lambda i, j: sp.factor(
            sp.diff(matrix[i, j], kappa, order).subs(kappa, 0) / sp.factorial(order)
        ),
    )


def main() -> int:
    ledger = CheckLedger("P251-0044-Beltrami-Bloch-sector")
    kappa = sp.symbols("kappa", real=True)
    frequency = sp.symbols("nu", real=True)
    nu0 = 5 * sp.sqrt(14) / 24
    for axis in range(3):
        slow, velocities, hessian, kks = matrices(axis, kappa)
        name = "xyz"[axis]
        for index, field in enumerate(slow):
            ledger.check(
                f"{name}: generator {index} exact shifted divergence",
                all(
                    value == 0
                    for value in bloch_divergence(field, kappa, axis).values()
                ),
            )
        for index, field in enumerate(velocities):
            ledger.check(
                f"{name}: velocity {index} exact shifted divergence",
                all(
                    value == 0
                    for value in bloch_divergence(field, kappa, axis).values()
                ),
            )
        ledger.check(
            f"{name}: full energy Hessian Hermitian",
            (hessian - hessian.conjugate().T).applyfunc(sp.simplify) == sp.zeros(2),
        )
        ledger.check(
            f"{name}: full KKS anti-Hermitian",
            (kks + kks.conjugate().T).applyfunc(sp.simplify) == sp.zeros(2),
        )
        ledger.check(
            f"{name}: frozen energy recovered",
            (
                hessian.subs(kappa, 0)
                - sp.diag(sp.Rational(7, 48), sp.Rational(25, 96))
            ).applyfunc(sp.simplify)
            == sp.zeros(2),
        )
        ledger.check(
            f"{name}: frozen KKS recovered",
            kks.subs(kappa, 0)
            == sp.Matrix([[0, -sp.Rational(1, 4)], [sp.Rational(1, 4), 0]]),
        )
        h_terms = [taylor_matrix(hessian, kappa, order) for order in range(3)]
        o_terms = [taylor_matrix(kks, kappa, order) for order in range(3)]
        h_series = sum(
            (kappa**j * value for j, value in enumerate(h_terms)), sp.zeros(2)
        )
        o_series = sum(
            (kappa**j * value for j, value in enumerate(o_terms)), sp.zeros(2)
        )
        characteristic = sp.expand((h_series - sp.I * frequency * o_series).det())
        p0 = characteristic.subs({kappa: 0, frequency: nu0})
        ledger.check(f"{name}: frozen optical root", sp.simplify(p0) == 0)
        p_nu = sp.diff(characteristic, frequency).subs({kappa: 0, frequency: nu0})
        first = sp.simplify(
            -sp.diff(characteristic, kappa).subs({kappa: 0, frequency: nu0}) / p_nu
        )
        second = sp.simplify(
            -(
                sp.diff(characteristic, kappa, 2)
                + 2 * first * sp.diff(characteristic, kappa, frequency)
                + first**2 * sp.diff(characteristic, frequency, 2)
            ).subs({kappa: 0, frequency: nu0})
            / (2 * p_nu)
        )
        squared_first = sp.simplify(2 * nu0 * first)
        squared_second = sp.simplify(first**2 + 2 * nu0 * second)
        print(f"AXIS {name} H_FULL = {hessian}")
        print(f"AXIS {name} OMEGA_FULL = {kks}")
        for order in range(3):
            print(f"AXIS {name} H{order} = {h_terms[order]}")
            print(f"AXIS {name} OMEGA{order} = {o_terms[order]}")
        print(
            f"AXIS {name} NU = {nu0} + ({first}) kappa + ({second}) kappa^2 + O(kappa^3)"
        )
        print(
            f"AXIS {name} NU_SQUARED = {nu0**2} + ({squared_first}) kappa + ({squared_second}) kappa^2 + O(kappa^3)"
        )
    print(
        "Physical units: H=rho b^2 Hhat, Omega=rho ell b Omegahat, kappa=ell k, omega=(b/ell) nu."
    )
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
