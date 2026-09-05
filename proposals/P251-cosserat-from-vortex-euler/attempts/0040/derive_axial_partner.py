"""Exact continuation to an axial displacement with incompressible return."""

import sympy as s

from fourier_orbit import (
    actual_tube, core_cage_generators, derivative, divergence, mul,
    orbit_matrices, scale, trig,
)


def main():
    a, b = s.symbols("a b", positive=True)
    f = mul(trig(0, kind="sin"), trig(1, kind="sin"))
    for k in (1, 2, 3):
        angle, _ = core_cage_generators(k)
        axial_sin, axial_cos = trig(2, k, "sin"), trig(2, k)
        partner = (scale(mul(derivative(f, 0), axial_sin), s.Rational(1, k)),
                   scale(mul(derivative(f, 1), axial_sin), s.Rational(1, k)),
                   scale(mul(f, axial_cos), -s.Rational(2, k**2)))
        tangents, hessian, kks = orbit_matrices(actual_tube(a, b), (angle, partner))
        print("k =", k)
        print("div generators =", [divergence(v) for v in (angle, partner)])
        print("div tangents =", [divergence(v) for v in tangents])
        print("H =", hessian)
        print("KKS =", kks)
        print("det H =", s.factor(hessian.det()))


if __name__ == "__main__":
    main()
