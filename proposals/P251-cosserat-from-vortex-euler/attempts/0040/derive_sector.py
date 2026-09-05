"""Exact algebraic candidate comparison; no numerical spectral verdict."""

import sympy as s

from fourier_orbit import actual_tube, core_cage_generators, divergence, orbit_matrices


def main():
    a, b = s.symbols("a b", positive=True)
    for k in (0, 1, 2, 3):
        generators = core_cage_generators(k)
        tangents, hessian, kks = orbit_matrices(actual_tube(a, b), generators)
        print("k =", k)
        print("div generators =", [divergence(v) for v in generators])
        print("div tangents =", [divergence(v) for v in tangents])
        print("H =", hessian)
        print("KKS =", kks)
        print("det H =", s.factor(hessian.det()))


if __name__ == "__main__":
    main()
