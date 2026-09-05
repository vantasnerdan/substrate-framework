"""Exact periodic Euler stress transfer; no spectral truncation or solver."""

import sympy as s

from substrate_framework import euler_fourier as f
from substrate_framework.verification import CheckLedger


def vadd(*vectors):
    return tuple(f.add(*(v[i] for v in vectors)) for i in range(3))


def vscale(v, c):
    return tuple(f.scale(x, c) for x in v)


def zero(v):
    return all(not f.add(x) for x in v)


def mean(v):
    return s.Matrix([s.simplify(x.get(f.ZERO, 0)) for x in v])


def shift(v, k):
    return tuple({tuple(p[i] + k[i] for i in range(3)): c for p, c in x.items()}
                 for x in v)


def project(v, k):
    return shift(f.leray(shift(v, k)), tuple(-x for x in k))


def curl_k(v, k):
    return shift(f.curl(shift(v, k)), tuple(-x for x in k))


def derivative_k(x, axis, k):
    return f.add(f.derivative(x, axis), f.scale(x, s.I * k[axis]))


def stress(u, v):
    return s.Matrix(3, 3, lambda i, j: f.add(f.mul(u[i], v[j]),
                                          f.mul(v[i], u[j])).get(f.ZERO, 0))


def main():
    checks = CheckLedger("P251-0116-periodic-Bloch-transfer")
    sn = [f.trig(i, kind="sin") for i in range(3)]
    cs = [f.trig(i) for i in range(3)]
    u2 = (f.scale(sn[1], -1), sn[0], f.add(cs[0], cs[1]))
    u = (f.add(sn[2], f.scale(sn[1], -1)), f.add(sn[0], cs[2]),
         f.add(cs[0], cs[1]))
    checks.check("three-wave fixture is exactly solenoidal", not f.divergence(u))
    checks.check("three-wave fixture has curl eigenvalue plus one", zero(vadd(f.curl(u), vscale(u, -1))))
    checks.check("background modes all lie on the unit curl shell",
                 all(sum(x*x for x in p) == 1 for x in u for p in x))

    force2 = (u2[2], {}, u2[0])
    pf2 = f.leray(force2)
    checks.check("two-wave pressure projection is independently explicit",
                 zero(vadd(pf2, vscale((cs[1], {}, f.scale(sn[1], -1)), -1))))
    j2 = f.curl(f.cross(u2, pf2))
    expected = (f.mul(sn[0], sn[1]), f.mul(cs[0], cs[1]), {})
    checks.check("two-wave stress adjoint has nonzero off-shell component",
                 zero(vadd(j2, vscale(expected, -1))))

    force = (u[2], {}, u[0])
    jfield = f.curl(f.cross(u, f.leray(force)))
    translations = [tuple(f.derivative(x, axis) for x in u) for axis in range(3)]
    gram = s.Matrix(3, 3, lambda i, j: f.inner(translations[i], translations[j]))
    checks.check("all three exact translation-response directions are independent",
                 gram == s.eye(3))
    off = (1, 1, 0)
    checks.check("ABC extension preserves the two-wave off-shell witness",
                 any(x.get(off, 0) != 0 for x in jfield))
    h = s.Matrix([f.inner(t, jfield) for t in translations])
    coefficients = gram.inv() * h
    jperp = vadd(jfield, *(vscale(t, -coefficients[i]) for i, t in enumerate(translations)))
    xi = f.curl(jperp)
    velocity = f.leray(f.cross(xi, u))
    mu = s.factor(f.inner(jperp, jperp))
    checks.check("Gram-projected generator is divergence free", not f.divergence(xi))
    checks.check("prepared COMPLETE velocity has exactly zero initial mean", mean(velocity) == s.zeros(3, 1))
    checks.check("actual stress equals positive adjoint Gram residual",
                 s.simplify(stress(u, velocity)[2, 0] - mu) == 0 and mu > 0)
    checks.check("stress equality also follows through full Leray pairing",
                 s.simplify(f.inner(force, velocity) - f.inner(jfield, jperp)) == 0)

    potentials = [vadd(*(vscale(translations[j], gram.inv()[j, i]) for j in range(3)))
                  for i in range(3)]
    responses = [f.leray(f.cross(f.curl(a), u)) for a in potentials]
    response_matrix = s.Matrix.hstack(*(mean(v) for v in responses))
    checks.check("three independently derived mean responses form the identity", response_matrix == s.eye(3))

    # Rational Bloch fiber is an actual seven-cell supercell mode. The
    # zero-amplitude Fourier mode projects onto n-perp, not onto all R3.
    k = (s.Rational(1, 7), 0, 0)
    n = s.Matrix([1, 0, 0])
    pn = s.eye(3) - n*n.T
    vk = project(f.cross(curl_k(jperp, k), u), k)
    rks = [project(f.cross(curl_k(a, k), u), k) for a in potentials[1:]]
    rk_matrix = s.Matrix.hstack(*(mean(v) for v in rks))
    transverse = rk_matrix[1:3, :]
    correction = transverse.inv() * mean(vk)[1:3, :]
    prepared_potential = vadd(jperp, *(vscale(potentials[i+1], -correction[i]) for i in range(2)))
    prepared_generator = curl_k(prepared_potential, k)
    wk = project(f.cross(prepared_generator, u), k)
    checks.check("finite rational-fiber reaction matrix is invertible", transverse.det() != 0)
    checks.check("finite-k initial momentum is exactly zero after Kelvin-generator preparation",
                 mean(wk) == s.zeros(3, 1))
    checks.check("finite-k generator and induced velocity satisfy Bloch incompressibility",
                 not f.divergence(shift(prepared_generator, k)) and not f.divergence(shift(wk, k)))

    convective = tuple(f.add(*(f.mul(u[j], derivative_k(wk[i], j, k)) for j in range(3)))
                       for i in range(3))
    deformation = f.transport(wk, u)
    acceleration = vscale(project(vadd(convective, deformation), k), -1)
    direct_mean = mean(acceleration)
    moment_mean = -s.I * pn * stress(u, wk) * s.Matrix(k)
    checks.check("full Euler transport and pressure give exact Bloch stress transfer",
                 (direct_mean - moment_mean).applyfunc(s.simplify) == s.zeros(3, 1))
    checks.check("zero initial momentum has genuinely nonzero transverse acceleration",
                 s.simplify(direct_mean[2]) != 0)
    checks.check("using identity instead of transverse mean projector is exposed",
                 mean(project(({f.ZERO: 1}, {}, {}), k)) == s.zeros(3, 1))
    print(f"Exact prototype stress residual mu = {mu}")
    print(f"Exact finite-k transverse initial acceleration = {s.factor(direct_mean[2])}")
    print("The compact EPS-domain Gram argument is analytic, not this periodic prototype fixture.")
    raise SystemExit(checks.finish())


if __name__ == "__main__":
    main()
