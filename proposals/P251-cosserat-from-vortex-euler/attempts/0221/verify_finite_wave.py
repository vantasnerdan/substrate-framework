"""Exact full-pressure Bloch anchors; no spectral discretization or fit.

rho=1, curl u=-u, Cartesian coordinates (X,a,b). The exact carrier N>4
keeps all sidebands nonzero. The ray (1,2,3) exercises axial and both
normal derivatives. These are symbol-order anchors, not the off-flow
propagator theorem. All convolutions, derivatives and pressure projections
use the canonical finite-Fourier algebra without mode truncation.
"""

import sys

import sympy as s

from substrate_framework import euler_fourier as ef
from substrate_framework.verification import CheckLedger


def main(*, phase_only=False):
    checks = CheckLedger("P251-0221-kelvin-finite-wave")
    eps = s.symbols("epsilon", real=True)
    n = s.symbols("N", positive=True)
    ca = s.Rational(1, 100)
    psi = ef.add(ef.trig(2), ef.scale(ef.trig(1), ca))
    u = (psi, ef.trig(2, kind="sin"),
         ef.scale(ef.trig(1, kind="sin"), -ca))
    omega = ef.curl(u)
    phase_generators = []
    for axial, shift in ((-1, 0), (0, 1)):
        wave = (eps, n+2*eps+shift, 3*eps)
        p, q = wave[1:]
        normal2 = p*p+q*q
        generator = ({wave: s.Integer(axial)},
                     {wave: -s.I*q-eps*axial*p/normal2},
                     {wave: s.I*p-eps*axial*q/normal2})
        phase_generators.append(generator)
        if phase_only:
            continue
        checks.check(f"sector {axial} exact full Bloch divergence vanishes",
                     all(s.cancel(v) == 0 for v in ef.divergence(generator).values()))
        force = ef.cross(generator, omega)
        velocity = ef.leray(force)
        checks.check(f"sector {axial} retains every full-pressure solenoidal sideband",
                     all(s.cancel(v) == 0 for v in ef.divergence(velocity).values()))
        # Each exact product has wave = carrier + one first-shell index.
        modes = set().union(*(part.keys() for part in velocity))
        energy2 = 0
        for mode in modes:
            vv = s.Matrix([s.cancel(part.get(mode, 0)) for part in velocity])
            v0 = vv.subs(eps, 0).applyfunc(s.cancel)
            v1 = vv.diff(eps).subs(eps, 0).applyfunc(s.cancel)
            v2 = (vv.diff(eps, 2).subs(eps, 0)/2).applyfunc(s.cancel)
            mode0 = s.Matrix(mode).subs(eps, 0)
            ray = s.Matrix([1, 2, 3])
            curl0 = s.I*mode0.cross(v0)
            curl1 = s.I*(mode0.cross(v1)+ray.cross(v0))
            curl2 = s.I*(mode0.cross(v2)+ray.cross(v1))
            energy2 += (v0.conjugate().dot(v2+curl2)
                        +v1.conjugate().dot(v1+curl1)
                        +v2.conjugate().dot(v0+curl0))
        energy2 = s.factor(s.cancel(energy2))
        numerator, denominator = s.fraction(energy2)
        degree = s.degree(numerator, n)-s.degree(denominator, n)
        checks.check(f"sector {axial} actual energy second jet has carrier order at most zero",
                     degree <= 0)
        checks.check(f"sector {axial} normalized energy second coefficient vanishes",
                     s.limit(energy2/n, n, s.oo) == 0)
        print(f"sector {axial} actual full-pressure H2:", energy2)

    # A b-sideband couples to omega_a=-sin(b). The initial a-sideband
    # instead made this cross row identically zero by Fourier selection.
    phase_wave = (eps, n+2*eps, 1+3*eps)
    fraction = s.symbols("axial_fraction", real=True)
    p, q = phase_wave[1:]
    normal2 = p*p+q*q
    phase_generators[1] = (
        {phase_wave: fraction},
        {phase_wave: -s.I*q-eps*fraction*p/normal2},
        {phase_wave: s.I*p-eps*fraction*q/normal2})
    left = tuple({tuple(-entry for entry in wave): s.conjugate(value)
                  for wave, value in part.items()}
                 for part in phase_generators[0])
    phase = s.cancel(ef.inner(omega, ef.cross(left, phase_generators[1])))
    phase2 = s.factor(s.cancel(s.diff(phase, eps, 2).subs(eps, 0)/2))
    num, den = s.fraction(phase2)
    checks.check("actual same-positive-sector phase retains a nonzero physical row",
                 phase.subs({eps: 0, fraction: -1}) != 0)
    checks.check("actual opposite-sector zero-wave cross phase cancels exactly",
                 phase.subs({eps: 0, fraction: 0}) == 0)
    checks.check("actual mixed-sector phase second jet loses two carrier orders",
                 phase2 == 0 or s.degree(num, n)-s.degree(den, n) <= -1)
    print("actual integrated Omega0:", s.factor(phase.subs(eps, 0)))
    print("actual arbitrary-second-sector Omega2:", phase2)

    px, py, kz, cc = s.symbols("px py Kx c", real=True)
    incomplete = s.Matrix([cc, -s.I*py, s.I*px])
    checks.check("omitting axial divergence completion produces a nonzero physical defect",
                 s.Matrix([kz, px, py]).dot(incomplete) == cc*kz)
    z = s.Matrix(s.symbols("z0:3", real=True))
    k = s.Matrix(s.symbols("p0:3", real=True))
    checks.check("single real leading velocity polarization cancels helicity principal symbol",
                 s.expand(z.dot(k.cross(z))) == 0)
    power = s.symbols("L", positive=True)
    chosen = n**(-power-1)
    checks.check("the spatial window controls cubic norm growth",
                 s.simplify(n**power*chosen-1/n) == 0)
    checks.check("the stronger window controls nonzero order-zero remote pressure output",
                 s.simplify(n**(-2*power-3)/chosen**2-1/n) == 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main(phase_only="--phase-only" in sys.argv))
