"""Actual solenoidal initial-data repair for the nonplanar SH response."""

import sympy as s

from substrate_framework import euler_fourier as f
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0144-matched-SH-data")
    gx, gy, kx, ky, kz, amplitude = s.symbols("gx gy kx ky kz amplitude", real=True)
    horizontal_wave = s.Matrix([gx+kx, gy+ky])
    norm2 = horizontal_wave.dot(horizontal_wave)
    vertical_velocity = -s.I*kz*amplitude
    horizontal_return = s.I*kz**2*horizontal_wave*amplitude/norm2
    ledger.check("physical axial preparation has an exact bounded microscopic divergence return",
                 s.simplify(s.I*horizontal_wave.dot(horizontal_return)
                            +s.I*kz*vertical_velocity) == 0)
    ledger.check("the return starts at second order and has no hidden first axial derivative",
                 horizontal_return.subs(kz, 0) == s.zeros(2, 1)
                 and horizontal_return.diff(kz).subs(kz, 0) == s.zeros(2, 1))
    pressure_multiplier = 1/(norm2+kz*kz)
    ledger.check("the actual full-pressure cell denominator has no first axial derivative",
                 s.diff(pressure_multiplier, kz).subs(kz, 0) == 0
                 and s.diff(pressure_multiplier, kz, 2).subs(kz, 0) != 0)
    wave = s.Matrix([kx, ky, kz])
    planar_wave = s.Matrix([kx, ky, 0])
    sh = s.Matrix([-ky, kx, 0])
    full_pressure = s.eye(3)-wave*wave.T/wave.dot(wave)
    planar_pressure = s.eye(3)-planar_wave*planar_wave.T/planar_wave.dot(planar_wave)
    ledger.check("the actual harmonic pressure projections agree on the SH mean observation",
                 all(s.simplify(x) == 0 for x in full_pressure*sh-sh)
                 and all(s.simplify(x) == 0 for x in planar_pressure*sh-sh))

    psi = f.mul(f.trig(0), f.trig(1))
    base = (f.scale(f.derivative(psi, 1), -1), f.derivative(psi, 0), {})
    # This pure-planar base is not Beltrami. Derive its actual pressure
    # from Delta p=-div[(v.grad)v], rather than copying a sign convention.
    acceleration = f.transport(base, base)
    pressure = {g: value/s.sympify(sum(x*x for x in g))
                for g, value in f.divergence(acceleration).items() if g != f.ZERO}
    ledger.check("the comparison pressure is derived from its actual Euler acceleration",
                 all(not f.add(acceleration[j], f.derivative(pressure, j)) for j in range(3)))
    source_row = f.derivative(pressure, 0)
    axial_row = f.scale(base[0], -1)
    transport_row = f.transport(base, ({}, {}, axial_row))[2]
    ledger.check("the matched axial row solves the actual stationary Euler transport equation",
                 not f.add(transport_row, f.scale(source_row, -1)))
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
