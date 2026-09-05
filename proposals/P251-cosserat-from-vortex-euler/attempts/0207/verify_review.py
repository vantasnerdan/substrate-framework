"""Independent canonical-input audit of the actual acoustic phase repair."""

import sympy as s

from substrate_framework import euler_fourier as ef
from substrate_framework.euler_displacement_preparation import (
    finite_displacement_cell,
    prepared_displacement,
    transverse_pair_average,
)
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0207-independent-acoustic-review")
    cell = finite_displacement_cell()
    direction = s.symbols("k0:3", real=True)
    displacement = s.symbols("D0:3", real=True)
    tau = s.symbols("tau", real=True)
    prepared = prepared_displacement(cell, direction, displacement, amplitude=tau)
    f = direction[1]*displacement[0]+direction[0]*displacement[1]

    def average(field):
        return ef.add({wave: transverse_pair_average(f*value, direction, displacement)
                       for wave, value in field.items()})

    raw = ef.add(*(ef.scale(cell.background[i], direction[0]*displacement[i])
                   for i in range(3)))
    au = ef.add(*(ef.scale(cell.background[i], direction[i]) for i in range(3)))
    internal = ef.add(prepared.material_rate[0], ef.scale(au, displacement[0]))
    checks.check("full accepted D field gives zero internal configuration phase row",
                 not average(internal))
    checks.check("the raw physical momentum gives the actual negative phase correction",
                 not ef.add(average(raw), ef.scale(ef.trig(2, kind="sin"),
                                                   -s.Rational(1, 1000))))
    checks.check("the full high-harmonic displacement return was included",
                 any(sum(q*q for q in wave) > 1 for wave in prepared.material_rate[0]))
    checks.check("configuration energy cancellation has the exact whole-law normalization",
                 transverse_pair_average(f*f, direction, displacement) == s.Rational(1, 5))
    angle, omega, time = s.symbols("theta omega t", real=True)
    g = s.Function("G")()
    displacement_return = g*(s.sin(angle-omega*time)/(2*omega)
                             +time*s.cos(angle-omega*time))
    velocity_return = g*s.cos(angle-omega*time)
    checks.check("actual configuration history solves Lin with its true velocity source",
                 s.simplify(s.diff(displacement_return, time)
                            +omega*s.diff(displacement_return, angle)-velocity_return) == 0)
    checks.check("the large initial configuration is retained rather than set to zero",
                 s.simplify(displacement_return.subs(time, 0)
                            -g*s.sin(angle)/(2*omega)) == 0)
    epsilon = s.symbols("epsilon", positive=True)
    weights = (-1/(3*epsilon), 2/(3*epsilon))
    frequencies = (epsilon, 2*epsilon)
    checks.check("positive actual frequencies admit signed control with exact zero phase moment",
                 s.simplify(sum(w/q for w, q in zip(weights, frequencies))) == 0)
    output = sum(w*s.sin(q*time) for w, q in zip(weights, frequencies))
    checks.check("that phase repair preserves the desired linear output and its true error",
                 s.simplify(s.series(output, epsilon, 0, 4).removeO()
                            -time+s.Rational(5, 6)*epsilon**2*time**3) == 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
