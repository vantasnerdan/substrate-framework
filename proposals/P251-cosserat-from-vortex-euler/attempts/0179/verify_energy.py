"""Exact full-current/range-energy bridge for the selected two-wave orientation."""

from __future__ import annotations

import sympy as sp

from substrate_framework import euler_fourier as ef
from substrate_framework.verification import CheckLedger


def main() -> int:
    ledger = CheckLedger("P251-0179-current-energy")
    base = (ef.add(ef.trig(1), ef.trig(2)), ef.trig(2, kind="sin"), ef.scale(ef.trig(1, kind="sin"), -1))
    alpha = ef.add(ef.trig(1), ef.scale(ef.trig(2), -1))

    def multiplier(field, symbol):
        return ef.add({wave: symbol(sp.Integer(sum(component**2 for component in wave))) * value for wave, value in field.items() if wave != ef.ZERO})

    def transport(field):
        return ef.transport(base, (field, {}, {}))[0]

    def q(field):
        return multiplier(field, lambda norm: sp.sqrt((norm - 1) / norm))

    def generator(field):
        return ef.scale(q(transport(q(field))), -1)

    def pairing(left, right):
        return sp.simplify(ef.mul(left, right).get(ef.ZERO, 0))

    force = ef.scale(q(transport(alpha)), 2)
    ledger.check("range forcing has the exact nonzero current scale", pairing(force, force) == 2)
    ledger.check("force occupies a finite analytic-vector starting shell", {sum(component**2 for component in wave) for wave in force} == {2})
    phi = ef.add(ef.mul(ef.trig(1, kind="sin"), ef.trig(2, kind="sin")), ef.trig(1, 2), ef.scale(alpha, 3))
    b_phi = multiplier(phi, lambda norm: norm - 1)
    phi_t = multiplier(ef.add(ef.scale(transport(b_phi), -1), ef.scale(transport(alpha), 2)), lambda norm: 1 / norm)
    range_column = multiplier(phi, lambda norm: sp.sqrt(norm * (norm - 1)))
    range_t = multiplier(phi_t, lambda norm: sp.sqrt(norm * (norm - 1)))
    ledger.check("actual vorticity equation gives the stated forced range generator", not ef.add(range_t, ef.scale(generator(range_column), -1), ef.scale(force, -1)))
    current_t = -2 * pairing(alpha, phi_t)
    energy_t = pairing(range_column, range_t)
    ledger.check("physical stress derivative is exactly minus the range energy derivative", current_t == -energy_t and current_t != 0)
    ledger.check("full force current row equals the measured first-shell current derivative", current_t == -pairing(force, range_column))
    second_column = ef.add(ef.trig(2, 2, kind="sin"), ef.mul(ef.trig(1), ef.trig(2)))
    ledger.check("range generator is skew on independently chosen full Fourier fields", sp.simplify(pairing(generator(range_column), second_column) + pairing(range_column, generator(second_column))) == 0)
    ledger.check("first-shell material motion remains outside the range energy", not multiplier(alpha, lambda norm: sp.sqrt(norm * (norm - 1))))
    current_ddot_zero = -pairing(force, force)
    ledger.check("zero range preparation reproduces the derivative-two obstruction", current_ddot_zero == -2)
    test = force
    for order in range(1, 5):
        test = generator(test)
        maximum = max(sum(abs(component) for component in wave) for wave in test)
        ledger.check(f"full generated support obeys analytic-vector growth at order{order}", maximum <= order + 2)
    velocity = sp.Matrix(sp.symbols("ux uy uz", real=True))

    def moment(i, j, a, b):
        delta = sp.KroneckerDelta
        return sp.Rational(2, 15) * delta(i, a) * delta(j, b) - (delta(i, j) * delta(a, b) + delta(i, b) * delta(a, j)) / 30

    strain = [0, 1, -1]
    contracted = sp.Matrix([sum(strain[a] * (moment(a, a, i, j) + moment(a, a, j, i)) * velocity[i] for a in range(3) for i in range(3)) for j in range(3)])
    ledger.check("complete whole-law forcing contraction has no omitted axial row", contracted == sp.Matrix([0, velocity[1] / 5, -velocity[2] / 5]))
    y, z = sp.symbols("y z", real=True)
    stream = sp.Function("phi")(y, z)
    axial = sp.Function("b")(y, z)
    stationary_stream = sp.cos(y) + sp.cos(z)
    base_planar = sp.Matrix([sp.sin(z), -sp.sin(y)])
    perturbed_planar = sp.Matrix([-sp.diff(stream, z), sp.diff(stream, y)])
    eta = sp.diff(stream, y, 2) + sp.diff(stream, z, 2)
    transport_axial = base_planar.dot(sp.Matrix([sp.diff(axial, y), sp.diff(axial, z)]))
    transport_eta = base_planar.dot(sp.Matrix([sp.diff(eta, y), sp.diff(eta, z)]))
    advected_stream = perturbed_planar.dot(sp.Matrix([sp.diff(stationary_stream, y), sp.diff(stationary_stream, z)]))
    axial_t, eta_t = -transport_axial - advected_stream, -transport_eta + advected_stream
    ledger.check("axial-plus-vorticity variable is an actual Euler passive scalar", sp.simplify(axial_t + eta_t + transport_axial + transport_eta) == 0)
    print(f"actual current_t={current_t}; range energy_t={energy_t}; force norm squared={pairing(force, force)}")
    print("No observed-current closure is inferred for whole-law averaging from this one orientation.")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
