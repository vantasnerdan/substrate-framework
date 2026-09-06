"""Exact algebraic checks for the P253/0075 pressure and decay ledger."""

from __future__ import annotations

import sympy as sp

from substrate_framework.euler_asymptotic_tails import (
    newton_pressure_kernel,
    pressure_stress_gradient,
    tail_evolution_order_ledger,
)


def main() -> None:
    x, y, z, scale = sp.symbols("x y z scale", positive=True)
    point = sp.ImmutableMatrix([x, y, z])
    r = sp.sqrt(point.dot(point))
    green = 1 / (4 * sp.pi * r)
    derived = sp.ImmutableMatrix(
        [[sp.diff(green, point[i], point[j]) for j in range(3)] for i in range(3)]
    )
    kernel = newton_pressure_kernel(point)
    assert all(
        sp.simplify(derived[i, j] - kernel[i, j]) == 0
        for i in range(3)
        for j in range(3)
    )
    assert sp.simplify(sp.trace(kernel)) == 0

    scaled_kernel = kernel.subs(
        {x: scale * x, y: scale * y, z: scale * z}, simultaneous=True
    )
    assert all(
        sp.simplify(scaled_kernel[i, j] - kernel[i, j] / scale**3) == 0
        for i in range(3)
        for j in range(3)
    )

    stress = sp.diag(1, 2, 4)
    gradient = pressure_stress_gradient(point, stress)
    scaled_gradient = gradient.subs(
        {x: scale * x, y: scale * y, z: scale * z}, simultaneous=True
    )
    assert all(
        sp.simplify(scaled_gradient[i] - gradient[i] / scale**4) == 0
        for i in range(3)
    )

    radius, cutoff = sp.symbols("radius cutoff", positive=True)
    assert sp.integrate(1 / radius, (radius, 1, cutoff)) == sp.log(cutoff)

    a = sp.ImmutableMatrix([0, 0, 1])
    homogeneous_velocity = a.cross(point) / r**3
    acceleration = sp.ImmutableMatrix(
        [
            sum(
                homogeneous_velocity[j]
                * sp.diff(homogeneous_velocity[i], point[j])
                for j in range(3)
            )
            for i in range(3)
        ]
    )
    scaled_acceleration = acceleration.subs(
        {x: scale * x, y: scale * y, z: scale * z}, simultaneous=True
    )
    assert all(
        sp.simplify(scaled_acceleration[i] - acceleration[i] / scale**5) == 0
        for i in range(3)
    )

    ledger = tail_evolution_order_ledger()
    assert (
        ledger.stress,
        ledger.pressure,
        ledger.pressure_gradient,
        ledger.transport,
        ledger.velocity_time_derivative,
    ) == (4, 3, 4, 5, 4)
    assert ledger.stress_first_moment_is_logarithmic

    print("P253/0075 exact pressure/order checks: all passed")


if __name__ == "__main__":
    main()
