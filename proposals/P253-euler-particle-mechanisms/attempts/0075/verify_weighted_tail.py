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

    # Derive the equatorial stationary residual of the complete fixed-frame
    # inverse.  Only the even Taylor jets that can survive at z=0 are needed;
    # the result is therefore an exact projection, not a harmonic truncation.
    x0, y0, z0, radius0 = sp.symbols("x0 y0 z0 radius0", real=True)
    radius0 = sp.sqrt(x0**2 + y0**2)
    tau0 = z0 / radius0
    c0, f0, f2, h2, g0, g2 = sp.symbols(
        "c0 f0 f2 h2 g0 g2", real=True
    )
    fjet = f0 + f2 * tau0**2 / 2
    gjet = g0 + g2 * tau0**2 / 2
    hprime = h2 * tau0
    radial_profile = (
        -c0 + tau0**2 * h2 + 2 * tau0 * hprime + fjet - g2 * tau0
    )
    angular_profile = -c0 + tau0 * hprime - fjet + g2 * tau0
    axial_profile = (
        -tau0 * g2 * tau0 - tau0 * h2 - gjet - hprime
    )
    radial_velocity = x0 * radial_profile / radius0**3
    angular_velocity = y0 * angular_profile / radius0**3
    axial_velocity = x0 * axial_profile / radius0**3
    velocity = sp.ImmutableMatrix(
        [
            radial_velocity * x0 / radius0
            - angular_velocity * y0 / radius0,
            radial_velocity * y0 / radius0
            + angular_velocity * x0 / radius0,
            axial_velocity,
        ]
    )
    coordinates = (x0, y0, z0)
    acceleration = sp.ImmutableMatrix(
        [
            sum(
                velocity[j] * sp.diff(velocity[i], coordinates[j])
                for j in range(3)
            )
            for i in range(3)
        ]
    )
    curl_acceleration = sp.ImmutableMatrix(
        [
            sp.diff(acceleration[2], y0) - sp.diff(acceleration[1], z0),
            sp.diff(acceleration[0], z0) - sp.diff(acceleration[2], x0),
            sp.diff(acceleration[1], x0) - sp.diff(acceleration[0], y0),
        ]
    )
    positive_x = sp.symbols("positive_x", positive=True)
    equatorial = sp.factor(
        curl_acceleration[1].subs({x0: positive_x, y0: 0, z0: 0})
        * positive_x**6
    )
    reduced = sp.factor(equatorial.subs({h2: -f0, f2: -f0}))
    target = g0 * (-10 * c0 + 7 * f0) + g2 * (-5 * c0 + 3 * f0)
    assert sp.simplify(reduced - target) == 0

    integral0, integral1 = sp.symbols("integral0 integral1", positive=True)
    aa = -10 * c0 + 7 * f0
    bb = -5 * c0 + 3 * f0
    residual_brace = (integral0 * aa - integral1 * bb) / (4 * sp.pi**2)
    positive_split = (
        (integral0 - integral1) * bb + integral0 * (aa - bb)
    ) / (4 * sp.pi**2)
    assert sp.simplify(residual_brace - positive_split) == 0
    assert sp.simplify(aa - bb - (-5 * c0 + 4 * f0)) == 0

    print("P253/0075 exact pressure/order checks: all passed")


if __name__ == "__main__":
    main()
