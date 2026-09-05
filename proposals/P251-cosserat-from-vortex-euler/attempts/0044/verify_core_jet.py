"""Independent real-space Bloch-return and physical-angle coordinate check."""

from __future__ import annotations

import sympy as sp

from substrate_framework.verification import CheckLedger


def main() -> int:
    ledger = CheckLedger("P251-0044-core-jet-observable")
    x, y, z = sp.symbols("x y z", real=True)
    kx, ky, kz = sp.symbols("kx ky kz", real=True)
    coordinates = (x, y, z)
    wave = sp.Matrix([kx, ky, kz])
    potentials = (
        sp.Matrix(
            [
                sp.sin(x) * sp.sin(z),
                sp.sin(y) * sp.sin(z),
                (sp.cos(x) + sp.cos(y)) * sp.cos(z),
            ]
        )
        / 2,
        sp.Matrix(
            [-sp.sin(x) * sp.cos(y) * sp.cos(z), sp.cos(x) * sp.sin(y) * sp.cos(z), 0]
        ),
    )

    def curl(vector):
        return sp.Matrix(
            [
                sp.diff(vector[(j + 2) % 3], coordinates[(j + 1) % 3])
                - sp.diff(vector[(j + 1) % 3], coordinates[(j + 2) % 3])
                for j in range(3)
            ]
        )

    def curl_bloch(vector):
        return curl(vector) + sp.I * wave.cross(vector)

    for index, potential in enumerate(potentials):
        generator = curl_bloch(potential)
        divergence = sum(
            sp.diff(generator[j], coordinates[j]) for j in range(3)
        ) + sp.I * wave.dot(generator)
        ledger.check(
            f"generator {index}: real-space Bloch divergence vanishes",
            sp.simplify(divergence) == 0,
        )
        local_angle = sp.simplify(curl_bloch(generator)[2] / 2)
        for section, sign in ((0, 1), (sp.pi, -1)):
            expected = sign * (1 + (kx**2 + ky**2) / 2) if index == 0 else 0
            ledger.check(
                f"generator {index}: exact core rotation at section {section}",
                sp.simplify(local_angle.subs({x: 0, y: 0, z: section}) - expected) == 0,
            )
    inertia, stiffness, i2, k2, coefficient = sp.symbols("I0 K0 I2 K2 r", real=True)
    transformed_i2 = i2 - 2 * coefficient * inertia
    transformed_k2 = k2 - 2 * coefficient * stiffness
    ledger.check(
        "physical-angle coordinate map preserves normalized gradient",
        sp.simplify(
            transformed_k2
            - stiffness * transformed_i2 / inertia
            - k2
            + stiffness * i2 / inertia
        )
        == 0,
    )
    relative_map = sp.cos(sp.pi * kz / 2) * (1 + (kx**2 + ky**2) / 2)
    ledger.check(
        "finite-section axial relative-angle correction",
        sp.diff(relative_map, kz, 2).subs({kx: 0, ky: 0, kz: 0}) / 2 == -(sp.pi**2) / 8,
    )
    print(
        "Local section rotation: plus/minus (1+k_perp^2/2) times its local Bloch amplitude."
    )
    print(
        "Relative angle at the pair midpoint: chi=2 cos(pi k_z/2)(1+k_perp^2/2) q_mid."
    )
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
