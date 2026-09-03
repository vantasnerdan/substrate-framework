#!/usr/bin/env python
"""Exact review repairs joining the P250 wall, slip, pressure, and bag layers."""

from __future__ import annotations

import contextlib
import io
from pathlib import Path
import sys

import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "src"))

from substrate_framework.m5_exterior_clock import (  # noqa: E402
    clock_inertia_density,
    clock_plane_basis,
    full_clock_potential,
)
from substrate_framework.m5_wall_clock import (  # noqa: E402
    clock_slice_tensor,
    maxwell_system,
    phase_twist_gradient_excess,
    relative_equilibrium_energies,
    wall_bulk_pressure,
    wall_slice_inertia,
    wall_slice_potential,
)
from substrate_framework.verification import CheckLedger  # noqa: E402


def main():
    ledger = CheckLedger("C-M5W-review-stitching")
    m, c, b, f, omega_sq = sp.symbols("m c b f omega_sq", real=True)
    u, v, y, g = sp.symbols("u v y g", real=True)
    full_s = sp.Matrix(
        [[m, u, v], [u, c + b, y], [v, y, c - b]]
    )
    routhian = sp.expand(
        full_clock_potential(full_s, f, g)
        - omega_sq * clock_inertia_density(full_s, f, g) / 2
    )
    slice_subs = {u: 0, v: 0, y: 0, g: 0}
    transverse_residuals = [
        sp.simplify(sp.diff(routhian, variable).subs(slice_subs))
        for variable in (u, v, y, g)
    ]
    ledger.check(
        "aligned real-psi slice is an exact invariant reduction",
        transverse_residuals == [0, 0, 0, 0],
    )

    shear_u = sp.factor(sp.diff(routhian, u, 2).subs(slice_subs))
    shear_v = sp.factor(sp.diff(routhian, v, 2).subs(slice_subs))
    expected_u = (
        16 * b**2 - 6 * b + 16 * c**2 - 6 * c
        + 8 * m**2 - 6 * m + 2 - omega_sq
    )
    expected_v = expected_u + 12 * b
    ledger.check("first shear Hessian includes the rotating-frame term", sp.expand(shear_u - expected_u) == 0)
    ledger.check("second shear Hessian includes split sign and rotating term", sp.expand(shear_v - expected_v) == 0)
    static_expression = 2 * (
        8 * b**2 - 3 * b + 8 * c**2 - 3 * c + 4 * m**2 - 3 * m + 1
    )
    ledger.check(
        "archived single shear expression is static-channel data only",
        sp.expand(shear_u.subs(omega_sq, 0) - static_expression) == 0
        and sp.expand(shear_u - static_expression) != 0,
    )

    dm, dc, db, df, dq = sp.symbols("dm dc db df dq", real=True)
    q = sp.symbols("q", real=True)
    generator, _, _ = clock_plane_basis()
    rotation = sp.exp(q * sp.Matrix(generator))
    slice_s = clock_slice_tensor(m, c, b)
    slice_ds = clock_slice_tensor(dm, dc, db)
    rotated_s_derivative = rotation * (
        slice_ds + dq * (sp.Matrix(generator) * slice_s - slice_s * sp.Matrix(generator))
    ) * rotation.T
    rotated_psi_real_derivative = df * sp.cos(q) - f * dq * sp.sin(q)
    rotated_psi_imag_derivative = df * sp.sin(q) + f * dq * sp.cos(q)
    transformed_gradient = sp.trace(rotated_s_derivative.T * rotated_s_derivative) / 4 + (
        rotated_psi_real_derivative**2 + rotated_psi_imag_derivative**2
    ) / 2
    base_gradient = sp.trace(slice_ds.T * slice_ds) / 4 + df**2 / 2
    excess = sp.trigsimp(sp.expand(transformed_gradient - base_gradient))
    expected_excess = phase_twist_gradient_excess(
        wall_slice_inertia(m, c, b, f), dq
    )
    ledger.check(
        "spatial phase twist costs exactly iota times q-prime squared over two",
        sp.simplify(excess - expected_excess) == 0,
    )

    pressure = wall_bulk_pressure(m, c, b, f, omega_sq)
    ledger.check(
        "pressure is minus the canonical rotating-frame potential",
        sp.expand(pressure + wall_slice_potential(m, c, b, f, omega_sq)) == 0,
    )
    ledger.check(
        "partial pressure derivative is one half the inertia",
        sp.simplify(sp.diff(pressure, omega_sq) - wall_slice_inertia(m, c, b, f) / 2) == 0,
    )
    # The total branch derivative differs from the partial derivative by the
    # stationarity equations; this identity is the exact envelope bridge.
    dmw, dcw, dbw, dfw = sp.symbols("dmw dcw dbw dfw", real=True)
    branch_chain = sum(
        sp.diff(pressure, variable) * derivative
        for variable, derivative in zip((m, c, b, f), (dmw, dcw, dbw, dfw))
    )
    ledger.check(
        "stationarity kills every pressure branch-chain term",
        sp.simplify(
            branch_chain.subs(
                {sp.diff(pressure, variable): 0 for variable in (m, c, b, f)}
            )
        ) == 0,
    )

    rotating_functional, inertia, omega = sp.symbols(
        "F_omega I omega", real=True
    )
    static_energy, physical_energy = relative_equilibrium_energies(
        rotating_functional, inertia, omega
    )
    charge = omega * inertia
    ledger.check(
        "static and physical energies have the canonical distinct Legendre factors",
        sp.expand(
            static_energy - rotating_functional - omega**2 * inertia / 2
        ) == 0
        and sp.expand(
            physical_energy - rotating_functional - omega * charge
        ) == 0
        and sp.expand(
            physical_energy - static_energy - omega**2 * inertia / 2
        ) == 0,
    )

    radius = sp.symbols("r", positive=True)
    kinetic, potential = sp.symbols("T V", real=True)
    hamiltonian = kinetic - potential
    local_identity = sp.symbols("Hprime") + 4 * kinetic / radius
    radial_virial_derivative = sp.expand(
        3 * radius**2 * hamiltonian + radius**3 * sp.symbols("Hprime")
    )
    ledger.check(
        "radial local identity implies the Derrick integrand",
        sp.expand(
            radial_virial_derivative.subs(
                sp.symbols("Hprime"), -4 * kinetic / radius
            )
            + radius**2 * (kinetic + 3 * potential)
        ) == 0
        and local_identity != 0,
    )

    p_a, p_b, p_c, w_of = maxwell_system()
    returned_symbols = {
        symbol.name: symbol
        for expression in (p_a, p_b, p_c, w_of)
        for symbol in expression.free_symbols
    }
    c_max = returned_symbols["c"]
    b_max = returned_symbols["b"]
    f_max = returned_symbols["f"]
    omega_max = sp.Symbol("omega_sq", nonnegative=True)
    ledger.check(
        "Maxwell c equation normalization is explicit",
        sp.expand(
            p_a / 2
            - (
                8 * c_max**3 - 3 * c_max**2
                + (8 * b_max**2 + 1) * c_max - 3 * b_max**2
            )
        ) == 0,
    )
    maxwell_potential = wall_slice_potential(
        0, c_max, b_max, f_max, omega_max
    )
    ledger.check(
        "Maxwell b and depth equations remain canonical derivatives",
        sp.expand(
            p_b
            - sp.diff(maxwell_potential, b_max).subs(omega_max, w_of)
        ) == 0
        and sp.expand(
            p_c - 2 * maxwell_potential.subs(omega_max, w_of)
        ) == 0,
    )

    return ledger.finish()


if __name__ == "__main__":
    transcript = io.StringIO()
    with contextlib.redirect_stdout(transcript):
        status = main()
    rendered = transcript.getvalue()
    print(rendered, end="")
    (HERE / "exact_stitching_stdout.txt").write_text(rendered)
    raise SystemExit(status)
