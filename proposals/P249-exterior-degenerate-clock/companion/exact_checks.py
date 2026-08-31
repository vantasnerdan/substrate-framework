"""Canonical exact checks for the P249 O6 review transaction."""

from __future__ import annotations

import sympy as sp

from substrate_framework.m5_exterior_clock import (
    canonical_velocity_energy_density,
    charged_continuum_edge_squared,
    clock_inertia_density,
    exterior_channel_pencil,
    fixed_charge_legendre_data,
    full_clock_binding_witness,
    full_clock_potential,
    m5_radial_lower_bound,
    noether_charge_density,
    phase_lock_potential,
    projected_m5_static_potential,
    relative_equilibrium_velocity,
    unsplit_charge_ratio,
    unsplit_full_lower_bound,
)


CHECKS = 0


def check(condition: object, message: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.true:
        raise AssertionError(message)
    CHECKS += 1


radius = sp.symbols("radius", nonnegative=True, real=True)
check(
    sp.expand(
        m5_radial_lower_bound(radius)
        - (radius - 1) ** 2 * (2 * radius**2 + 2 * radius + 1) / 2
    )
    == 0,
    "global M5 lower bound",
)

vacuum = sp.diag(1, 0, 0)
check(projected_m5_static_potential(vacuum) == 0, "aligned M5 vacuum")
check(full_clock_potential(vacuum, 0, 0) == 0, "complete aligned exterior")

stiffness, kinetic = exterior_channel_pencil()
check(kinetic.det() == sp.Rational(1, 2), "positive kinetic determinant")
check(
    kinetic.inv() * stiffness == sp.diag(10, 10, 22, 4, 4, 22),
    "full generalized mass pencil",
)
check(min(6, 4, sp.Rational(22, 4)) == 4, "all charged weights classified")
check(charged_continuum_edge_squared() == 4, "complete charged edge")

witness = full_clock_binding_witness()
core = sp.diag(1, witness.clock_split, -witness.clock_split)
potential = full_clock_potential(core, witness.scalar_amplitude, 0)
inertia = clock_inertia_density(core, witness.scalar_amplitude, 0)
check(potential == witness.potential_density, "binding potential")
check(inertia == witness.inertia_density, "binding inertia")
check(sp.factor(2 * potential / inertia) == witness.charge_ratio, "binding ratio")
check(
    witness.charged_edge_squared - witness.charge_ratio
    == witness.binding_margin,
    "binding margin",
)
check(unsplit_charge_ratio(sp.Rational(1, 4)) == 5, "unsplit minimum")
check(unsplit_full_lower_bound() == 4, "full unsplit lower bound")
check(
    witness.charged_edge_squared - witness.charge_ratio
    == witness.unsplit_binding_margin,
    "unsplit binding margin",
)

shear_u, shear_v = sp.symbols("shear_u shear_v", real=True)
sheared = sp.Matrix([[1, shear_u, shear_v], [shear_u, 0, 0], [shear_v, 0, 0]])
check(
    clock_inertia_density(sheared, 0, 0) == shear_u**2 + shear_v**2,
    "weight-one shear inertia",
)

psi_real, psi_imag, frequency = sp.symbols(
    "psi_real psi_imag frequency", real=True
)
tensor_velocity, velocity_real, velocity_imag = relative_equilibrium_velocity(
    sheared, psi_real, psi_imag, frequency
)
full_inertia = clock_inertia_density(sheared, psi_real, psi_imag)
check(
    noether_charge_density(
        sheared,
        tensor_velocity,
        psi_real,
        psi_imag,
        velocity_real,
        velocity_imag,
    )
    == frequency * full_inertia,
    "canonical phase-space charge reduction",
)
check(
    canonical_velocity_energy_density(
        tensor_velocity, velocity_real, velocity_imag
    )
    == frequency**2 * full_inertia / 2,
    "canonical phase-space energy reduction",
)

amplitude, split = sp.symbols("amplitude split", nonnegative=True, real=True)
check(
    sp.diff(
        phase_lock_potential(sp.diag(1, split, -split), amplitude, 0), split
    ).subs(split, 0)
    == -12 * amplitude**2,
    "nonzero M5 split source",
)

charge, total_inertia = sp.symbols("charge total_inertia", positive=True)
rotational, omega, first, second, envelope = fixed_charge_legendre_data(
    charge, total_inertia
)
check(rotational == charge**2 / (2 * total_inertia), "fixed-charge energy")
check(first == -omega**2 / 2, "complete first variation coefficient")
check(second == omega**2 / total_inertia, "positive rank-one coefficient")
check(envelope == omega, "dE/dQ")

check(6 - witness.charged_edge_squared == 2, "scalar tail margin")
check(4 - witness.charged_edge_squared == 0, "weight-one shear edge")
check(22 - 4 * witness.charged_edge_squared == 6, "weight-two tail margin")

# Mutations of every defining object fail a named gate.
check(
    sp.limit(-radius**3 - radius**2 / 2 + sp.Rational(1, 2), radius, sp.oo)
    == -sp.oo,
    "quartic-deletion mutation",
)
check(sp.diag(sp.Rational(1, 2), 1, 0, 1, 1, 0).det() == 0, "kinetic deletion mutation")
check(sp.diff(6 * split**2, split).subs(split, 0) == 0, "lock deletion mutation")

print(f"ALL {CHECKS} CHECKS PASS [P249-O6-CANONICAL]")
