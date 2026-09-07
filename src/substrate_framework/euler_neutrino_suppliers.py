"""Exact supplier ledgers for P253 neutral-current and flavor routes."""

from __future__ import annotations

import sympy as sp


def two_flavor_transition_probability(theta, delta_frequency, time):
    """Fixed-momentum two-state transition probability."""

    return sp.sin(2 * theta) ** 2 * sp.sin(delta_frequency * time / 2) ** 2


def packet_transition_probability(theta, coherence):
    """Momentum-averaged probability in terms of C(t)=<exp(-i Delta nu t)>."""

    return sp.sin(2 * theta) ** 2 * (1 - sp.re(coherence)) / 2


def relativistic_phase_leading(delta_mass_squared, speed, distance, energy, action):
    """Leading high-energy phase Delta(m^2)c^3 L/(2 E S)."""

    return delta_mass_squared * speed**3 * distance / (2 * energy * action)


def material_frame_hamiltonian(material_derivative_U, U_adjoint):
    """Pure-frame Hermitian generator H_U=i(D_t U)U^*."""

    return sp.I * material_derivative_U * U_adjoint


def helicity_flux(velocity, vorticity, pressure, speed_squared):
    """Return h and coefficient form of h*u+(p-|u|^2/2)*omega.

    Vector inputs may be SymPy matrices. The returned tuple keeps the exact
    pressure coefficient visible for independent convention checks.
    """

    h = (velocity.dot(vorticity) if hasattr(velocity, "dot") else velocity * vorticity)
    return h, h * velocity + (pressure - speed_squared / 2) * vorticity


def current_parity(tag_parity: str) -> dict[str, str]:
    """O(3) parity of (tag, tag*u) for polar velocity u."""

    if tag_parity == "scalar":
        return {"density": "scalar", "spatial_current": "polar"}
    if tag_parity == "pseudoscalar":
        return {"density": "pseudoscalar", "spatial_current": "axial"}
    raise ValueError("tag_parity must be 'scalar' or 'pseudoscalar'")
