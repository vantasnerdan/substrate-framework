"""Couple the full electronic tangent to independent ionic inertia and curvature.

Ry time, physical +g*u electronic perturbation, negative gradient ionic force.
No damping, screened-vertex substitution, counterterm, or equilibrium-noise model.
"""
from __future__ import annotations

import numpy as np

from .pw_electronic_dynamics import PWElectronicDynamics


class PWCoupledDynamics:
    """State is [six Cartesian u, six p, full electronic spectral amplitudes]."""

    def __init__(self, electronic: PWElectronicDynamics, kernel, frozen_curvature, masses):
        self.electronic = electronic
        self.kernel = kernel
        self.frozen = np.asarray(frozen_curvature, dtype=complex)
        self.masses = np.asarray(masses, dtype=float)
        if (self.frozen.shape != (6, 6) or not np.all(np.isfinite(self.frozen))
                or self.masses.shape != (6,) or not np.all(np.isfinite(self.masses))
                or np.any(self.masses <= 0)):
            raise ValueError("finite independent 6x6 frozen curvature and positive six-coordinate masses required")
        self.state_size = 12 + electronic.state_size

    def zero_state(self):
        return np.zeros(self.state_size, dtype=complex)

    def derivative(self, state):
        state = np.asarray(state, dtype=complex)
        if state.shape != (self.state_size,):
            raise ValueError("coupled state shape differs")
        u, p, transition = state[:6], state[6:12], state[12:]
        observable = self.electronic.observables(transition)
        potential = self.kernel(observable.density)
        result = np.empty_like(state)
        result[:6] = p / self.masses
        result[6:12] = -self.frozen @ u - observable.ionic_gradient
        result[12:] = self.electronic.derivative(transition, potential, u)
        return result

    def derivative_real(self, tau, state):
        """Interleaved real/imaginary doubles; no complex-state information discarded."""
        return self.derivative(np.ascontiguousarray(state).view(np.complex128)).view(np.float64)

    def absolute_tolerances(self, *, velocity=1e-12, displacement=1e-12, electronic=1e-10):
        """Occupation-weighted absolute scales, not a cutoff or an error guarantee.

        SciPy uses its standard RMS error norm. Spectral amplitudes are unweighted;
        electronic atol is a scale for sqrt(weighted occupation)*amplitude. Exact
        zero-weight sectors keep all state coordinates with unrestrictive finite
        atol. Relative tolerance is separately supplied by the caller.
        """
        tolerance = np.empty(self.state_size)
        tolerance[:6] = displacement
        tolerance[6:12] = self.masses * velocity
        for block in self.electronic.state_layout:
            occupation = self.electronic.model.channels[block.channel_index].weighted_occupations
            scale = np.sqrt(occupation)
            value = np.divide(electronic, scale, out=np.full_like(scale, 1e200), where=scale > 0)
            tolerance[12 + block.start:12 + block.stop] = np.broadcast_to(value, block.shape).ravel()
        if not np.all(np.isfinite(tolerance)) or np.any(tolerance <= 0):
            raise ValueError("positive finite component tolerances required")
        return np.repeat(tolerance, 2)

    def mechanical_observables(self, state):
        """Raw per-q ionic energy/power only; NOT conserved total coupled energy."""
        u, p = state[:6], state[6:12]
        value = self.electronic.observables(state[12:])
        velocity = p / self.masses
        return value, {
            "ionic_kinetic_ry": float(0.5 * np.vdot(p, velocity).real),
            "ionic_frozen_quadratic_ry": complex(0.5 * np.vdot(u, self.frozen @ u)),
            "electronic_power_into_ions_ry_per_tau": float(-np.vdot(velocity, value.ionic_gradient).real),
            "force_norm_ry_per_bohr": float(np.linalg.norm(-self.frozen @ u - value.ionic_gradient)),
        }
