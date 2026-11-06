"""Streamed full-PW evolution of the PreparedReference effective scalar.

Time is inverse Ry (hbar=1). Fixed original complex D is never diagonalized.
Self-consistent exponential midpoint is second order when its nonlinear solve
converges; it is not an Euler step or an exact energy-conserving integrator.
Each stage completes ALL sector reductions before any updated stage is used.
"""
from dataclasses import dataclass
from contextlib import contextmanager
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import LinearOperator, expm_multiply

from .prepared_reference import PreparedReference


@dataclass(frozen=True)
class Stage:
    density: np.ndarray
    gradient: np.ndarray
    potential: np.ndarray
    u: np.ndarray
    force: np.ndarray
    energy: float | None
    number: complex


def real_physical(value, name, tolerance):
    """Accept measured imaginary roundoff only, never repair an operator."""
    value = np.asarray(value)
    scale = 1 + np.max(np.abs(value), initial=0.)
    residual = np.max(np.abs(value.imag), initial=0.)
    if not np.isfinite(value).all() or residual > tolerance * scale:
        raise ValueError(name + ' is nonfinite or has non-roundoff imaginary part')
    return value.real


def exponential_action(action, frame, interval):
    """SciPy scaling/Taylor exponential on all PW rows, with exact adjoint.

    traceA=0 is a legal arbitrary centering choice, NOT a claim Tr(H)=0.
    Hermiticity is the supplied physical action contract; no symmetrization.
    SciPy's machine-precision action target is not a material error certificate.
    """
    n = frame.shape[0]
    def forward(v):
        one = v.ndim == 1
        result = -1j * interval * action(v[:, None] if one else v)
        return result[:, 0] if one else result
    def adjoint(v):
        one = v.ndim == 1
        result = 1j * interval * action(v[:, None] if one else v)
        return result[:, 0] if one else result
    operator = LinearOperator((n, n), matvec=forward, matmat=forward,
                              rmatvec=adjoint, rmatmat=adjoint, dtype=np.complex128)
    return expm_multiply(operator, frame, traceA=0.)


class StreamedMeanEvolution:
    """Factory(slot) is a context manager yielding (DyadCoupling, actual H0).

    The factory must authenticate the slot basis/gauge, rebuild H0 at that
    actual slot, and release per-sector operators on exit. No slot frames or
    full-PW Hamiltonian matrices are cached by this consumer.
    """
    def __init__(self, reference, factory, *, model_identity, expected_slots,
                 grid_measure, reality_tolerance=1e-10, max_h0_actions=None):
        if not isinstance(reference, PreparedReference):
            raise TypeError('existing PreparedReference required')
        self.reference, self.factory = reference, factory
        self.identity = model_identity
        self.expected_slots = tuple(expected_slots)
        self.measure = float(grid_measure)
        self.reality_tolerance = float(reality_tolerance)
        if (not model_identity or not self.expected_slots
                or len(set(self.expected_slots)) != len(self.expected_slots)
                or self.measure <= 0 or not np.isfinite(self.measure)
                or not 0 < self.reality_tolerance < 1):
            raise ValueError('explicit complete coverage, measure and tolerance required')
        if np.any(reference.held_density) or np.any(reference.held_gradient):
            raise ValueError('full-material streamed state admits no held electrons')
        if max_h0_actions is not None and (type(max_h0_actions) is not int or max_h0_actions < 1):
            raise ValueError('positive explicit H0 action cap required')
        self.max_h0_actions = max_h0_actions
        self.counts = dict(operator_constructions=0, h0_calls=0, h0_vector_columns=0,
                           moment_reductions=0, kernel_calls=0)

    def validate(self, state):
        if (state.identity != self.identity
                or tuple(s.identity for s in state.slots) != self.expected_slots
                or state.u.shape != self.reference.masses.shape):
            raise ValueError('state model, ordered coverage or ionic coordinates differ')
        for s in state.slots:
            if s.weight != 1 / self.reference.ncells or s.component >= self.reference.nbar.shape[1]:
                raise ValueError('explicit independent spin/full-grid weights required')
        for component in range(self.reference.nbar.shape[1]):
            if sum(s.component == component for s in state.slots) != self.reference.ncells:
                raise ValueError('each spin component must have every uniform mesh slot')

    @contextmanager
    def operators(self, slot):
        with self.factory(slot) as (coupling, h0):
            if (coupling.npw != slot.shape[0] or coupling.nr != self.reference.nbar.shape[0]
                    or coupling.ndof != self.reference.ndof
                    or coupling.grid_measure_bohr3 != self.measure):
                raise ValueError('actual operator differs from complete slot/grid coordinates')
            self.counts['operator_constructions'] += 1
            def counted_h0(vectors):
                if self.max_h0_actions is not None and self.counts['h0_calls'] >= self.max_h0_actions:
                    raise RuntimeError('full-PW H0 action cap exhausted; no successor published')
                self.counts['h0_calls'] += 1
                self.counts['h0_vector_columns'] += 1 if vectors.ndim == 1 else vectors.shape[1]
                return h0(vectors)
            yield coupling, counted_h0

    def force(self, gradient, u):
        r = self.reference
        return r.ionic_force + r.support_force - gradient - r.curvature @ u

    def stage(self, density, gradient, u, *, one_body=None, p=None, number=0j):
        r = self.reference
        density = real_physical(density, 'global density', self.reality_tolerance)
        gradient = real_physical(gradient, 'positive electronic gradient', self.reality_tolerance)
        potential = real_physical(r.kernel(density - r.nbar), 'K delta n', self.reality_tolerance)
        if potential.shape != r.nbar.shape:
            raise ValueError('kernel returned wrong spin/grid shape')
        energy = None
        if one_body is not None:
            raw = (one_body + self.measure * np.sum((density-r.nbar)*potential)/2
                   + u @ (gradient-r.ionic_force-r.support_force)
                   + u @ r.curvature @ u/2 + np.sum(p*p/r.masses)/2)
            energy = float(real_physical(raw, 'effective scalar', self.reality_tolerance))
        self.counts['kernel_calls'] += 1
        return Stage(density, gradient, potential, np.array(u, copy=True),
                     self.force(gradient, u), energy, complex(number))

    @staticmethod
    def action(coupling, h0, slot, stage, vectors):
        return (h0(vectors) + coupling.potential_action(stage.potential[:, slot.component], vectors)
                + coupling.ionic_action(stage.u, vectors))

    def reduce(self, state, *, frozen_stage=None, interval=0., u=None, energy=False):
        self.validate(state)
        r = self.reference
        density, gradient = np.zeros_like(r.nbar, dtype=complex), np.zeros(r.ndof, complex)
        one_body, number = 0j, 0j
        for slot in state.slots:
            psi, d = state.arrays(slot)
            with self.operators(slot) as (coupling, h0):
                if frozen_stage is not None:
                    psi = exponential_action(lambda v: self.action(coupling, h0, slot, frozen_stage, v), psi, interval)
                if hasattr(coupling, 'prepared_moments'):
                    n, g = coupling.prepared_moments(psi, d)
                else:
                    n, g = coupling.moments(psi, psi, d)
                self.counts['moment_reductions'] += 1
                density[:, slot.component] += slot.weight*n
                gradient += slot.weight*g
                # Number via physical density avoids an unnecessary wide Gram.
                number += slot.weight*self.measure*np.sum(n)
                if energy:
                    one_body += slot.weight*np.einsum('ij,ji->', d, psi.conj().T @ h0(psi))
            del psi, d
        return self.stage(density, gradient, state.u if u is None else u,
                          one_body=one_body if energy else None, p=state.p, number=number)

    def derivative(self, state, directory):
        """Useful exact initial RHS; explicitly NOT a time trajectory."""
        stage = self.reduce(state, energy=True)
        from .streamed_material_state import save_array
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=False)
        bindings = []
        for index, slot in enumerate(state.slots):
            psi, _ = state.arrays(slot)
            with self.operators(slot) as (coupling, h0):
                value = -1j*self.action(coupling, h0, slot, stage, psi)
                bindings.append(dict(slot=slot.metadata(), derivative=save_array(
                    directory / ('dpsi-%04d.npy' % index), value)))
        return stage, dict(scope='initial derivative only; no finite-time/coherence result',
                           dpsi=bindings, du=(state.p/self.reference.masses).tolist(),
                           dp=stage.force.tolist())

    def step(self, state, directory, *, dt, max_iterations, tolerance):
        if (not np.isfinite(dt) or dt <= 0 or type(max_iterations) is not int
                or max_iterations < 1 or not 0 < tolerance < 1):
            raise ValueError('positive time, bounded iterations and explicit tolerance required')
        r = self.reference
        stage = self.reduce(state)
        residuals = []
        for iteration in range(max_iterations):
            target_u = state.u + dt*state.p/(2*r.masses) + dt*dt*stage.force/(4*r.masses)
            updated = self.reduce(state, frozen_stage=stage, interval=dt/2, u=target_u)
            residual = max(float(np.max(np.abs(a-b)/(1+np.abs(a)))) for a, b in (
                (updated.potential, stage.potential), (updated.gradient, stage.gradient),
                (updated.u, stage.u)))
            residuals.append(residual)
            if residual <= tolerance:
                break
            stage = updated
        else:
            raise RuntimeError('global midpoint failed to converge; no successor published: '+str(residuals))
        # Use the tested frozen stage, not a fresh unsolved mixture. Its update
        # residual bounds the nonlinear midpoint defect; time refinement still
        # must dominate that defect (typically tolerance << dt^3).
        p = state.p + dt*stage.force
        u = state.u + dt*(state.p+p)/(2*r.masses)
        def frames():
            for slot in state.slots:
                psi, _ = state.arrays(slot)
                with self.operators(slot) as (coupling, h0):
                    yield exponential_action(lambda v: self.action(coupling, h0, slot, stage, v), psi, dt)
            evidence['work_counts'] = dict(self.counts)
        evidence = dict(method='globally self-consistent exponential midpoint',
                        midpoint_iterations=iteration+1, midpoint_residuals=residuals,
                        midpoint_tolerance=tolerance, dt_ry_inverse=dt,
                        prior_boundary=str(state.path),
                        exponential='scipy.expm_multiply full-PW action; machine precision target, not independently certified',
                        scope='bounded mean evolution, not coherence or equilibrium evidence')
        return state.successor(directory, frames(), time=state.time+dt, u=u, p=p, evidence=evidence)

    def evolve(self, state, directory, *, physical_interval, steps, max_steps,
               max_iterations, tolerance, refinement_label):
        if (type(steps) is not int or type(max_steps) is not int or not 1 <= steps <= max_steps
                or not np.isfinite(physical_interval) or physical_interval <= 0 or not refinement_label):
            raise ValueError('bounded steps, independent positive physical interval and refinement label required')
        records = []
        for index in range(steps):
            state = self.step(state, Path(directory)/('step-%06d' % (index+1)),
                              dt=physical_interval/steps, max_iterations=max_iterations, tolerance=tolerance)
            stage = self.reduce(state, energy=True)
            records.append(dict(time_ry_inverse=state.time, energy_ry=stage.energy,
                                number=[stage.number.real, stage.number.imag],
                                boundary=str(state.path), refinement_label=refinement_label))
        return state, records
