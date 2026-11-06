"""All-slot barriers and same-scalar physical joint mean evaluations."""
from dataclasses import dataclass, fields
from contextlib import contextmanager

import numpy as np

from .material_joint_field_scalar import coordinate_correction
from .material_physical_operators import PhysicalFrozenAction, PhysicalFieldVertices
from .qe_material_hamiltonian import MaterialHamiltonian
from .streamed_material_evolution import real_physical


@dataclass(frozen=True)
class JointStage:
    density: np.ndarray
    gradient: np.ndarray
    field_gradient_electronic: np.ndarray
    potential: np.ndarray
    u: np.ndarray
    p: np.ndarray
    coordinate: np.ndarray
    momentum: np.ndarray
    shift: np.ndarray
    velocity: np.ndarray
    force: np.ndarray
    field_velocity: np.ndarray
    field_force: np.ndarray
    electronic_current: np.ndarray
    ionic_current: np.ndarray
    energy: float
    number: float

    def save(self, directory, budget, *, metadata):
        directory.mkdir(parents=True, exist_ok=False)
        arrays = {f.name: budget.array(directory/(f.name+'.npy'), getattr(self, f.name))
                  for f in fields(self) if isinstance(getattr(self, f.name), np.ndarray)}
        return dict(metadata, arrays=arrays, energy_total_ry=self.energy, number_per_cell=self.number)

    @classmethod
    def restore(cls, record):
        from .streamed_material_state import authenticate
        values = {name: np.load(authenticate(binding), allow_pickle=False)
                  for name, binding in record['arrays'].items()}
        return cls(**values, energy=record['energy_total_ry'], number=record['number_per_cell'])


class JointReduction:
    """Factory yields physical coupling, H0 object (or action), field-only vertex.

    Material H0/g MUST share factors; the generic action interface exists for
    independent finite oracles, not legacy material scalar adoption.
    """
    def __init__(self, reference, factory, coefficients, field, charges, *, model_identity,
                 expected_slots, grid_measure, budget, exponential, reality_tolerance=1e-10):
        self.reference, self.factory = reference, factory
        self.coefficients, self.field = coefficients, field
        self.charges = np.asarray(charges, float)
        self.identity, self.expected_slots = model_identity, tuple(expected_slots)
        self.measure, self.budget = float(grid_measure), budget
        self.exponential, self.reality_tolerance = exponential, reality_tolerance
        self.reductions = self.exponentials = 0
        if (np.any(reference.held_density) or np.any(reference.held_gradient)
                or coefficients.mixed.shape != (reference.ndof, 3)
                or self.charges.shape != (reference.ndof//3,) or np.any(self.charges <= 0)
                or not np.isfinite(self.charges).all() or not 0 < reality_tolerance < 1
                or not np.isfinite(self.measure) or self.measure <= 0 or not model_identity
                or not self.expected_slots or len(set(self.expected_slots)) != len(self.expected_slots)):
            raise ValueError('complete joint physical scalar with no held electrons required')

    def real(self, value, name):
        return real_physical(value, name, self.reality_tolerance)

    def validate(self, state):
        r = self.reference
        if (state.identity != self.identity or tuple(s.identity for s in state.slots) != self.expected_slots
                or state.u.shape != r.masses.shape or state.coordinate.shape != (self.field.size,)):
            raise ValueError('joint model/coverage/coordinates differ')
        for component in range(r.nbar.shape[1]):
            if sum(s.component == component for s in state.slots) != r.ncells:
                raise ValueError('every original independent spin/grid slot required')
        if any(s.weight != 1/r.ncells or not 0 <= s.component < r.nbar.shape[1] for s in state.slots):
            raise ValueError('original uniform explicit slot weights required')

    @contextmanager
    def operators(self, slot):
        self.budget.check()
        with self.factory(slot) as (coupling, h0, vertex):
            if (coupling.npw != slot.shape[0] or coupling.nr != self.reference.nbar.shape[0]
                    or coupling.ndof != self.reference.ndof or coupling.grid_measure_bohr3 != self.measure):
                raise ValueError('physical operator and original slot/grid differ')
            if isinstance(h0, MaterialHamiltonian):
                if coupling.factors is not h0.factors or not isinstance(vertex, PhysicalFieldVertices) or vertex.hamiltonian is not h0:
                    raise ValueError('H0/g/L must share the identical physical factors')
                if not np.isclose(self.field.volume_bohr3,
                        self.reference.ncells*coupling.cell_volume_bohr3, rtol=1e-12, atol=0):
                    raise ValueError('Maxwell volume differs from original matter volume')
            yield coupling, h0, vertex
        self.budget.check()

    def counted(self, action):
        def apply(vectors):
            self.budget.action(vectors)
            return action(vectors)
        return apply

    def frozen_action(self, coupling, h0, vertex, slot, stage):
        if isinstance(h0, MaterialHamiltonian):
            action = PhysicalFrozenAction(h0, coupling, potential=stage.potential[:, slot.component],
                                          direction=stage.u, shift=stage.shift).apply
        else:
            optical = vertex.prepare(stage.shift)
            def action(v):
                return (h0(v)+coupling.potential_action(stage.potential[:, slot.component], v)
                        +coupling.ionic_action(stage.u, v)+optical.apply(v))
        return self.counted(action)

    def propagate(self, coupling, h0, vertex, slot, stage, psi, interval):
        self.exponentials += 1
        return self.exponential(self.frozen_action(coupling, h0, vertex, slot, stage), psi, interval)

    def stage(self, density, gradient, ell, one_body, number, u, p, coordinate, momentum):
        r = self.reference
        density, gradient, ell = (self.real(v, name) for v, name in
            ((density, 'density'), (gradient, 'ionic gradient'), (ell, 'field gradient')))
        potential = self.real(r.kernel(np.asarray(density-r.nbar, complex)), 'K delta n')
        if potential.shape != r.nbar.shape:
            raise ValueError('kernel grid/spin shape differs')
        shift = self.field.cartesian_shift(coordinate)
        # Physical H0 and g already contain DeltaV0/Delta_g. DeltaC occurs ONLY here.
        correction = coordinate_correction(self.coefficients, u, p, shift, r.masses,
                                            self.charges, np.zeros(r.ndof), ell)
        base = (one_body+self.measure*np.sum((density-r.nbar)*potential)/2
                +u @ (gradient-r.ionic_force-r.support_force)+u @ r.curvature @ u/2
                +np.sum(p*p/r.masses)/2)
        energy = float(self.real(r.ncells*(base+correction['energy'])+self.field.energy(momentum), 'energy'))
        force = self.real(r.ionic_force+r.support_force-gradient-r.curvature @ u
                          +correction['force_correction'], 'force')
        field_force = self.real(-r.ncells*self.field.polarizations @ correction['field_gradient'], 'Maxwell force')
        return JointStage(density, gradient, ell, potential, np.array(u), np.array(p),
            np.array(coordinate), np.array(momentum), shift, correction['velocity'], force,
            momentum/self.field.mass, field_force,
            self.real(r.ncells*correction['electronic_current_hartree'], 'electronic current'),
            self.real(r.ncells*correction['ionic_current_hartree'], 'ionic current'),
            energy, float(self.real(number, 'number')))

    def reduce(self, state, *, frozen_stage=None, interval=0., coordinates=None):
        self.validate(state)
        r = self.reference
        density = np.zeros_like(r.nbar, dtype=complex)
        gradient, ell = np.zeros(r.ndof, complex), np.zeros(3, complex)
        one_body, number = 0j, 0j
        for slot in state.slots:
            psi, d = state.arrays(slot)
            with self.operators(slot) as (coupling, h0, vertex):
                if frozen_stage is not None:
                    psi = self.propagate(coupling, h0, vertex, slot, frozen_stage, psi, interval)
                n, g = (coupling.prepared_moments(psi, d) if hasattr(coupling, 'prepared_moments')
                        else coupling.moments(psi, psi, d))
                l = vertex.moments(psi, d)['field_gradient']
                bare = PhysicalFrozenAction(h0).apply if isinstance(h0, MaterialHamiltonian) else h0
                e = np.einsum('ij,ji->', d, psi.conj().T @ self.counted(bare)(psi))
                density[:, slot.component] += slot.weight*n
                gradient += slot.weight*g
                ell += slot.weight*l
                one_body += slot.weight*e
                number += slot.weight*self.measure*np.sum(n)
        self.reductions += 1
        # No kernel, ionic or Maxwell update can observe a partial reduction.
        coords = (state.u, state.p, state.coordinate, state.momentum) if coordinates is None else coordinates
        return self.stage(density, gradient, ell, one_body, number, *coords)
