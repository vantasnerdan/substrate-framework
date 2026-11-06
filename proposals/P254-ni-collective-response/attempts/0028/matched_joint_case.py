"""Case-specific physical joint construction from completed matched coefficients.

No coefficient recalculation: these archives ALREADY contain physical onsite
curvature. DeltaC is exactly zero, not the old loaded A25 correction. The fixed
geometric approximation and nonstationary initial density remain unchanged.
"""
from contextlib import contextmanager
from dataclasses import dataclass
from hashlib import sha256
import json

import matched_reference_io as io
import numpy as np
from substrate_framework.material_control_reference import bound_json, bound_arrays
from substrate_framework.material_joint_field_scalar import JointFieldCoefficients, conditional_field_gaussian
from substrate_framework.ionic_gaussian_preparation import prepare_local_trap_vacuum
from substrate_framework.material_modes import AMU_RY
from substrate_framework.physical_field_environment import UniformMaxwellField
from substrate_framework.prepared_reference import PreparedReference
from substrate_framework.streamed_material_evolution import real_physical


def slot_identity(slot):
    return str(slot.ispin)+':'+','.join(map(str, slot.mesh_index))


def model_identity(specification):
    return sha256(json.dumps(specification, sort_keys=True).encode()).hexdigest()


def joint_blocks(arrays, *, authority, reality_tolerance=1e-10):
    """Only the new join: physical common C plus ZERO additional spatial delta."""
    masses = np.repeat(np.asarray(arrays['atomic_masses_amu'])*AMU_RY, 3)
    charges = np.asarray(arrays['valence_charges'])
    d = len(masses)
    coefficients = JointFieldCoefficients(np.zeros((d, d)), arrays['contact'],
        np.asarray(arrays['mixed_ionic_field']).reshape(d, 3), authority)
    field = UniformMaxwellField(64*float(abs(np.linalg.det(arrays['cell_bohr']))),
                                np.eye(3), authority+'; original BvK4^3, not cavity')
    ionic = prepare_local_trap_vacuum(arrays['common_force_lattice_ry_bohr2'],
        arrays['atomic_masses_amu'], state_authority=authority,
        reality_relative_tolerance=reality_tolerance)
    photon = conditional_field_gaussian(field, coefficients.contact, masses, charges, 64,
        state_authority=authority, reality_relative_tolerance=reality_tolerance)
    return coefficients, field, ionic, photon, masses, charges


@dataclass
class MatchedJointCase:
    case: str
    control: object
    background: object
    ledger: dict
    rows: tuple
    arrays: dict
    aggregate: dict
    aggregate_binding: dict
    source_states: dict
    source_metadata: dict
    specification: dict
    identity: str
    reference: object
    coefficients: object
    field: object
    ionic_preparation: object
    field_preparation: dict
    charges: np.ndarray
    reality_tolerance: float

    @property
    def grid(self):
        return self.control.grid

    @property
    def ordered(self):
        return tuple(row[1] for row in self.rows)

    @property
    def expected_slots(self):
        return tuple(slot_identity(s) for s in self.ordered)

    def make_mean(self, budget, *, block_size=8):
        from substrate_framework.material_physical_operators import physical_hamiltonian, physical_coupling, PhysicalFieldVertices
        from substrate_framework.streamed_joint_mean import StreamedJointMean
        slots = {slot_identity(s): s for s in self.ordered}
        @contextmanager
        def factory(record):
            slot = slots[record.identity]
            frame = self.grid.load_frame(slot)
            if (frame.basis.identity != record.basis_identity or frame.gauge_identity != record.gauge_identity
                    or frame.coefficients.shape != record.shape or slot.source_binding['sha256'] != record.source_identity):
                raise ValueError('original raw frame/basis/gauge/source differs')
            basis = frame.basis
            del frame
            h0 = physical_hamiltonian(self.background, basis, slot.ispin, block_size)
            yield physical_coupling(self.background.inputs, basis, h0.factors, block_size), h0, PhysicalFieldVertices(h0)
        ref = self.control.reference
        return StreamedJointMean(self.reference, factory, self.coefficients, self.field, self.charges,
            model_identity=self.identity, expected_slots=self.expected_slots,
            grid_measure=abs(np.linalg.det(ref.cell_bohr))/int(np.prod(ref.fft_smooth)),
            budget=budget, reality_tolerance=self.reality_tolerance)

    def make_history(self, final_binding, budget, *, block_size=8):
        from substrate_framework.material_joint_response import MaterialJointHistory
        return MaterialJointHistory(final_binding, model_identity=self.identity,
            background=self.background, grid=self.grid, ordered=self.ordered,
            slot_identity=slot_identity, budget=budget, block_size=block_size)

    def bind_response(self, history, q):
        q = np.asarray(q, float)
        step = np.rint(4*q).astype(int)
        if q.shape != (3,) or not np.allclose(4*q, step, rtol=0, atol=1e-10):
            raise ValueError('original commensurate4^3 matter q required')
        curvature = self.arrays['common_q_grid_ry_bohr2'][tuple(step % 4)]
        return history.bind_response(q=q, coefficients=self.coefficients,
            common_curvature=curvature.reshape(len(self.reference.masses), -1),
            masses=self.reference.masses, charges=self.charges,
            polarizations=self.field.polarizations, field_mass=self.field.mass, ncells=64,
            density_shape=self.reference.nbar.shape, scalar_identity=self.identity,
            coverage_identity=self.specification['context']['population_ledger']['sha256'])


def load_case(case, *, reality_tolerance=1e-10):
    from matched_joint_sources import registration, source_catalog, validate_aggregate, validate_atoms
    reg_binding, reg, implementation = registration()
    if case not in reg['cases']:
        raise ValueError('one registered original case required')
    _, original_registration = io.registration()
    control, background, ledger, rows, context = io.load_case(case, original_registration)
    aggregate_binding = reg['cases'][case]
    aggregate = bound_json(aggregate_binding)
    rows = tuple(sorted(rows, key=lambda row: (row[1].ispin, row[1].mesh_index)))
    validate_aggregate(aggregate, context, rows)
    arrays = bound_arrays(aggregate['archive'])
    validate_atoms(control, background, arrays)
    states, metadata = source_catalog(aggregate, context, rows, ledger)
    # Group only storage order, never weights or physical slot identities.
    rows = tuple(sorted(rows, key=lambda row: (row[1].ispin, row[1].source_ik, row[1].mesh_index)))
    spec = dict(case=case, registration=reg_binding, implementation=implementation,
        context=context, aggregate=aggregate_binding,
        source_states={str(k): v for k, v in states.items()},
        scalar=reg['construction']['scalar'], preparation=reg['construction']['preparation'],
        reality_tolerance=reality_tolerance, initial_means='u=p=q0=pi0=0; explicit zero support',
        approximation=aggregate['approximation'], coefficient_defect_upper_bound=None,
        scientific_limits=aggregate['scientific_limits'], positive_time_symmetry=False)
    identity = model_identity(spec)
    coefficients, field, ionic, photon, masses, charges = joint_blocks(arrays,
        authority=identity, reality_tolerance=reality_tolerance)
    d = len(masses)
    real = lambda value, name: real_physical(value, name, reality_tolerance)
    reference = PreparedReference(real(arrays['saved_density_bohr3'], 'saved nbar'),
        real(arrays['ionic_force_ry_bohr'], 'bare ionic force').reshape(d),
        np.zeros_like(arrays['saved_density_bohr3'], dtype=float), np.zeros(d),
        real(arrays['common_q_grid_ry_bohr2'][0, 0, 0], 'physical common Gamma C').reshape(d, d),
        masses, np.zeros(d), background.apply_kernel, 64, aggregate_binding['sha256'],
        aggregate_binding['sha256']+':finite ionic Taylor prescription', context['population_ledger']['sha256'],
        'explicit zero support; original nonstationarity retained',
        control.reference.source['sha256']+':original-atom-order', context['preparation']['sha256'])
    return MatchedJointCase(case, control, background, ledger, rows, arrays, aggregate,
        aggregate_binding, states, metadata, spec, identity, reference, coefficients,
        field, ionic, photon, charges, reality_tolerance)
