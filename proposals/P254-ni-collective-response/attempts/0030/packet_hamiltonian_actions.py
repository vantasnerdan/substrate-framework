"""Original-packet physical actions, with no initialization or action-bank replay.

``restore_packet_hamiltonian(packet_path)`` restores one authenticated original
sector. Both context methods return ``(ActionFactors, work_dict)``. Arrays are
unweighted, full-PW by all original columns, with independently applied Hilbert
adjoints. H means the literal saved t0 matched Hamiltonian, including its saved
K(delta n) potential, ionic coordinate and Cartesian field shift; it is not a
saved eigenvalue diagonal and makes no stationarity assertion. Ionic force is
minus the complex-linear physical ionic Hamiltonian derivative.

Only saved background/control, native basis, small state coordinates and the
already computed initial stage are restored. No Gram, population, kernel or
trajectory is recomputed. Actions stay in memory and are not cached by context
or written to disk. Callers own resource limits and release returned actions.
"""
from dataclasses import dataclass, field
from pathlib import Path
import resource
from time import perf_counter, process_time

import numpy as np

from full_noise_probe_inputs import (
    HERE, _identity, _same_binding, _read_json, _force_actions,
    material_io, slot_identity, validate_aggregate, bound_arrays, bound_json,
    ActionFactors, transport_real_frame, load_wavefunctions, DiskMean,
    authenticate, file_binding,
)
from substrate_framework.material_physical_operators import (
    build_physical_factors, physical_coupling, PhysicalFrozenAction,
)
from substrate_framework.qe_material_hamiltonian import MaterialHamiltonian
from substrate_framework.streamed_joint_state import DiskJointMean


_ORIGINAL_PACKET_SHA256 = '6639122665f0e31b77b2a04a996bb452b874818ff51955983af37e61c23b46b0'
# Existing immutable full initial reduction: only its t0 arrays are consumed.
_SAVED_MEAN = dict(
    path=str(HERE.parent / '0028' / 'matched-working-precision-mean-complete-001.json'),
    bytes=10063,
    sha256='3048fc8e630facf9403a8f4d58a5cf2d9b960416a736df8e4b7b235a12966b75',
)


def _load_array(binding):
    value = np.load(authenticate(binding), mmap_mode='r', allow_pickle=False)
    if not np.isfinite(value).all():
        raise ValueError('finite authenticated saved array required')
    return value


def _direction_spec(direction):
    # JSON has no complex scalar type; neither leg may discard the imaginary part.
    return dict(real=direction.real.tolist(), imaginary=direction.imag.tolist(),
                layout='original atom-major nat by Cartesian xyz',
                operator_convention='sum_i direction_i F_i, without conjugation')


def _no_replay_work():
    return dict(native_calls=0, trajectory_actions=0, initial_kernel_reductions=0,
        gram_recomputations=0, population_refills=0, preparation_replays=0,
        h37_action_replays=0, action_disk_bytes=0, quadrature_applied=False,
        original_inputs_unchanged=True)


def _hamiltonian_actions(prepared, psi, block_size):
    # Unlike _force_actions there is deliberately NO force minus sign here.
    action = np.empty(psi.shape, dtype=np.complex128)
    adjoint = np.empty(psi.shape, dtype=np.complex128)
    forward_seconds = adjoint_seconds = 0.
    blocks = 0
    for start in range(0, psi.shape[1], block_size):
        stop = min(start + block_size, psi.shape[1])
        t = perf_counter()
        action[:, start:stop] = prepared.apply(psi[:, start:stop])
        forward_seconds += perf_counter() - t
        t = perf_counter()
        adjoint[:, start:stop] = prepared.adjoint_apply(psi[:, start:stop])
        adjoint_seconds += perf_counter() - t
        blocks += 1
    action.setflags(write=False)
    adjoint.setflags(write=False)
    return action, adjoint, dict(h0_forward_actions=1, h0_adjoint_actions=1,
        h0_actions=2, force_forward_actions=0, force_adjoint_actions=0,
        forward_block_calls=blocks, adjoint_block_calls=blocks,
        forward_action_columns=psi.shape[1], adjoint_action_columns=psi.shape[1],
        forward_seconds=forward_seconds, adjoint_seconds=adjoint_seconds,
        block_size=block_size)


@dataclass(frozen=True)
class PacketHamiltonianContext:
    """One original sector; shared physical factories, no retained action arrays.

    ``operator_work`` describes restoration ONLY. Each method's returned work
    describes that invocation ONLY, with the exact operator specification and
    shared restoration work separately identified (not charged twice).
    """
    psi: np.ndarray
    gram: np.ndarray
    density_coefficients: np.ndarray
    active_coordinates: np.ndarray
    factor_identity: str
    pw_basis_identity: str
    state_authority: str
    source: dict
    frame: dict
    selected: np.ndarray
    quadrature_weight: float
    operator_work: dict
    initial_potential: np.ndarray
    initial_direction: np.ndarray
    initial_shift: np.ndarray
    hamiltonian_specification: dict
    block_size: int
    _coupling: object = field(repr=False)
    _hamiltonian: object = field(repr=False)

    def _result(self, action, adjoint, specification, work, started, cpu_started):
        identity = _identity(specification)
        result = ActionFactors(action, adjoint, identity,
            identity + ':actual-Hilbert-adjoint', specification['units'],
            self.pw_basis_identity, self.state_authority)
        result.validate(self.psi.shape, self.pw_basis_identity)
        work.update(_no_replay_work())
        work.update(operator_specification=specification, operator_identity=identity,
            source_authority=self.state_authority, full_pw_rows=self.psi.shape[0],
            raw_columns=self.psi.shape[1], elapsed_seconds=perf_counter()-started,
            process_cpu_seconds=process_time()-cpu_started,
            peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
            shared_restoration_work=self.operator_work,
            restoration_repeated=False, physical_factor_rebuilds=0)
        return result, work

    def hamiltonian_actions(self):
        """Return (H(t0) Psi, H(t0)^dagger Psi) factors and invocation work."""
        started, cpu_started = perf_counter(), process_time()
        action, adjoint, work = _hamiltonian_actions(
            self._hamiltonian, self.psi, self.block_size)
        work['direction_preparation_seconds'] = 0.
        return self._result(action, adjoint, self.hamiltonian_specification,
                            work, started, cpu_started)

    def ionic_force_actions(self, direction):
        """Return factors for F(d)=-sum_i d_i dH/dR_i, including complex d.

        Direction is full original Cartesian (nat,3) or atom-major (3*nat,).
        Adjoint calls use BlochDirectionalAction.adjoint_apply, not F(d) again.
        """
        started, cpu_started = perf_counter(), process_time()
        direction = np.asarray(direction, dtype=np.complex128)
        nat = self._coupling.ndof // 3
        if direction.shape == (3*nat,):
            direction = direction.reshape(nat, 3)
        if direction.shape != (nat, 3) or not np.isfinite(direction).all():
            raise ValueError('finite full original nat by3 or atom-major flat direction required')
        t = perf_counter()
        prepared = self._coupling.ionic.prepare(direction)
        preparation_seconds = perf_counter()-t
        specification = dict(kind='minus-current-physical-ionic-Hamiltonian-derivative',
            derivative_direction=_direction_spec(prepared.direction), force_sign=-1,
            derivative_scale=1, external_q=[0, 0, 0], time_inverse_ry=0,
            units='Ry/bohr', basis_identity=self.pw_basis_identity,
            source=self.source['source'],
            physical_sources=list(self._coupling.factors.sources),
            matched_model_identity=self.source['model_identity'],
            original_state_authority=self.state_authority,
            atom_names=self.source['atom_names'], quadrature_applied=False,
            adjoint_api='BlochDirectionalAction.adjoint_apply')
        action, adjoint, work = _force_actions(prepared, self.psi, self.block_size)
        work.update(h0_actions=0, h0_forward_actions=0, h0_adjoint_actions=0,
                    direction_preparation_seconds=preparation_seconds)
        return self._result(action, adjoint, specification, work, started, cpu_started)


def _initial_stage(source, saved, original, record):
    """Join the saved t0 scalar stage to this exact original, without reduction."""
    mean = bound_json(_SAVED_MEAN)
    if (mean['status'] != 'complete' or mean['mode'] != 'evolve'
            or mean['case'] != source['case']
            or not _same_binding(mean['initialization'], source['initial_result'])
            or not _same_binding(mean['initial_boundary'], saved['joint_boundary'])):
        raise ValueError('saved t0 scalar stage belongs to another initialization')
    joint = DiskJointMean(authenticate(saved['joint_boundary']),
                          expected_identity=source['model_identity'], verify=False)
    if (joint.time != 0 or not _same_binding(joint.origin, source['original_boundary'])
            or len(joint.slots) != len(original.slots)):
        raise ValueError('saved t0 joint/original boundary lineage differs')
    for a, b in zip(joint.slots, original.slots):
        if (a.metadata() != b.metadata() or not _same_binding(a.frame, b.frame)
                or not _same_binding(a.coefficient, b.coefficient)):
            raise ValueError('t0 joint state is not the unchanged full original state')
    stage = mean['initial_reduction']
    if (stage['kind'] != 'full-initial-reduction' or stage['time_ry_inverse'] != 0
            or Path(stage['boundary']).resolve() != Path(saved['joint_boundary']['path']).resolve()):
        raise ValueError('literal saved t0 reduction required')
    arrays = {name: _load_array(stage['arrays'][name]) for name in
              ('potential', 'u', 'p', 'coordinate', 'momentum', 'shift')}
    for name in ('u', 'p', 'coordinate', 'momentum'):
        if not np.array_equal(arrays[name], getattr(joint, name)):
            raise ValueError('saved stage and joint coordinate differ: '+name)
    if (not np.array_equal(joint.u, original.u) or not np.array_equal(joint.p, original.p)
            or arrays['potential'].ndim != 2 or arrays['potential'].shape[1] != 2
            or record.component >= arrays['potential'].shape[1]
            or arrays['shift'].shape != (3,)):
        raise ValueError('original ionic/spin/field stage dimensions differ')
    provenance = dict(receipt=dict(_SAVED_MEAN), joint_boundary=saved['joint_boundary'],
        kind=stage['kind'], time_inverse_ry=0,
        arrays={name: stage['arrays'][name] for name in arrays},
        potential_convention='literal saved K(conditional density minus saved density)',
        direction_convention='literal saved original u_bohr',
        field_convention='literal saved Cartesian shift, not field coordinate')
    return arrays['potential'][:, record.component], arrays['u'], arrays['shift'], provenance


def restore_packet_hamiltonian(packet_path):
    """Restore the contract's original packet context, performing zero actions.

    S,D,R are consumed from its authenticated coordinates NPZ verbatim. Native
    transport is used only to recover/check its physical basis and full raw frame.
    H0 and Gamma coupling share one physical factor bank by object identity.
    """
    started, cpu_started = perf_counter(), process_time()
    packet_binding, packet = _read_json(packet_path)
    if (packet_binding['sha256'] != _ORIGINAL_PACKET_SHA256
            or packet['schema'] != 'P254-full-noise-probe-input-v1'
            or packet['status'] != 'complete' or packet['time_inverse_ry'] != 0
            or packet['external_q'] != [0, 0, 0] or packet['quadrature_applied'] is not False):
        raise ValueError('contract-authenticated complete original t0 Gamma packet required')
    source = packet['source']
    initial = bound_json(source['initial_result'])
    if (initial['schema'] != 'P254-matched-joint-case-v1'
            or initial['status'] != 'complete' or initial['mode'] != 'initialize'
            or initial['material_trajectory'] is not False or initial['case'] != source['case']):
        raise ValueError('original completed matched initialization required')
    saved = initial['saved']
    for name in ('model_identity', 'original_model_identity'):
        if saved[name] != source[name]:
            raise ValueError('packet original model authority differs: '+name)
    for name in ('original_boundary', 'preparation'):
        if not _same_binding(saved[name], source[name]):
            raise ValueError('packet original saved input differs: '+name)
    authenticate(saved['preparation'])
    original = DiskMean(authenticate(source['original_boundary']),
                        expected_identity=source['original_model_identity'], verify=False)
    if (original.time != 0 or len(original.slots) != 128
            or original.evidence['context_identity'] != source['model_identity']
            or original.evidence['case'] != source['case']):
        raise ValueError('literal original matched full128 t0 boundary required')
    index = source['original_record']
    if type(index) is not int or not 0 <= index < len(original.slots):
        raise ValueError('exact original record index required')
    record = original.slots[index]
    if (_identity(record.metadata()) != _identity(source['original_metadata'])
            or not _same_binding(record.frame, packet['frame'])
            or not _same_binding(record.coefficient, source['original_raw_d'])
            or record.gauge_identity != packet['factor_identity']
            or record.basis_identity != packet['pw_basis_identity']
            or record.weight != packet['quadrature_weight'] or record.shape[1] != 346):
        raise ValueError('packet frame/D/basis/gauge/weight differs from original record')
    psi = _load_array(packet['frame'])
    raw_d = _load_array(record.coefficient)
    coordinates = bound_arrays(packet['coordinates'])
    if set(coordinates) != {'gram', 'density_coefficients', 'active_coordinates', 'selected'}:
        raise ValueError('exact original saved S,D,R,selected archive required')
    gram, d, active, selected = (coordinates[name] for name in
        ('gram', 'density_coefficients', 'active_coordinates', 'selected'))
    if (psi.shape != record.shape or psi.dtype != np.complex128
            or gram.shape != raw_d.shape or d.shape != raw_d.shape
            or raw_d.shape != (record.shape[1], record.shape[1])
            or not np.array_equal(d, raw_d) or selected.dtype.kind not in 'iu'
            or selected.ndim != 1 or not np.array_equal(selected, source['selected'])
            or active.shape != (record.shape[1], selected.size)
            or any(not np.isfinite(v).all() for v in (gram, d, active))):
        raise ValueError('saved full original S,D,R coordinates differ')
    physical_spec = saved['specification']
    if not _same_binding(physical_spec['aggregate'], source['matched_scalar']):
        raise ValueError('packet matched scalar differs from original specification')
    _, registration = material_io.registration()
    control, background, ledger, rows, context = material_io.load_case(source['case'], registration)
    if (context != physical_spec['context'] or control.provenance() != source['control_provenance']
            or not _same_binding(context['background_receipt'], source['saved_background_receipt'])
            or not _same_binding(context['population_ledger'], source['population_ledger'])
            or not _same_binding(bound_json(context['background_receipt'])['archive'],
                                 source['saved_background_archive'])):
        raise ValueError('restored physical background/control lineage differs from packet')
    rows = tuple(sorted(rows, key=lambda row: (row[1].ispin, row[1].mesh_index)))
    aggregate = bound_json(source['matched_scalar'])
    validate_aggregate(aggregate, context, rows)
    if not _same_binding(aggregate['archive'], source['matched_scalar_archive']):
        raise ValueError('matched scalar archive differs from packet')
    slots = {slot_identity(slot): slot for _, slot, _ in rows}
    if len(slots) != 128 or set(slots) != {r.identity for r in original.slots}:
        raise ValueError('current and original full sector identities differ')
    slot = slots[record.identity]
    if (slot.ispin != source['spin'] or record.component != slot.ispin-1
            or slot.weight != record.weight or slot.mesh_index != tuple(source['mesh_index'])
            or slot.source_ik != source['source_ik']
            or not _same_binding(slot.source_binding, source['source'])
            or record.source_identity != slot.source_binding['sha256']
            or slot.orbital_conjugation != source['orbital_conjugation']
            or slot.operation.index != source['operation']
            or slot.reciprocal_shift.tolist() != source['reciprocal_shift']
            or slot.k_fractional.tolist() != source['k_fractional']
            or slot.gauge_policy != source['gauge_policy']):
        raise ValueError('original physical spin/slot/source/transport lineage differs')
    state_binding = source['source_state']
    if not _same_binding(physical_spec['source_states'][str((slot.ispin, slot.source_ik))], state_binding):
        raise ValueError('original matched source-state binding differs')
    state = bound_json(state_binding)
    sc = state['context']
    if (state['status'] != 'complete' or any(sc.get(k) != v for k, v in context.items())
            or not _same_binding(sc['original_source'], slot.source_binding)
            or (sc['spin'], sc['source_ik']) != (slot.ispin, slot.source_ik)
            or not _same_binding(state['archive'], source['source_archive'])):
        raise ValueError('source-state original physical authority differs')
    source_arrays = bound_arrays(state['archive'])
    expected_s, expected_d = (source_arrays[name] for name in ('gram', 'density_coefficients'))
    if slot.orbital_conjugation:
        expected_s, expected_d = expected_s.conj(), expected_d.conj()
    if not np.array_equal(gram, expected_s) or not np.array_equal(d, expected_d):
        raise ValueError('packet raw Gram/D differs from admitted transported source')
    authority = _identity(dict(initial=source['initial_result'], preparation=saved['preparation'],
        boundary=source['original_boundary'], record=record.metadata(), frame=record.frame,
        raw_d=record.coefficient, state=state_binding, archive=state['archive']))
    if authority != packet['state_authority']:
        raise ValueError('packet original state authority hash differs')
    del source_arrays, expected_s, expected_d, ledger
    wave = load_wavefunctions(slot.source_binding)
    native = {name: getattr(wave, name) for name in
              ('ik', 'ispin', 'ngw', 'igwx', 'npol', 'nbnd', 'gamma_only', 'scalef', 'abi')}
    if (native != source['native_wavefunction']
            or any(type(native[name]) is not int for name in ('ik', 'ispin', 'ngw', 'igwx', 'npol', 'nbnd'))):
        raise ValueError('native wavefunction metadata changed or lost exact integer types')
    transported = transport_real_frame(control.grid, slot, wave)
    basis = transported.basis
    if (basis.identity != record.basis_identity or transported.gauge_identity != record.gauge_identity
            or transported.coefficients.shape != record.shape
            or not np.array_equal(transported.coefficients, psi)
            or dict(basis.source_admission) != source['source_admission']):
        raise ValueError('native original full frame/basis/gauge/admission differs')
    del transported, wave
    potential, direction, shift, stage_source = _initial_stage(source, saved, original, record)
    if potential.shape != (background.fft.nnr,) or direction.shape != (3*control.reference.nat,):
        raise ValueError('saved t0 stage does not fit original physical grid/atoms')
    restoration_seconds = perf_counter()-started
    t = perf_counter()
    factors = build_physical_factors(background.inputs, basis)
    if _identity(list(factors.sources)) != _identity(source['operator']['physical_sources']):
        raise ValueError('current physical operator sources differ from authenticated packet')
    h0 = MaterialHamiltonian(background, factors, slot.ispin, 8)
    coupling = physical_coupling(background.inputs, basis, factors, block_size=8)
    frozen = PhysicalFrozenAction(h0, coupling, potential=potential,
        direction=direction, shift=shift, fft_block_size=8, projector_block_size=64)
    factor_seconds = perf_counter()-t
    specification = dict(kind='literal-t0-matched-physical-Hamiltonian', units='Ry',
        external_q=[0, 0, 0], time_inverse_ry=0, force_sign=1,
        basis_identity=basis.identity, source=dict(slot.source_binding),
        physical_sources=list(factors.sources), original_state_authority=authority,
        matched_model_identity=source['model_identity'], saved_stage=stage_source,
        scalar='T+saved-vrs+physical-UPF-nonlocal+saved-Kdelta-n+g(saved-u)+L(saved-shift)',
        quadrature_applied=False, eigenvalue_labels_used=False,
        adjoint_api='PhysicalFrozenAction.adjoint_apply')
    restored_source = dict(source, packet=packet_binding, coordinates=packet['coordinates'],
        initial_stage=stage_source, native_wavefunction=native,
        atom_names=list(control.reference.atom_names),
        hamiltonian_implementation=file_binding(Path(__file__)),
        hamiltonian_specification=specification)
    work = _no_replay_work()
    work.update(restoration_seconds=restoration_seconds, physical_factor_seconds=factor_seconds,
        elapsed_seconds=perf_counter()-started, process_cpu_seconds=process_time()-cpu_started,
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
        h0_actions=0, force_forward_actions=0, force_adjoint_actions=0,
        full_pw_rows=psi.shape[0], raw_columns=psi.shape[1],
        background_restorations=1, physical_factor_builds=1,
        hamiltonian_factory_builds=1, coupling_factory_builds=1,
        saved_t0_potential_nonzero=bool(np.any(potential != 0)),
        saved_t0_ionic_direction_nonzero=bool(np.any(direction != 0)),
        saved_t0_field_shift_nonzero=bool(np.any(shift != 0)))
    for array in (gram, d, active, selected):
        array.setflags(write=False)
    return PacketHamiltonianContext(psi, gram, d, active, packet['factor_identity'],
        basis.identity, authority, restored_source, packet['frame'], selected,
        record.weight, work, potential, direction, shift, specification, 8, coupling, frozen)
