"""Original t=0 full-PW force factors, without H0 or preparation replay.

``prepare_full_noise_probe`` returns the numeric interface consumed by the
retained/omitted algebra. Its arrays are unweighted: the original sector measure
is metadata, never a force multiplier. Only one original record is acted on.
"""
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import resource
import sys
from time import perf_counter, process_time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT / 'src'), str(HERE.parent / '0028')]

import numpy as np
import matched_reference_io as material_io
from matched_joint_case import slot_identity
from matched_joint_sources import source_catalog, validate_aggregate, validate_atoms
from substrate_framework.material_control_reference import bound_arrays, bound_json
from substrate_framework.material_physical_operators import build_physical_factors, physical_coupling
from substrate_framework.physical_factor_preparation import ActionFactors
from substrate_framework.qe_real_bloch_grid import transport_real_frame
from substrate_framework.qe_saved_native import load_wavefunctions
from substrate_framework.streamed_material_state import DiskMean, authenticate, file_binding


@dataclass(frozen=True)
class FullNoiseProbeInputs:
    psi: np.ndarray
    gram: np.ndarray
    density_coefficients: np.ndarray
    active_coordinates: np.ndarray
    selected: np.ndarray
    probes: dict
    factor_identity: str
    pw_basis_identity: str
    state_authority: str
    frame: dict
    source: dict
    quadrature_weight: float
    operator_work: dict

    def algebra_arguments(self):
        """Keyword arguments only; do not construct another full action bank."""
        return dict(psi=self.psi, gram=self.gram,
                    density_coefficients=self.density_coefficients,
                    active_coordinates=self.active_coordinates, probes=self.probes,
                    factor_identity=self.factor_identity,
                    pw_basis_identity=self.pw_basis_identity,
                    state_authority=self.state_authority)


def _identity(value):
    return sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def _same_binding(left, right):
    return (Path(left['path']).resolve() == Path(right['path']).resolve()
            and int(left['bytes']) == int(right['bytes']) and left['sha256'] == right['sha256'])


def _read_json(path):
    # The completed union manifest exceeds bound_json's 32 MiB receipt limit.
    binding = file_binding(path)
    with authenticate(binding).open() as stream:
        value = json.load(stream)  # native integer metadata stays integer
    return binding, value


def _saved_selection(manifest, original, record, bands):
    specification = manifest['specification']
    if (manifest['schema'] != 'p254-physical-activeband-cache-v1'
            or manifest['status'] != 'complete' or specification['bands'] != bands
            or not _same_binding(specification['original_boundary'], original)):
        raise ValueError('complete matching saved-band manifest and original boundary required')
    slots = [(key, item) for key, item in manifest['items'].items() if key.startswith('slot/')]
    identities = [item['metadata']['identity'] for _, item in slots]
    if len(identities) != 128 or len(set(identities)) != 128:
        raise ValueError('saved coefficients must identify all128 original sectors exactly once')
    matches = [(key, item) for key, item in slots if item['metadata']['identity'] == record.identity]
    if len(matches) != 1:
        raise ValueError('original identity has no unique saved active selection')
    key, item = matches[0]
    metadata = item['metadata']
    if (metadata['basis_identity'] != record.basis_identity
            or metadata['frame_identity'] != record.gauge_identity
            or not _same_binding(metadata['frame'], record.frame)
            or not _same_binding(metadata['raw_d'], record.coefficient)):
        raise ValueError('saved active selection belongs to another raw frame/D/basis/gauge')
    raw = metadata['selected']
    if (len(raw) != bands or any(type(i) is not int for i in raw)
            or len(set(raw)) != bands or min(raw) < 0 or max(raw) >= record.shape[1]):
        raise ValueError('saved selected indices must be distinct exact in-range integers')
    selected = np.asarray(raw, dtype=np.int64)
    arrays = bound_arrays(item['archive'])
    if not np.array_equal(arrays['selected'], selected) or arrays['selected'].dtype.kind not in 'iu':
        raise ValueError('saved selected metadata disagrees with authenticated coefficient archive')
    return key, item, selected, set(identities)


def _coordinates(gram, selected):
    """The A29 symmetric-Gram convention, with the SAVED selection, once."""
    if gram.ndim != 2 or gram.shape[0] != gram.shape[1] or not np.isfinite(gram).all():
        raise ValueError('finite original square admitted Gram required')
    tolerance = 128 * np.finfo(float).eps * len(gram) * (1 + np.max(np.abs(gram)))
    if np.max(np.abs(gram - gram.conj().T)) > tolerance:
        raise ValueError('original Gram is not Hermitian; no repair permitted')
    values, vectors = np.linalg.eigh(gram)
    if values.min() <= 0:
        raise ValueError('original Gram is not positive definite; no truncation permitted')
    # Same inverse square root as A29; no Psi^dagger Psi and no Q allocation.
    inverse = (vectors / np.sqrt(values)) @ vectors.conj().T
    return inverse[:, selected]


def _force_actions(prepared, psi, block_size=8):
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
    # The operator API is the POSITIVE derivative. Force is its negative,
    # independently on both legs; no spin, cell, or quadrature scaling.
    action *= -1
    adjoint *= -1
    action.setflags(write=False)
    adjoint.setflags(write=False)
    return action, adjoint, dict(force_forward_actions=1, force_adjoint_actions=1,
        forward_block_calls=blocks, adjoint_block_calls=blocks,
        forward_action_columns=psi.shape[1], adjoint_action_columns=psi.shape[1],
        forward_seconds=forward_seconds, adjoint_seconds=adjoint_seconds,
        block_size=block_size)


def _implementation_bindings():
    paths = {str(HERE / 'full_noise_probe_inputs.py'), str(HERE / 'export_full_noise_probe.py')}
    for name, module in tuple(sys.modules.items()):
        path = getattr(module, '__file__', None)
        if path and path.endswith('.py') and (name.startswith('substrate_framework.')
                or name.startswith('matched_') or name in ('material_reference_output', 'material_pilot_io')):
            paths.add(str(Path(path).resolve()))
    return [file_binding(path) for path in sorted(paths)]


def prepare_full_noise_probe(initial_result, coefficient_manifest, original_record, atom, axis, bands=12):
    """Restore one original sector and evaluate its actual force/adjoint pair.

    No Native backend, H0 apply, density/kernel reduction, population selection,
    ionic/field preparation constructor, or whole projection exporter is called.
    The caller owns resource limits and any subsequent covariance contraction.
    """
    started, cpu_started = perf_counter(), process_time()
    if type(original_record) is not int or original_record < 0:
        raise ValueError('nonnegative exact original record index required')
    if type(atom) is not int or axis not in ('x', 'y', 'z') or bands != 12:
        raise ValueError('one original atom, Cartesian axis and saved12-band selection required')
    initial_binding, initial = _read_json(initial_result)
    if (initial['schema'] != 'P254-matched-joint-case-v1' or initial['status'] != 'complete'
            or initial['mode'] != 'initialize' or initial['material_trajectory'] is not False):
        raise ValueError('completed original matched initialization required')
    saved = initial['saved']
    original_binding = saved['original_boundary']
    original = DiskMean(authenticate(original_binding),
        expected_identity=saved['original_model_identity'], verify=False)
    if (original.time != 0 or np.any(original.u) or np.any(original.p)
            or original.evidence['context_identity'] != saved['model_identity']
            or original.evidence['case'] != initial['case'] or len(original.slots) != 128
            or original_record >= len(original.slots)):
        raise ValueError('literal original matched complete128 t=0 boundary required')
    record = original.slots[original_record]
    authenticate(record.frame)
    authenticate(record.coefficient)
    authenticate(saved['preparation'])
    psi, raw_d = original.arrays(record)
    coefficient_binding, manifest = _read_json(coefficient_manifest)
    key, coefficient_item, selected, coefficient_identities = _saved_selection(
        manifest, original_binding, record, bands)
    coefficient_spec = manifest['specification']
    physical_spec = saved['specification']
    # Execution implementation changed after initialization/A29. It is bound
    # separately, not allowed to replace the unchanged physical model identity.
    for name in ('case', 'context', 'aggregate', 'source_states', 'scalar', 'preparation'):
        if coefficient_spec['context'][name] != physical_spec[name]:
            raise ValueError('A29 and original matched physical specification differ: ' + name)
    coefficient_case_identity = coefficient_spec['case_identity']
    selection_recipe = coefficient_spec['selection_recipe']
    del manifest
    restore_start = perf_counter()
    _, registration = material_io.registration()
    control, background, ledger, rows, context = material_io.load_case(initial['case'], registration)
    if context != physical_spec['context']:
        raise ValueError('restored saved background/control differs from original matched scalar')
    # Authenticate the matched scalar, but do not run its ionic/field Gaussian
    # constructors. Original product preparation is already bound above.
    aggregate = bound_json(physical_spec['aggregate'])
    rows = tuple(sorted(rows, key=lambda row: (row[1].ispin, row[1].mesh_index)))
    validate_aggregate(aggregate, context, rows)
    aggregate_arrays = bound_arrays(aggregate['archive'])
    validate_atoms(control, background, aggregate_arrays)
    del aggregate_arrays
    states, _ = source_catalog(aggregate, context, rows, ledger)
    slots = {slot_identity(slot): slot for _, slot, _ in rows}
    if set(slots) != coefficient_identities or set(slots) != {r.identity for r in original.slots}:
        raise ValueError('original/coefficient/current physical slot identity sets differ')
    slot = slots[record.identity]
    if (not 1 <= atom <= control.reference.nat or record.component != slot.ispin - 1
            or record.weight != slot.weight or record.source_identity != slot.source_binding['sha256']):
        raise ValueError('original atom/spin/source/quadrature differs from current physical slot')
    metadata = coefficient_item['metadata']
    state_binding = states[(slot.ispin, slot.source_ik)]
    state = bound_json(state_binding)
    if (not _same_binding(state_binding, physical_spec['source_states'][str((slot.ispin, slot.source_ik))])
            or not _same_binding(state_binding, metadata['source_state'])
            or not _same_binding(state['archive'], metadata['source_archive'])
            or not _same_binding(slot.source_binding, metadata['source'])
            or metadata['orbital_conjugation'] != slot.orbital_conjugation
            or metadata['operation'] != slot.operation.index
            or metadata['reciprocal_shift'] != slot.reciprocal_shift.tolist()):
        raise ValueError('saved coefficient source/conjugation/operation lineage differs')
    source_arrays = bound_arrays(state['archive'])
    gram = source_arrays['gram'].conj() if slot.orbital_conjugation else source_arrays['gram']
    expected_d = (source_arrays['density_coefficients'].conj() if slot.orbital_conjugation
                  else source_arrays['density_coefficients'])
    if gram.shape != raw_d.shape or not np.array_equal(raw_d, expected_d):
        raise ValueError('original raw D differs from admitted source-state transport')
    active_coordinates = _coordinates(gram, selected)
    wave = load_wavefunctions(slot.source_binding)
    native = {name: getattr(wave, name) for name in
              ('ik', 'ispin', 'ngw', 'igwx', 'npol', 'nbnd', 'gamma_only', 'scalef', 'abi')}
    frame = transport_real_frame(control.grid, slot, wave)
    if (frame.basis.identity != record.basis_identity or frame.gauge_identity != record.gauge_identity
            or frame.coefficients.shape != record.shape or not np.array_equal(frame.coefficients, psi)
            or dict(frame.basis.source_admission) != metadata['source_admission']):
        raise ValueError('actual native original full frame/basis/gauge/admission differs')
    basis = frame.basis
    del frame, wave, source_arrays, expected_d
    restore_seconds = perf_counter() - restore_start
    t = perf_counter()
    factors = build_physical_factors(background.inputs, basis)
    coupling = physical_coupling(background.inputs, basis, factors)
    direction = np.zeros((control.reference.nat, 3))
    direction[atom - 1, 'xyz'.index(axis)] = 1.
    prepared = coupling.ionic.prepare(direction)
    operator_spec = dict(kind='minus-current-physical-ionic-Hamiltonian-derivative',
        atom=atom, species=control.reference.atom_names[atom - 1], axis=axis,
        external_q=[0, 0, 0], time_inverse_ry=0, units='Ry/bohr',
        basis_identity=basis.identity, source=dict(slot.source_binding),
        physical_sources=list(factors.sources), derivative_direction=direction.tolist(),
        force_sign=-1, derivative_scale=1, quadrature_applied=False,
        adjoint_api='BlochDirectionalAction.adjoint_apply')
    operator_identity = _identity(operator_spec)
    factor_seconds = perf_counter() - t
    action, adjoint, work = _force_actions(prepared, psi)
    authority = _identity(dict(initial=initial_binding, preparation=saved['preparation'],
        boundary=original_binding, record=record.metadata(), frame=record.frame,
        raw_d=record.coefficient, state=state_binding, archive=state['archive']))
    probe_name = '%s%d/%sGamma' % (control.reference.atom_names[atom - 1], atom, axis)
    probe = ActionFactors(action, adjoint, operator_identity,
        operator_identity + ':actual-Hilbert-adjoint', 'Ry/bohr', basis.identity, authority)
    probe.validate(psi.shape, basis.identity)
    source = dict(case=initial['case'], model_identity=saved['model_identity'],
        original_model_identity=original.identity, initial_result=initial_binding,
        original_boundary=original_binding, preparation=saved['preparation'],
        matched_scalar=physical_spec['aggregate'], matched_scalar_archive=aggregate['archive'],
        saved_background_receipt=context['background_receipt'],
        saved_background_archive=bound_json(context['background_receipt'])['archive'],
        control_provenance=control.provenance(), population_ledger=context['population_ledger'],
        original_record=original_record, original_metadata=record.metadata(),
        original_raw_d=record.coefficient, spin=slot.ispin, mesh_index=list(slot.mesh_index),
        k_fractional=slot.k_fractional.tolist(), source_ik=slot.source_ik,
        source=dict(slot.source_binding), source_state=state_binding, source_archive=state['archive'],
        orbital_conjugation=slot.orbital_conjugation, operation=slot.operation.index,
        reciprocal_shift=slot.reciprocal_shift.tolist(), gauge_policy=slot.gauge_policy,
        source_admission=dict(basis.source_admission), native_wavefunction=native,
        coefficient_manifest=coefficient_binding, coefficient_case_identity=coefficient_case_identity,
        coefficient_item=key, coefficient_archive=coefficient_item['archive'],
        selection_recipe=selection_recipe, selected=selected.tolist(), operator=operator_spec,
        implementation=_implementation_bindings(),
        state_convention='rho=Psi D Psi^dagger; exactly empty orthogonal finite-PW complement',
        measure_convention='literal original sector weight; no implicit spin or cell multiplier')
    work.update(restoration_seconds=restore_seconds, physical_factor_seconds=factor_seconds,
        elapsed_seconds=perf_counter() - started, process_cpu_seconds=process_time() - cpu_started,
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        h0_actions=0, native_calls=0, trajectory_actions=0, initial_kernel_reductions=0,
        gram_recomputations=0, population_refills=0, preparation_replays=0,
        original_inputs_unchanged=True, full_pw_rows=psi.shape[0], raw_columns=psi.shape[1])
    for array in (gram, active_coordinates, selected):
        array.setflags(write=False)
    return FullNoiseProbeInputs(psi, gram, raw_d, active_coordinates, selected,
        {probe_name: probe}, record.gauge_identity, basis.identity, authority,
        record.frame, source, record.weight, work)
