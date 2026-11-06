"""Actual saved raw-state boundaries and exhaustive source-bound read-back.

The source Gram/D are restored, never reconstructed/refilled. Only the admitted
initial orbital transport is used. No H0 action or coefficient recomputation.
"""
from itertools import groupby
from pathlib import Path
import numpy as np
import matched_reference_io as io
from matched_joint_case import slot_identity
from substrate_framework.material_control_reference import bound_json, bound_arrays
from substrate_framework.material_reference_state import restore_source
from substrate_framework.qe_saved_native import load_wavefunctions
from substrate_framework.qe_bloch_frame import admit_source
from substrate_framework.qe_real_bloch_grid import transport_real_state
from substrate_framework.streamed_material_state import DiskMean, file_binding, authenticate
from substrate_framework.streamed_joint_state import DiskJointMean


def original_sectors(context, budget):
    key = lambda row: (row[1].ispin, row[1].source_ik)
    for source_key, group in groupby(context.rows, key=key):
        rows = tuple(group)
        state_record = bound_json(context.source_states[source_key])
        arrays = bound_arrays(state_record['archive'])
        if any(not np.array_equal(arrays['populations'], row[2]) for row in rows):
            raise ValueError('original source populations differ from exact slot ledger')
        slot = rows[0][1]
        source = load_wavefunctions(slot.source_binding)
        state = restore_source(arrays, source_identity=slot.source_binding['sha256'])
        admission = admit_source(context.control.reference, slot, source)
        for _, slot, _ in rows:
            budget.check()
            frame, target = transport_real_state(context.grid, slot, source, state, source_admission=admission)
            recorded = context.source_metadata[slot_identity(slot)]
            if (frame.basis.identity != recorded['basis_identity']
                    or frame.gauge_identity != recorded['frame_identity']
                    or frame.coefficients.shape != (recorded['rows'], recorded['columns'])):
                raise ValueError('original complete raw frame identity differs from aggregate lineage')
            metadata = dict(identity=slot_identity(slot), component=slot.ispin-1,
                weight=slot.weight, basis_identity=frame.basis.identity,
                gauge_identity=frame.gauge_identity, source_identity=state.source_identity)
            yield metadata, frame.coefficients, target.density_coefficients
            del frame, target
        del source, state, arrays


def preparation_arrays(context):
    arrays = context.ionic_preparation.archive_arrays()
    arrays.update({'field_'+name: np.asarray(value) for name, value in context.field_preparation.items()
                   if isinstance(value, np.ndarray)})
    return arrays


def initialize(context, directory, budget):
    directory = Path(directory)
    def sectors():
        for metadata, psi, d in original_sectors(context, budget):
            # Conservative: charge D for every slot even when DiskMean byte-deduplicates it.
            budget.reserve(psi.nbytes+d.nbytes+8192)
            yield metadata, psi, d
    ndof = len(context.reference.masses)
    evidence = dict(case=context.case, specification=context.specification,
        original_state_receipts={str(k): v for k, v in context.source_states.items()},
        no_gram_inverse=True, no_refilling=True, full_raw_columns=True,
        initial_orbital_transport_only=True)
    budget.reserve(len(__import__('json').dumps(evidence).encode())+65536)
    original = DiskMean.create(directory/'original', sectors(), identity=context.identity,
        time=0., u=np.zeros(ndof), p=np.zeros(ndof), evidence=evidence)
    state = DiskJointMean.wrap_original(file_binding(original.path), directory/'initial',
        identity=context.identity, coordinate=context.field_preparation['coordinate_mean'],
        momentum=context.field_preparation['momentum_mean'], expected_slots=context.expected_slots,
        evidence=evidence, budget=budget)
    arrays = preparation_arrays(context)
    size = sum(a.nbytes for a in arrays.values())+65536
    budget.reserve(size)
    prep = io.save_arrays(directory/'preparation.npz', arrays, size)
    return dict(case=context.case, model_identity=context.identity,
        original_boundary=file_binding(original.path), joint_boundary=file_binding(state.path),
        preparation=prep, specification=context.specification)


def read_back(context, record, budget):
    if (record['case'] != context.case or record['model_identity'] != context.identity
            or record['specification'] != context.specification):
        raise ValueError('saved boundary is not this exact original case/model')
    original = DiskMean(authenticate(record['original_boundary']), expected_identity=context.identity)
    state = DiskJointMean(authenticate(record['joint_boundary']), expected_identity=context.identity)
    if (state.time != 0 or original.time != 0 or state.history
            or state.origin != record['original_boundary']
            or tuple(s.identity for s in state.slots) != context.expected_slots
            or tuple(s.identity for s in original.slots) != context.expected_slots):
        raise ValueError('complete original t0 lineage differs')
    for name, expected in (('u', np.zeros_like(context.reference.masses)),
                           ('p', np.zeros_like(context.reference.masses)),
                           ('coordinate', context.field_preparation['coordinate_mean']),
                           ('momentum', context.field_preparation['momentum_mean'])):
        if not np.array_equal(getattr(state, name), expected):
            raise ValueError('original canonical initial means differ: '+name)
    for name in ('u', 'p'):
        if not np.array_equal(getattr(state, name), getattr(original, name)):
            raise ValueError('original ionic means differ')
    checked = 0
    for (metadata, psi, d), saved, raw in zip(original_sectors(context, budget), state.slots, original.slots):
        expected = dict(metadata, shape=psi.shape)
        if saved.metadata() != expected or raw.metadata() != expected or saved.frame != raw.frame or saved.coefficient != raw.coefficient:
            raise ValueError('raw slot metadata/reference differs')
        actual_psi, actual_d = state.arrays(saved)
        if not np.array_equal(actual_psi, psi) or not np.array_equal(actual_d, d):
            raise ValueError('stored frame or original complex D changed')
        checked += 1
    if checked != 128:
        raise ValueError('incomplete original frame read-back')
    saved_preparation = bound_arrays(record['preparation'])
    expected_preparation = preparation_arrays(context)
    if (set(saved_preparation) != set(expected_preparation)
            or any(not np.array_equal(saved_preparation[k], v) for k, v in expected_preparation.items())):
        raise ValueError('actual original ionic/field preparation read-back differs')
    return dict(checked_original_slots=checked, raw_frames_and_complex_D_exact=True,
                canonical_means_exact=True, preparation_arrays_exact=True,
                material_trajectory=False, coefficient_recomputation=False)
