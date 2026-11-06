"""Exact saved-mean inputs for the Native mixed scalar readout.

Only the registered completed initial/final reductions and actual stored
midpoint JointStages are supported. An internal boundary without its own saved
reduction, or any other interior time, is unavailable: a frozen midpoint is not
an endpoint, and canonical means are never interpolated. This adapter restores
source authority using the existing saved-mean accounting APIs; it performs no
frame, H0, density, gradient, field-vertex or kernel contraction.
"""
from types import SimpleNamespace

import numpy as np

from substrate_framework.streamed_material_state import file_binding
from joint_mean_energy_sources import (
    bound_json, parameter_identities, read_stage, same_binding)
from joint_mean_energy_authority import restore_current_sources


class MeanSourceUnavailable(ValueError):
    """A requested physical mean has no authenticated saved read-time source."""

    def __init__(self, message, *, requested_time, available):
        super().__init__(message)
        self.requested_time = requested_time
        self.available = available


def ledger_binding(dynamics):
    """Require the explicitly selected authenticated mean ledger."""
    supplied = getattr(dynamics, 'energy_mean_readout', None)
    if supplied is None:
        raise ValueError('an explicit authenticated current saved-mean energy ledger is required')
    # Do not silently rebind changed bytes at an explicitly supplied source.
    bound_json(supplied)
    return dict(supplied)


def _validate_stage(stage, dynamics):
    case = dynamics.case
    shapes = dict(density=case.reference.nbar.shape, potential=case.reference.nbar.shape,
        gradient=case.reference.masses.shape, field_gradient_electronic=(3,),
        u=case.reference.masses.shape, p=case.reference.masses.shape,
        coordinate=(case.field.size,), momentum=(case.field.size,), shift=(3,))
    for name, shape in shapes.items():
        value = np.asarray(getattr(stage, name))
        if value.shape != shape or not np.isfinite(value).all():
            raise ValueError('actual saved mean shape/nonfinite value: '+name)
    if not np.array_equal(stage.shift, stage.coordinate @ case.field.polarizations):
        raise ValueError('saved mean shift differs from its actual canonical mode coordinate')


def saved_mean(state, dynamics, ledger_source):
    """Return the genuine read-time JointStage, grid measure and source evidence."""
    ledger = bound_json(ledger_source)
    if (ledger.get('schema') != 'P254-saved-coupled-mean-energy-readout-v1'
            or ledger.get('status') != 'complete'):
        raise ValueError('complete original saved-mean energy ledger required')
    registration = ledger['provenance']['registration']
    reg, model, saved, provenance = restore_current_sources(
        registration, ledger['provenance']['current_source_compatibility_proof'])
    if (not same_binding(registration, provenance['registration'])
            or ledger['provenance'] != provenance):
        raise ValueError('ledger and restored original model/source digest authority differ')
    case, baseline = dynamics.case, dynamics.baseline
    if case.identity != model.identity:
        raise MeanSourceUnavailable(
            'Selected saved mean ledger has no source for requested Native model '+case.identity+
            '; registered source model is '+model.identity+'. No surrogate baseline is substituted.',
            requested_time=state.time,
            available=dict(model_identity=model.identity, ledger=ledger_source))
    if (case.identity != reg['source_model']['identity']
            or case.specification != reg['source_model']['specification']
            or dynamics.block_size != reg['source_model']['block_size']):
        raise ValueError('mixed dynamics changed the registered original saved model authority')
    ref = case.control.reference
    measure = float(abs(np.linalg.det(ref.cell_bohr))/int(np.prod(ref.fft_smooth)))
    authority = SimpleNamespace(reference=case.reference, coefficients=case.coefficients,
        charges=case.charges, field=case.field, measure=measure,
        reality_tolerance=case.reality_tolerance)
    if parameter_identities(authority) != parameter_identities(model):
        raise ValueError('mixed dynamics changed original saved scalar parameters/grid measure')
    if (not same_binding(dynamics.final_binding, reg['slot_sources'][1]['boundary'])
            or not same_binding(file_binding(baseline.initial.path), reg['slot_sources'][0]['boundary'])
            or not same_binding(baseline.initial.origin, saved['original_boundary'])
            or tuple(s.identity for s in baseline.initial.slots) != model.expected_slots):
        raise ValueError('mixed history differs from the original registered completed mean')
    original_boundary = bound_json(reg['slot_sources'][0]['boundary'])
    if len(original_boundary['slots']) != dynamics.count:
        raise ValueError('original complete source coverage differs')
    for original, slot in zip(original_boundary['slots'], baseline.initial.slots):
        if (original['coefficient'] != slot.coefficient or original['frame'] != slot.frame
                or original['identity'] != slot.identity):
            raise ValueError('mixed history replaced an original raw D/frame/source')
    available = dict(boundaries=[dict(time=s['time'], boundary=s['boundary'])
                                for s in reg['slot_sources']],
        actual_midpoints=[dict(time=r['time_ry_inverse'], interval=i)
                          for i, r in enumerate(baseline.records)])
    saved_boundary = baseline._saved_frames(state.time)
    stage_record = None
    for index, source in enumerate(reg['slot_sources']):
        if state.time != source['time']:
            continue
        if (saved_boundary is None
                or not same_binding(file_binding(saved_boundary.path), source['boundary'])):
            raise ValueError('requested endpoint is not the registered actual saved boundary')
        stage, unused_saved_bare_scalar, evidence = read_stage(reg, model, saved, index)
        label = ('initial', 'final')[index]
        if (ledger[label]['stage'] != reg['stages'][index]
                or ledger[label]['model_identity'] != model.identity
                or ledger[label]['source_identity'] != reg['source_model']['source_identity']
                or not same_binding(ledger[label]['boundary'], source['boundary'])):
            raise ValueError('saved ledger endpoint stage/model/digest differs')
        stage_record = reg['stages'][index]
        semantics = 'actual completed '+label+' boundary reduction; not frozen-midpoint substitution'
        break
    else:
        index, prior, stage = baseline._stage(state.time)
        record = baseline.records[index]
        if saved_boundary is not None or state.time != record['time_ry_inverse']:
            raise MeanSourceUnavailable(
                'No saved complete mean density/potential/u/p/a/Pi/g/ell reduction at requested '
                'time; only registered initial/final boundaries and exact actual stored midpoint '
                'stages are available. No endpoint substitution or canonical interpolation.',
                requested_time=state.time, available=available)
        if (record.get('kind') != 'actual-frozen-midpoint'
                or not same_binding(record['prior_boundary_binding'], file_binding(prior.path))
                or record['start_time_ry_inverse'] != prior.time
                or state.time != prior.time+record['dt_ry_inverse']/2):
            raise ValueError('actual matching saved midpoint source/time differs')
        stage_record = record
        evidence = dict(time_ry_inverse=state.time, interval=index,
            prior_boundary=record['prior_boundary_binding'],
            completed_mean=dynamics.final_binding, stage=record,
            source_identity=reg['source_model']['source_identity'], model_identity=model.identity,
            actual_reads='actual matching stored JointStage; no mean contraction or interpolation')
        semantics = 'actual stored midpoint JointStage at its exact saved time only'
    _validate_stage(stage, dynamics)
    # The restored model is authority only. Its lazy operator factory is never
    # opened; numerical means above come exclusively from the actual saved stage.
    return stage, measure, dict(ledger=ledger_source, provenance=provenance,
        read_time_semantics=semantics, stage_record=stage_record,
        actual_stage=evidence, available_read_times=available,
        unused_saved_mean_bare_scalar='not used or recomputed; endpoint accounting API authenticates existing scalar receipts')
