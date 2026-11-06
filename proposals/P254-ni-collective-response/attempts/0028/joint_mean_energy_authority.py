"""Restore saved mean authority using an explicitly selected current proof.

Historical physical execution and current source admission remain distinct.
Original registration, mean, stage, scalar and implementation bindings are never
rewritten. This module does not load a frame or apply a physical operator.
"""
from hashlib import sha256
import json

from joint_mean_energy_sources import bound_json, same_binding
from substrate_framework.streamed_material_state import authenticate


def restore_current_sources(registration_binding, qualification_binding):
    """Return the authentic saved model under the admitted current execution."""
    from matched_joint_case import load_case
    from run_joint_q_material import bind_saved_model
    from run_mixed_q_material import admit_oracle, existing_dependencies

    reg = bound_json(registration_binding)
    if reg.get('schema') != 'P254-saved-coupled-mean-energy-registration-v1':
        raise ValueError('registered saved coupled mean energy required')
    original = reg['source_model']
    mean = bound_json(reg['source_mean'])
    initialization = bound_json(original['initialization'])
    saved = initialization['saved']
    if (mean.get('status') != 'complete' or mean.get('mode') != 'evolve'
            or initialization.get('status') != 'complete'
            or initialization.get('mode') != 'initialize'
            or mean['case'] != initialization['case'] or saved['case'] != mean['case']
            or not same_binding(mean['initialization'], original['initialization'])
            or saved['model_identity'] != original['identity']
            or saved['specification'] != original['specification']):
        raise ValueError('original complete mean/initialization/model differs')

    # Use the existing current-source qualification, not a re-stamped historical
    # proof. Its own validator authenticates its actual current source closure.
    qualification = bound_json(qualification_binding)
    admit_oracle(qualification, existing_dependencies())
    retained_binding = qualification['retained_physical_operator_proof']
    retained = bound_json(retained_binding)
    if (retained.get('status') != 'passed'
            or retained.get('schema') != 'P254-joint-partitioned-q-oracle-v1'
            or retained['model_identity'] != original['identity']
            or not same_binding(retained['initialization'], original['initialization'])):
        raise ValueError('retained original physical-scalar proof differs')
    continuation = qualification['source_compatibility']
    if continuation != retained['source_compatibility']:
        raise ValueError('current and retained physical source compatibility differ')
    before, after = continuation['before'], continuation['after']
    for binding in original['implementation'] + reg['implementation_sources']:
        if binding['path'] == before['path']:
            if not same_binding(binding, before):
                raise ValueError('registered history-owner scope differs from admitted original')
            authenticate(after)
        else:
            authenticate(binding)
    source_digest = sha256(json.dumps({b['path']: b['sha256'] for b in
        original['implementation']}, sort_keys=True).encode()).hexdigest()
    source_digest = sha256((source_digest + json.dumps(
        original['column_limits'], sort_keys=True)).encode()).hexdigest()
    if source_digest != original['source_identity']:
        raise ValueError('original implementation/column-limit source identity differs')
    for binding in reg['archive_sources']:
        authenticate(binding)

    context = bind_saved_model(load_case(mean['case'],
        reality_tolerance=original['reality_tolerance']), saved, continuation)
    model = context.make_mean(None, block_size=original['block_size'])
    if model.identity != original['identity']:
        raise ValueError('restored model authority differs')
    if len(reg['stages']) != 2 or len(reg['slot_sources']) != 2 or len(mean['branches']) != 1:
        raise ValueError('exact registered initial/final pair required')
    if reg['stages'][0] != mean['initial_reduction']:
        raise ValueError('registered initial stage differs from original mean')
    if (not same_binding(mean['initial_boundary'], reg['slot_sources'][0]['boundary'])
            or not same_binding(saved['joint_boundary'], reg['slot_sources'][0]['boundary'])
            or not same_binding(mean['branches'][0]['boundary'], reg['slot_sources'][1]['boundary'])
            or not same_binding(retained['history'], reg['slot_sources'][1]['boundary'])
            or reg['stages'][1]['energy_total_ry'] != mean['branches'][0]['diagnostics']['energy_total_ry']):
        raise ValueError('original saved mean boundary/energy provenance differs')
    provenance = dict(registration=dict(registration_binding), source_mean=reg['source_mean'],
        initialization=original['initialization'], original_model=original,
        current_source_compatibility_proof=dict(qualification_binding),
        retained_physical_operator_proof=retained_binding,
        source_compatibility=continuation,
        authenticated_parameter_archives=reg['archive_sources'],
        authenticated_scalar_implementation=reg['implementation_sources'])
    return reg, model, saved, provenance
