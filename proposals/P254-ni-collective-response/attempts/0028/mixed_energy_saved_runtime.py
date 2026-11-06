"""Assemble the existing physical energy provider around a saved raw boundary.

No raw seed is initialized, no response RHS is evaluated and no propagation is
entered.  Work lives in the independent consumer namespace; original producer
clock, source receipts and partially completed energy work remain unchanged.
"""
from copy import copy
from pathlib import Path

import numpy as np

from mixed_material_dynamics import NativeMixedDynamics
from mixed_response_storage import COMPONENTS, MappedMixedState
from run_mixed_q_material import (admit_oracle, bind_saved_model, checked_json,
    existing_dependencies, load_case, validate_physical_preparation)
from substrate_framework.partitioned_joint_slots import digest
from substrate_framework.partitioned_q_work import QWork
from substrate_framework.streamed_material_state import authenticate, file_binding


def build_saved_runtime(admission, *, budget, cpus, work_directory):
    """Return genuine dynamics/work/evidence for the admitted physical state."""
    producer = Path(admission['producer_directory']).resolve()
    work_directory = Path(work_directory).resolve()
    if work_directory == producer or work_directory.is_relative_to(producer):
        raise ValueError('independent energy work must not mutate the producer namespace')
    request = admission['owner']['request']
    if admission['source_archive'] is None:
        for source in request['source_bindings']:
            authenticate(source)
    else:
        # Boundary admission already checked every exact original member. The
        # same archive bytes retain that evidence; current code is bound apart.
        authenticate(admission['source_archive'])
    initialization_binding, initialized = checked_json(
        request['initialization']['path'], request['initialization']['sha256'])
    if (initialized.get('status') != 'complete' or initialized.get('mode') != 'initialize'
            or initialized.get('case') != request['case']
            or initialized['saved']['case'] != request['case']):
        raise ValueError('original selected matched initialization required')
    saved = initialized['saved']
    history_binding, _ = checked_json(request['history']['path'], request['history']['sha256'])
    proof_binding, proof = checked_json(request['execution_proof']['path'],
                                      request['execution_proof']['sha256'])
    admit_oracle(proof, existing_dependencies())
    block_size = request['block_size']
    case = bind_saved_model(load_case(request['case']), saved, proof.get('source_compatibility'))
    original = case.make_history(history_binding, budget, block_size=block_size)
    histories, geometries = {}, {}
    q = np.asarray(request['q'])
    for component, vector in (('plus', q), ('minus', -q), ('mixed', np.zeros(3))):
        key = tuple(vector)
        if key not in geometries:
            history = original if not geometries else copy(original)
            case.bind_response(history, vector)
            geometries[key] = history
        histories[component] = geometries[key]
    validate_physical_preparation(case, histories, saved)
    # A readout must use actual saved mean frames, never silently generate a
    # trajectory/interior approximation just because its time is in range.
    time = admission['state']['time']
    for history in histories.values():
        history.clear_interval()
        if history._saved_frames(time) is None:
            raise ValueError('saved physical boundary lacks genuine saved mean frames')
    dynamics = NativeMixedDynamics(case, histories, budget=budget, block_size=block_size)
    with MappedMixedState(admission['state_binding'], budget=budget) as state:
        dynamics._validate_state(state)
        if state.identity != case.identity or state.time != time:
            raise ValueError('saved raw boundary differs from the original matched provider')
    identity = digest(dict(schema='P254-saved-mixed-energy-work-v1',
        original_producer=admission['owner_binding'], state=admission['state_binding'],
        native_physics=dynamics.source_identity, consumer_runtime=file_binding(__file__)))
    work = QWork(work_directory, budget, cpus, identity)
    evidence = dict(schema='P254-saved-mixed-energy-runtime-v1',
        original_owner=admission['owner_binding'], state=admission['state_binding'],
        boundary_kind=admission['boundary_kind'], matched_case=request['case'],
        original_producer_source_archive=admission['source_archive'],
        source_policy='Historical source bytes are original-owner-bound; imported consumer physics and current source qualification are authenticated separately.',
        model_identity=case.identity, initialization=initialization_binding,
        completed_mean=history_binding, execution_proof=proof_binding,
        native_physics_source_identity=dynamics.source_identity,
        work_identity=work.identity, work_directory=str(work.directory),
        original_baseline=dict(initial=file_binding(original.initial.path),
            final=file_binding(original.final.path), stages=original.records),
        current_runtime_source=file_binding(__file__), propagation_performed=False,
        raw_initialization_performed=False, original_owner_modified=False,
        producer_energy_cache_reused=False,
        cache_policy='Separate source-bound consumer work only. Existing original partial energy/cache work is retained, not adopted or relabelled as current; any repeated readout work stays charged.',
        prepared_provenance='Actual endpoint-specific prepared binding and scalar receipts are published by mixed_energy_readout, not inferred here.')
    return dynamics, work, evidence
