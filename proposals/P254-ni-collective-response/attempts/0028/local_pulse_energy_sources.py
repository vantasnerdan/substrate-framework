"""Authenticate the saved matched preparation and actual Gamma seed metadata."""
from pathlib import Path

def restore_registration(path):
    """Authenticate every original archive before consuming any parameter."""
    from joint_mean_energy_sources import bound_json, same_binding
    from substrate_framework.streamed_material_state import authenticate, file_binding
    binding = file_binding(path)
    reg = bound_json(binding)
    if reg['schema'] != 'P254-local-pulse-energy-preregistration-v1':
        raise ValueError('original local-pulse registration required')
    expected_cases = {'loaded-ni31h6', 'pristine-ni32', 'empty-vacancy-ni31', 'dispersed-ni32h6'}
    if len(reg['cases']) != 4 or {c['case'] for c in reg['cases']} != expected_cases:
        raise ValueError('all four distinct original registered cases required')
    for source in reg['sources']:
        authenticate(source)
    case_registration = bound_json(reg['case_registration'])
    model_proof = bound_json(reg['matched_model_source'])
    adoption = bound_json(reg['actual_gamma_source_adoption'])
    if reg['matched_preparation_declaration'] != 'u=p=q0=pi0=0; explicit zero support':
        raise ValueError('registered centered preparation differs')
    if (model_proof['status'] != 'complete' or model_proof['mode'] != 'read-back'
            or model_proof['case'] != 'loaded-ni31h6'
            or model_proof['proof']['canonical_means_exact'] is not True
            or model_proof['proof']['preparation_arrays_exact'] is not True):
        raise ValueError('source-bound loaded centered-preparation read-back required')
    aggregates = {}
    for case in reg['cases']:
        aggregate = bound_json(case['aggregate'])
        authenticate(case['archive'])
        if (aggregate['status'] != 'complete' or aggregate['complete_case'] is not True
                or aggregate['context']['case'] != case['case']
                or aggregate['slot_count'] != 128
                or not same_binding(aggregate['archive'], case['archive'])
                or not same_binding(case_registration['cases'][case['case']], case['aggregate'])):
            raise ValueError('original complete case/archive lineage differs: '+case['case'])
        aggregates[case['case']] = aggregate
    return reg, aggregates, case_registration, model_proof, adoption, binding


def centered_preparation(reg, case_registration, model_proof):
    """Corroborate the declaration with genuine loaded initial-boundary means.

    Other cases use the registered common centered construction; this is not
    represented as a new read-back of their many-electron preparations.
    """
    import numpy as np
    from joint_mean_energy_sources import bound_json, same_binding
    init = bound_json(model_proof['initialization'])
    if init['status'] != 'complete' or init['mode'] != 'initialize' or init['case'] != 'loaded-ni31h6':
        raise ValueError('original loaded initialization required')
    saved = init['saved']
    spec = saved['specification']
    loaded = next(c for c in reg['cases'] if c['case'] == 'loaded-ni31h6')
    if (spec['initial_means'] != reg['matched_preparation_declaration']
            or not same_binding(spec['registration'], reg['case_registration'])
            or not same_binding(spec['aggregate'], loaded['aggregate'])
            or spec['preparation'] != case_registration['construction']['preparation']):
        raise ValueError('loaded preparation does not bind common matched declaration')
    boundary = bound_json(saved['joint_boundary'])
    original = bound_json(saved['original_boundary'])
    if (boundary['schema'] != 'streamed-joint-mean-v1' or boundary['complete'] is not True
            or boundary['model_identity'] != saved['model_identity']
            or boundary['time_ry_inverse'] != 0 or original['time_ry_inverse'] != 0
            or not same_binding(boundary['original_electronic_boundary'], saved['original_boundary'])):
        raise ValueError('genuine loaded zero-time boundary identity differs')
    actual = {key: boundary[key] for key in ('u_bohr', 'p_canonical', 'field_coordinate', 'field_momentum')}
    for name, value in actual.items():
        array = np.asarray(value)
        if array.ndim != 1 or not array.size or not np.isfinite(array).all() or np.any(array != 0):
            raise ValueError('saved loaded initial mean is not centered: '+name)
    if original['u_bohr'] != actual['u_bohr'] or original['p_canonical'] != actual['p_canonical']:
        raise ValueError('original electronic/ionic boundary canonical means differ')
    return dict(declaration=reg['matched_preparation_declaration'],
        case_registration=reg['case_registration'], construction=case_registration['construction'],
        actual_loaded_readback=reg['matched_model_source'], initialization=model_proof['initialization'],
        loaded_joint_boundary=saved['joint_boundary'], loaded_original_boundary=saved['original_boundary'],
        actual_loaded_mean_fields=actual,
        evidence_scope='Common source-authenticated centered matched construction for all four cases; actual saved loaded zero-time means additionally checked. Not equilibrium; no new read-back of other cases or original orbital payloads.')
def saved_gamma_comparison(reg, adoption, initialization_binding, loaded, fourier_function):
    import numpy as np
    from joint_mean_energy_sources import bound_json, same_binding
    from substrate_framework.partitioned_joint_slots import digest
    from substrate_framework.streamed_material_state import file_binding
    request = reg['actual_gamma_seed_request']
    identity = reg['current_source_identity']
    if (adoption['status'] != 'metadata-adoption-complete' or adoption['current_identity'] != identity
            or adoption['physical_input_unchanged'] is not True):
        raise ValueError('adoption002 source identity/physical input changed')
    # The actual admitted checkpoint is adoption002.new_seed, never an old seed
    # found by a recursive search, and never an invented response/zero array.
    seed = bound_json(adoption['new_seed'])
    root = Path(adoption['new_seed']['path']).parents[1]
    request_binding = file_binding(root/'request.json')
    actual_request = bound_json(request_binding)
    scientific = bound_json(adoption['new_scientific_request'])
    impulse = reg['impulse']
    if (actual_request != request or root.name != digest(request)
            or request['identity'] != identity or request['time'] != 0
            or request['kind'] != 'coherent-initial' or request['dynamics'] != 'screened-joint'
            or request['key']['q'] != [0., 0., 0.]
            or request['key']['amplitude'] != impulse['amplitude'] or request['key']['cell'] != impulse['cell']
            or seed['schema'] != 'streamed-joint-response-v1' or seed['complete'] is not True
            or seed['identity'] != identity or seed['time'] != 0 or seed['kind'] != request['kind']
            or seed['q'] != request['key']['q'] or seed['provenance'] != request['provenance']
            or scientific['identity'] != identity or scientific['case'] != loaded['case']
            or not same_binding(scientific['initialization'], initialization_binding)
            or not same_binding(scientific['baseline'], request['provenance']['baseline'])
            or any(scientific[key] != impulse[key] for key in ('atom', 'axis', 'amplitude', 'cell'))
            or scientific['q'] != request['key']['q']):
        raise ValueError('actual adopted Gamma seed/request/initialization binding differs')
    ndof = 3*len(loaded['original_parameter_fields']['atomic_masses_amu']['values'])
    direction = np.asarray(request['key']['direction'])
    coordinate_index = 3*loaded['struck_atom']['zero_based_array_index']+('x', 'y', 'z').index(impulse['axis'])
    if (direction.shape != (ndof,) or direction[coordinate_index] != 1
            or np.count_nonzero(direction) != 1):
        raise ValueError('actual saved Gamma source is not the registered original atom/axis')
    gamma_coefficient = next(row['coefficient'] for row in
        fourier_function(impulse['mesh'], impulse['cell'], impulse['amplitude'])
        if row['mesh_index'] == (0, 0, 0))
    small, checks = {}, {}
    for name, size in (('u', ndof), ('p', ndof), ('coordinate', 3), ('momentum', 3)):
        encoded = np.asarray(seed['small'][name])
        if encoded.shape != (size, 2) or not np.isfinite(encoded).all():
            raise ValueError('actual saved seed small field dimension/nonfinite value: '+name)
        value = encoded[:, 0]+1j*encoded[:, 1]
        small[name] = value
        checks[name] = bool(np.array_equal(value, gamma_coefficient*direction)) if name == 'p' else bool(np.all(value == 0))
    if not all(checks.values()):
        raise ValueError('actual saved Gamma canonical seed metadata does not match the derived pulse')
    gamma_dose = float(np.sum(np.abs(small['p'])**2))/(2*loaded['struck_atom']['mass_Ry'])
    return dict(status='matched', adoption_binding=reg['actual_gamma_source_adoption'],
        actual_checkpoint_binding=adoption['new_seed'], actual_seed_request_binding=request_binding,
        actual_scientific_request_binding=adoption['new_scientific_request'],
        original_adoption=adoption, actual_seed_request=actual_request, actual_scientific_request=scientific,
        actual_checkpoint=seed, checks=checks,
        expected_Gamma_momentum_coefficient=[gamma_coefficient.real, gamma_coefficient.imag],
        quadratic_share_from_actual_saved_p_metadata_ry=gamma_dose,
        ledger_Gamma_share_ry=loaded['Fourier']['one_Gamma_share']['initial_quadratic_input_ry'],
        discrepancy_ry=gamma_dose-loaded['Fourier']['one_Gamma_share']['initial_quadratic_input_ry'],
        scope='Actual adopted q0 checkpoint JSON metadata only; leg payload bindings retained, not authenticated/read here. No fabricated frames; no q0 propagation result or whole-local response claim.')
