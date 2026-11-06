"""Read actual saved force gradients/RHS; never replay retired source moments."""
import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE)]

import numpy as np
import clamped_scientific_readout as saved
import saved_force_origins as algebra
from saved_canonical_response_readout import geometry, same_binding
from mixed_response_storage import MappedMixedState
from substrate_framework.material_control_reference import bound_json
from substrate_framework.streamed_material_state import authenticate, file_binding


def vector(binding):
    value = np.load(authenticate(binding), allow_pickle=False)
    if value.shape != (111,) or value.dtype != np.complex128 or not np.isfinite(value).all():
        raise ValueError('actual finite111 complex128 force array required')
    return value


def read_forces(registration_binding, oracle_binding):
    registration, oracle = bound_json(registration_binding), bound_json(oracle_binding)
    if (registration['schema'] != 'P254-saved-direct-force-origins-registration-v1'
            or oracle['schema'] != 'P254-saved-direct-force-origins-oracle-v1'
            or oracle['status'] != 'passed'):
        raise ValueError('registered saved sources and executed force-origin oracle required')
    for source in oracle['source_bindings']+registration['source_equations']:
        authenticate(source)
    if not any(source == file_binding(algebra.__file__) for source in oracle['source_bindings']):
        raise ValueError('force oracle does not bind this actual force kernel')
    inputs = registration['source_inputs']
    energy, prepared = bound_json(inputs['energy']), bound_json(inputs['prepared'])
    metadata = bound_json(inputs['metadata'])
    g = geometry(metadata)
    if (energy['schema'] != 'P254-native-mixed-material-energy-readout-v1'
            or energy['status'] != 'complete' or energy['complete'] is not True
            or not same_binding(energy['input'], inputs['mixed_endpoint'])
            or not same_binding(energy['prepared'], inputs['prepared'])
            or prepared['schema'] != 'P254-native-mixed-material-prepared-v1'
            or prepared['complete'] is not True
            or prepared['identity'] != energy['raw_provenance']['native_context']
            or prepared['identity']['source']['model'] != g['model']
            or not same_binding(prepared['identity']['input'], inputs['mixed_endpoint'])
            or metadata['identity']['extra_delta_curvature'] != 'exact zero'):
        raise ValueError('genuine same-model endpoint/prepared global moments required')
    retirement = bound_json(registration['retirement'])
    if not same_binding(retirement['barrier'], inputs['prepared']):
        raise ValueError('retired source payloads are not replaced by this prepared barrier')
    records = prepared['moment_records']
    if len(records) != prepared['moment_work']['count'] or len(records) != 128:
        raise ValueError('all original source moment receipts required')
    for index, (record, binding) in enumerate(zip(records, registration['original_moment_receipts'], strict=True)):
        if record['index'] != index or record['complete'] is not True or bound_json(binding) != record:
            raise ValueError('original complete source receipt/order differs')
    for source in prepared['identity']['source']['adapter_sources'].values():
        authenticate(source)
    with np.load(authenticate(metadata['identity']['archive']), allow_pickle=False) as archive:
        c = archive['common_q_grid_ry_bohr2'][0, 0, 0].reshape(111, 111)
        mixed = archive['mixed_ionic_field'].reshape(111, 3)
    sources = [file_binding(path) for path in (Path(__file__), Path(algebra.__file__))]
    result, checks = {}, {}
    with MappedMixedState(inputs['mixed_endpoint']) as state:
        if (state.kind != 'mixed-rk4-physical-endpoint' or state.identity != g['model']
                or state.time != energy['time_ry_inverse'] or np.any(state.q)):
            raise ValueError('actual same-time physical Gamma endpoint required, not RK node/accumulator')
        for direction in ('plus', 'minus', 'mixed'):
            small = state.small[direction]
            gradient = vector(prepared['arrays'][direction+'-gradient'])
            cached = np.sqrt(g['ncells'])*vector(prepared['arrays'][direction+'-p'])
            values = algebra.force_components(c, mixed, small.u, small.coordinate, gradient, ncells=g['ncells'])
            difference = values['total']-cached
            scale = np.sqrt(g['ncells'])*(np.abs(c) @ np.abs(small.u)+np.abs(gradient)+np.abs(mixed) @ np.abs(small.coordinate))
            tolerance = 1024*np.finfo(float).eps*scale
            if np.any(np.abs(difference) > tolerance):
                raise ValueError('direct origins disagree with genuine saved physical p-RHS')
            result[direction] = {name: value.reshape(37, 3) for name, value in values.items()}
            checks[direction] = dict(difference_from_cached_physical_rhs=difference,
                roundoff_comparison_threshold=tolerance, interpretation='arithmetic comparison only, not physical accuracy or convergence')
    control = bound_json(inputs['clamped_control'])
    if (control['model_identity'] != g['model'] or control['time_ry_inverse'] != state.time
            or control['original_seed'] != dict(atom=1,axis='x',amplitude=1.,cell=[0,0,0])
            or not same_binding(control['source_inputs']['mixed_initial'],
                bound_json(inputs['canonical_report'])['source_inputs']['mixed_initial'])):
        raise ValueError('actual same-time/seed frozen-electron counterfactual required')
    frozen = {}
    for direction in ('plus', 'minus', 'mixed'):
        small = control['endpoint'][direction]
        u = np.asarray(small['u']['real'])+1j*np.asarray(small['u']['imag'])
        a = np.asarray(small['coordinate']['real'])+1j*np.asarray(small['coordinate']['imag'])
        values = algebra.force_components(c,mixed,u,a,np.zeros(111,complex),ncells=g['ncells'])
        frozen[direction] = {name: value.reshape(37,3) for name,value in values.items()}
    original = bound_json(inputs['right_coherent'])
    seed = bound_json(original['seed'])
    if (original['model_identity'] != g['model'] or any(original['q'])
            or original['observable'] != 'original-H37x-total-force'
            or seed['atom'] != 1 or seed['axis'] != 'x' or seed['amplitude'] != 1 or seed['cell'] != [0,0,0]
            or bound_json(original['context'])['terminal_time'] != state.time):
        raise ValueError('independent original H37x coherent publication differs')
    old = complex(original['readout']['real'],original['readout']['imag'])
    h = result['plus']['total'][36,0]
    if sources != [file_binding(source['path']) for source in sources]:
        raise ValueError('force consumer changed during evaluation')
    return dict(schema='P254-saved-direct-force-origins-v1',status='complete',
        registration=registration_binding,numerical_oracle=oracle_binding,source_inputs=inputs,
        model_identity=g['model'],time_ry_inverse=state.time,q=[0,0,0],
        original_atom_indices=g['atom_ids'],atom_names=g['species'],units='Ry/bohr collective Fourier force coefficient',
        direct_force_origins=result,frozen_electron_counterfactual=frozen,
        coupled_minus_frozen={direction:result[direction]['total']-frozen[direction]['total'] for direction in result},
        independent_H37x=dict(actual_direct_force=h,original_dual_force=old,difference=h-old),
        saved_rhs_comparisons=checks,source_bindings=sources,
        source_retirement=registration['retirement'],
        interpretation='Direct retained GLOBAL electronic gradient aggregate, not an electronic residual or replayed retired source arrays. Independent plus/minus/mixed physical response coefficients and each origin retained; unexcited baseline not added. Frozen-electron gradients are zero by that explicitly different control model. Raw counterfactual differences and independent dual comparison do not certify material accuracy, stationary screening, finite-dose force/energy, quantum noise, heat, spatial coherence or lattice ranking.')


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('registration','oracle'):
        parser.add_argument('--'+name,type=Path,required=True)
        parser.add_argument('--'+name+'-sha256',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(); bindings={}
    for name in ('registration','oracle'):
        bindings[name]=file_binding(getattr(args,name))
        if bindings[name]['sha256'] != getattr(args,name+'_sha256'):
            raise ValueError('explicit '+name+' digest differs')
    report=read_forces(bindings['registration'],bindings['oracle'])
    with args.output.open('x') as stream:
        json.dump(saved.encode(report),stream,indent=2,allow_nan=False); stream.write('\n')
    print(json.dumps(file_binding(args.output)))
