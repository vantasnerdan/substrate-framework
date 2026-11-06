"""Original four-scalar initial force-response jets, never finite-time propagation."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE)]

import numpy as np
import clamped_scientific_readout as saved
import initial_force_response as jet
from substrate_framework.material_control_reference import bound_json
from substrate_framework.streamed_material_state import authenticate, file_binding


def one_case(bindings):
    original = bound_json(bindings['original_result'])
    if original['status'] != 'complete' or original['inputs'] != bindings['metadata']:
        raise ValueError('actual completed original scalar must bind these metadata bytes')
    inputs = bound_json(bindings['metadata'])
    identity = inputs['identity']
    model = sha256(json.dumps(identity['matched_joint_specification'], sort_keys=True).encode()).hexdigest()
    if (inputs['schema'] != 'P254-matched-clamped-inputs-v1'
            or model != identity['matched_joint_identity'] or model != bindings['model_identity']
            or identity['control'] != 'periodic' or identity['mass']['variant'] != 'original'
            or identity['extra_delta_curvature'] != 'exact zero'
            or identity['phase_space_order'] != saved.ION_ORDER or identity['field_order'] != saved.FIELD_ORDER):
        raise ValueError('original matched periodic scalar, canonical order and zero added DeltaC required')
    aggregate = bound_json(bindings['aggregate'])
    if aggregate['archive'] != identity['archive']:
        raise ValueError('original scalar archive differs')
    ids = saved.real_metadata(identity['atom_indices'], 'original atom IDs', integer=True)
    names = saved.decoded(identity['atom_names'])
    masses = saved.real_metadata(identity['preparation']['mass_ry'], 'Cartesian Ry masses')
    cells = saved.real_metadata(identity['cell_indices'], 'original cell indices', integer=True)
    mesh = cells.max(axis=0)+1
    if (not np.array_equal(ids, np.arange(1, len(ids)+1)) or names[0] != 'Ni'
            or masses.shape != (3*len(ids),) or not np.array_equal(cells, np.indices(tuple(mesh)).reshape(3, -1).T)):
        raise ValueError('original atom/xyz ordering, Ni1 seed and full cell grid required')
    with np.load(authenticate(bindings['archive']), allow_pickle=False) as archive:
        c = archive['common_q_grid_ry_bohr2'].reshape(*mesh, len(masses), len(masses))
    source_bindings = [file_binding(path) for path in (Path(__file__), Path(jet.__file__), Path(saved.__file__))]
    result = jet.full_grid(c, masses, mesh=mesh, cell=np.zeros(3, int), amplitude=1., coordinate=0)
    force = result['force_rate_cells'].reshape(len(cells), len(ids), 3)
    summaries = {}
    for species in ('all', 'Ni', 'H'):
        mask = np.ones(len(ids), bool) if species == 'all' else names == species
        if not np.any(mask):
            summaries[species] = dict(applicable=False, reason='species absent in original model')
            continue
        weights = np.sum(np.abs(force[:, mask, :])**2, axis=(1, 2))
        total = float(np.sum(weights))
        summaries[species] = dict(applicable=True, original_atom_indices=ids[mask],
            squared_force_rate_weight_by_cell=weights, total_squared_force_rate_weight=total,
            off_seed_cell_fraction=None if total == 0 else float(np.sum(weights[1:])/total),
            weight_units='Ry^4/bohr^2; magnitude diagnostic, not energy or correlation')
    if source_bindings != [file_binding(record['path']) for record in source_bindings]:
        raise ValueError('consumer source changed during initial-jet evaluation')
    return dict(case=identity['case'], model_identity=model, source_inputs=bindings,
        original_atom_indices=ids, atom_names=names,
        positions_bohr=identity['positions_bohr'], cell_bohr=identity['cell_bohr'], masses_ry=masses,
        seed=dict(atom=1, axis='x', amplitude=1., cell=[0, 0, 0], initial_electronic='zero',
                  initial_displacement='zero', initial_field_coordinate_and_momentum='zero'),
        initial_force_response=0., response_derivative=result, species=summaries,
        source_bindings=source_bindings,
        initial_mechanism='Original momentum-only seed gives zero initial electronic tangent and electronic RHS, zero displacement and field-coordinate variation, but u_prime=p/M. Existing canonical total-force probe has zero p/Pi coefficients. Therefore its first derivative is the physical matched C(q) contraction; later electronic/field feedback is NOT removed or bounded by this result.')


def read_initial_response(registration_binding, oracle_binding):
    registration, oracle = bound_json(registration_binding), bound_json(oracle_binding)
    if (registration['schema'] != 'P254-exact-initial-force-response-registration-v1'
            or oracle['schema'] != 'P254-exact-initial-force-response-oracle-v1'
            or oracle['status'] != 'passed'):
        raise ValueError('registered original sources and executed initial-jet oracle required')
    for source in oracle['source_bindings']:
        authenticate(source)
    if not any(source == file_binding(jet.__file__) for source in oracle['source_bindings']):
        raise ValueError('executed oracle does not bind the actual initial-jet kernel')
    for source in registration['source_equations']:
        authenticate(source)
    cases = {name: one_case(bindings) for name, bindings in registration['source_inputs'].items()}
    return dict(schema='P254-exact-initial-force-response-v1', status='complete',
        registration=registration_binding, numerical_oracle=oracle_binding, cases=cases,
        scalar_source_equations=registration['source_equations'],
        interpretation='Exact initial jet of the SAME coupled canonical-force observable for all original q and cells. Not a frozen-electron finite-time result. Original nonstationarity and raw complex defects retained. No remainder bound, positive-time spatial propagation, quantum covariance, coherence window, finite-dose redistribution, heat, material accuracy or lattice ranking.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('registration', 'oracle'):
        parser.add_argument('--'+name, type=Path, required=True)
        parser.add_argument('--'+name+'-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    bindings = {}
    for name in ('registration', 'oracle'):
        bindings[name] = file_binding(getattr(args, name))
        if bindings[name]['sha256'] != getattr(args, name+'_sha256'):
            raise ValueError('explicit '+name+' digest differs')
    report = read_initial_response(bindings['registration'], bindings['oracle'])
    with args.output.open('x') as stream:
        json.dump(saved.encode(report), stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(file_binding(args.output)))
