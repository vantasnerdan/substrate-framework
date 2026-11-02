"""Native coherent response of an original ionic impulse from a saved dual.

No electronic legs are replayed: their original perturbation is exactly zero.
The Hilbert pairing is dual.p.conjugate() @ seed.p, not a covariance contraction.
The affected spatial oracle compares this to independent discrete forward RK4
and dense physical-density evolution, including nonzero-cell complex phases.
"""
import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE)]

import numpy as np

from substrate_framework.material_control_reference import bound_json
from substrate_framework.streamed_material_state import file_binding
from joint_spatial_provenance import same_binding
import joint_spatial_dual_authority

SCHEMA = 'P254-source-local-dual-coherent-readout-v2'
SEED_SCHEMA = 'P254-original-local-ionic-momentum-seed-v1'
CONVENTION = 'vdot(original-t0-dual.p, amplitude*exp(-2pi i q.cell)*unit_atom_axis/sqrt(N)); all other initial perturbations zero'


def canonical_impulse(context, q, ncells, seed):
    """Pure scalar contraction; exact complex convention exercised by oracle."""
    fields = {'schema', 'initial_perturbation', 'case', 'model_identity', 'initialization',
              'atom', 'axis', 'amplitude', 'cell'}
    if (set(seed) != fields or seed.get('schema') != SEED_SCHEMA
            or seed.get('initial_perturbation') != 'ionic-momentum-only'
            or type(ncells) is not int or ncells < 1 or type(seed['atom']) is not int
            or seed['atom'] < 1 or seed['axis'] not in ('x', 'y', 'z')
            or type(seed['amplitude']) not in (int, float) or not np.isfinite(seed['amplitude'])):
        raise ValueError('explicit original ionic-momentum-only seed required')
    cell, q = np.asarray(seed['cell']), np.asarray(q, float)
    if (cell.shape != (3,) or not np.isfinite(cell).all() or not np.array_equal(cell, np.rint(cell))
            or q.shape != (3,) or not np.isfinite(q).all()):
        raise ValueError('original integer kick cell and finite q required')
    encoded = np.asarray(context['small']['p'], float)
    if encoded.ndim != 2 or encoded.shape[1] != 2 or len(encoded) % 3 or not np.isfinite(encoded).all():
        raise ValueError('native finite complex canonical dual momentum required')
    momentum = encoded[:, 0]+1j*encoded[:, 1]
    coordinate = 3*(seed['atom']-1)+'xyz'.index(seed['axis'])
    if coordinate >= len(momentum):
        raise ValueError('original seed atom outside canonical dual geometry')
    impulse = seed['amplitude']*np.exp(-2j*np.pi*q@cell)/np.sqrt(ncells)
    return momentum[coordinate].conjugate()*impulse


def build_readout(owner_binding, context_binding, seed_binding, oracle_binding, *, snapshots=()):
    """Read an independently complete pullback; covariance is not an input."""
    context, request, p, saved, mesh, slots = joint_spatial_dual_authority.complete_pullback(
        owner_binding, context_binding, snapshots)
    seed, proof = bound_json(seed_binding), bound_json(oracle_binding)
    sources = [file_binding(Path(__file__).resolve()),
               file_binding(Path(joint_spatial_dual_authority.__file__).resolve())]
    if (proof.get('schema') != 'P254-joint-spatial-readout-oracle-v1'
            or proof.get('status') != 'passed'
            or not all(any(same_binding(source, b) for b in proof['source_bindings'])
                       for source in sources)
            or not all(name in proof['checks'] for name in (
                'source_local_dual_discrete_forward_q0', 'source_local_dual_discrete_forward_q1',
                'source_local_dual_discrete_forward_q2'))):
        raise ValueError('executed source-bound dual/coherent spatial oracle required')
    if (not same_binding(seed['initialization'], p['initialization'])
            or seed['case'] != p['case'] or seed['model_identity'] != p['model_identity']):
        raise ValueError('original seed initialization/case/scalar differs')
    selectors = p['request']['selectors']
    value = canonical_impulse(context, request['q'], int(np.prod(mesh)), seed)
    return dict(schema=SCHEMA, status='complete', complete=True, owner=owner_binding,
        seed=seed_binding, numerical_oracle=oracle_binding, consumer_sources=sources,
        source_snapshots=list(snapshots), context=context_binding, context_request=context['request'],
        convention=CONVENTION, execution_scope='complete one-sided pullback; covariance not consumed',
        q=request['q'], case=p['case'], model_identity=p['model_identity'], scalar=p['scalar'],
        physical_specification=saved['specification'], initialization=p['initialization'],
        baseline=p['baseline'], matched_aggregate=p['matched_aggregate'],
        dynamics=request['dynamics'], observable=p['observable'],
        atom=selectors['atom'], axis=selectors['axis'], time_ry_inverse=context['terminal_time'],
        time_grid=request['times'], all_explicit_edges=p['all_explicit_edges'],
        fixed_original_D=p['fixed_original_D'], original_cell_count=int(np.prod(mesh)),
        original_source_slots=slots, source_bindings=p['source_bindings'],
        readout=dict(real=float(value.real), imag=float(value.imag)),
        numerical_accuracy_certified=False, whole_spatial_reconstruction=False)


def validate_readout(record, seed):
    if record.get('schema') != SCHEMA or record.get('status') != 'complete' or record.get('complete') is not True:
        raise ValueError('native complete dual coherent readout required')
    original_seed = bound_json(record['seed'])
    if original_seed != seed:
        raise ValueError('dual coherent original seed differs')
    rebuilt = build_readout(record['owner'], record['context'], record['seed'],
                            record['numerical_oracle'], snapshots=record['source_snapshots'])
    if rebuilt != record:
        raise ValueError('dual coherent publication differs from original context/seed/proof')
    return dict(readout=record['readout'], amplitude=seed['amplitude'], cell=seed['cell'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('owner', 'context', 'seed', 'oracle'):
        parser.add_argument('--'+name, type=Path, required=True)
        parser.add_argument('--'+name+'-sha256', required=True)
    parser.add_argument('--source-snapshots', type=Path)
    parser.add_argument('--source-snapshots-sha256')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    bindings = []
    for name in ('owner', 'context', 'seed', 'oracle'):
        binding = file_binding(getattr(args, name))
        if binding['sha256'] != getattr(args, name+'_sha256'):
            raise ValueError('explicit '+name+' digest differs')
        bindings.append(binding)
    snapshots = []
    if (args.source_snapshots is None) != (args.source_snapshots_sha256 is None):
        raise ValueError('source snapshot file and explicit digest must be supplied together')
    if args.source_snapshots is not None:
        binding = file_binding(args.source_snapshots)
        if binding['sha256'] != args.source_snapshots_sha256:
            raise ValueError('explicit source snapshot manifest digest differs')
        snapshots = bound_json(binding)
        if not isinstance(snapshots, list):
            raise ValueError('source snapshot manifest must contain the original/snapshot binding list')
    result = build_readout(*bindings, snapshots=snapshots)
    if result['original_cell_count'] != 64 or result['original_source_slots'] != 128:
        raise ValueError('material CLI requires original64-q/128-spin source geometry')
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')


if __name__ == '__main__':
    main()
