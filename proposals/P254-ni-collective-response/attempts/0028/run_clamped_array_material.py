#!/usr/bin/env python3
"""Incremental actual saved-input clamped-electron quantum CONTROL, not joint dynamics."""
from pathlib import Path
import argparse
import io
import json
import os
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
for directory in (ROOT/'src', HERE, HERE.parent/'0025', HERE.parent/'0024'):
    sys.path.insert(0,str(directory))

import numpy as np
from material_run_contract import binding, bound_json, recover_block
from material_reference_output import load_arrays
from run_full_field_moments_material import pack, unpack
from substrate_framework.material_clamped_array import ClampedArray


def source_bindings():
    paths = [HERE/'clamped-array-registration.json', Path(__file__),
             HERE/'run_clamped_array_oracle.py',
             ROOT/'src/substrate_framework/material_clamped_array.py',
             HERE/'run_full_field_moments_material.py',
             HERE.parent/'0025/material_run_contract.py',
             HERE.parent/'0025/material_reference_output.py']
    return [binding(p) for p in paths]


def authenticated_arrays(receipt):
    archive = Path(receipt['archive']['path'])
    actual = recover_block(archive.parent,archive.stem,receipt['recipe_sha256'],receipt['metadata'])
    if actual != receipt:
        raise ValueError('saved array receipt authentication differs')
    return load_arrays(receipt)


def load_material(reference_path, joint_path, ionic_path):
    records = {name:binding(path) for name,path in (
        ('reference',reference_path),('joint',joint_path),('ionic',ionic_path))}
    registration = bound_json(binding(HERE/'clamped-array-registration.json'))
    for name,record in records.items():
        if record != registration['settled_inputs'][name]:
            raise ValueError('input differs from settled registration: '+name)
    reference = bound_json(records['reference'])
    joint = unpack(bound_json(records['joint']))
    ionic = bound_json(records['ionic'])
    metadata = reference['metadata']
    if (not metadata['complete_electronic_coverage'] or len(metadata['slots']) != 128
            or len(metadata['source_records']) != 20 or joint['status'] != 'complete'
            or joint['periodic_cells'] != 64 or joint['full_reference'] != records['reference']
            or joint['ionic_preparation'] != records['ionic'] or ionic['status'] != 'complete'
            or not ionic['covariance_published']
            or ionic['receipt']['metadata']['reference'] != records['reference']
            or joint['field_preparation_status'] != 'conditional-preparation-defined'):
        raise ValueError('complete settled 20-source/128-slot scalar and product preparations required')
    arrays, prep = authenticated_arrays(reference), authenticated_arrays(ionic['receipt'])
    means = bound_json(ionic['receipt']['metadata']['initial_means'])
    if np.any(means['u_bohr']) or np.any(means['p_canonical']):
        raise ValueError('settled zero-mean sudden release required')
    if arrays['common_force_lattice_ry_bohr2'].shape != (4,4,4,37,3,37,3):
        raise ValueError('full 64-cell/111-coordinate original lattice required')
    if not np.array_equal(prep['masses_ry'],joint['original_atoms']['masses_ry']):
        raise ValueError('original mass identity differs')
    field = joint['maxwell']
    coefficients = joint['coefficients']
    initial = joint['initial_join']
    # Saved field momentum RHS = -N ell in these original Cartesian polarizations.
    # The original settled aggregate uses the complete identity polarization basis.
    if not np.array_equal(field['polarizations'],np.eye(3)):
        raise ValueError('settled material requires its original Cartesian field basis')
    model = ClampedArray(arrays['common_force_lattice_ry_bohr2'],prep['masses_ry'],
        joint['original_atoms']['valence_charges_e'],coefficients['delta_curvature'],
        coefficients['contact'],coefficients['mixed'],initial['force_per_cell_ry_bohr'],
        -np.asarray(initial['field_momentum_rhs'])/64,np.asarray(field['polarizations']),
        field['mass'],prep['local_ordered'],joint['field_preparation']['ordered_covariance'],
        prep['omega_ry'])
    if np.any(joint['field_preparation']['coordinate_mean']) or np.any(joint['field_preparation']['momentum_mean']):
        raise ValueError('settled field must have its recorded centered means')
    identity = dict(original_atoms=joint['original_atoms'],positions_bohr=arrays['positions_bohr'],
        cell_bohr=arrays['cell_bohr'],cell_indices=model.indices,
        phase_space_order='u[atom-original-order,xyz],p[atom-original-order,xyz]',
        field_order='a_xyz,Pi_xyz',initial_force_per_cell=initial['force_per_cell_ry_bohr'],
        original_force_per_cell=initial['original_force_per_cell_ry_bohr'])
    return model, records, identity


class BoundedOutputs:
    def __init__(self, path, capacity):
        self.path, self.capacity, self.used = path.resolve(), capacity, 0
        if self.path == ROOT or ROOT in self.path.parents:
            raise ValueError('saved output arrays must be outside Git/repository')
        if capacity <= 0:
            raise ValueError('positive total output byte capacity required')
        self.path.mkdir(parents=True,exist_ok=False)

    def put(self, name, payload):
        if self.used+len(payload) > self.capacity:
            raise ValueError('explicit output capacity exhausted; earlier completed times retained')
        with (self.path/name).open('xb') as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        self.used += len(payload)
        return binding(self.path/name)

    def json(self, name, value):
        return self.put(name,(json.dumps(pack(value),indent=2,allow_nan=False)+'\n').encode())

    def arrays(self, name, arrays):
        stream = io.BytesIO()
        np.savez(stream,**arrays)
        return self.put(name,stream.getvalue())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full-reference',type=Path,default=HERE/'full-common-reference-004/full-common-reference.json')
    parser.add_argument('--joint',type=Path,default=HERE/'joint-field-aggregate-001/result.json')
    parser.add_argument('--ionic',type=Path,default=HERE/'ionic-gaussian-material-002/result.json')
    parser.add_argument('--oracle',type=Path,required=True)
    parser.add_argument('--times',nargs='+',type=float,required=True,help='explicit times in inverse Ry')
    parser.add_argument('--reference-time',type=float,default=0.,help='second time of ordered connected covariance, inverse Ry')
    parser.add_argument('--kick-cell',nargs=3,type=int,default=(0,0,0))
    parser.add_argument('--kick-atom',type=int,required=True,help='original one-based atom identity')
    parser.add_argument('--kick-axis',choices=('x','y','z'),required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--output-max-bytes',type=int,default=268435456)
    args = parser.parse_args()
    if not np.all(np.isfinite(args.times+[args.reference_time])):
        raise ValueError('all explicit times must be finite')
    sources = source_bindings()
    oracle_binding = binding(args.oracle)
    oracle = bound_json(oracle_binding)
    if (oracle.get('status') != 'passed' or oracle.get('schema') != 'P254-clamped-array-oracle-v1'
            or oracle.get('source_bindings') != sources):
        raise ValueError('independent passed oracle must bind current sources and registration')
    model,inputs,identity = load_material(args.full_reference,args.joint,args.ionic)
    atom_ids = list(identity['original_atoms']['indices_one_based'])
    if args.kick_atom not in atom_ids:
        raise ValueError('kick original atom identity absent')
    coordinate = 3*atom_ids.index(args.kick_atom)+'xyz'.index(args.kick_axis)
    local,common = model.kicks(args.kick_cell,coordinate)
    ix = [coordinate,model.d+coordinate]
    # Output is bounded O(N*d + N*readout^2), not a full real-space covariance.
    out = BoundedOutputs(args.output,args.output_max_bytes)
    out.json('inputs.json',dict(schema='P254-clamped-array-inputs-v1',source_bindings=sources,
        inputs=inputs,oracle=oracle_binding,identity=identity,times_inverse_ry=args.times,
        reference_time_inverse_ry=args.reference_time,kick_cell=args.kick_cell,
        kick_atom_original_one_based=args.kick_atom,kick_axis=args.kick_axis,
        kick_added_energy_ry=model.omega[coordinate],readout_phase_space_indices=ix,
        coefficient_diagnostics=model.coefficient_diagnostics(),
        scope='Full periodic ionic plus ONE q0 field quantum CONTROL; all electrons explicitly held. Not full joint dynamics.',
        covariance='connected ordered; identical in baseline and both coherent-kick states',
        interpretation='Phase locking is mean motion, not evidence for entanglement or a quantum-correlation advantage. No damping/bath/FDT/nuclear prediction/irreversibility.',
        energy_zero='Frozen electronic constant omitted. Absolute total electronic energy not supplied.'))
    origin = model.evolve(0.)
    second = origin if args.reference_time == 0 else model.evolve(args.reference_time)
    initial_energy = {}
    for name,kick in (('baseline',None),('local',local),('common',common)):
        initial_energy[name] = model.diagnostics(origin,model.mean(origin,kick)['q'])['energy_ry']
    completed = []
    for number,time in enumerate(args.times):
        evolution = origin if time == 0 else model.evolve(time)
        baseline = model.mean(evolution)
        arrays = dict(baseline_ions=baseline['ions'],baseline_field=baseline['field'])
        diagnostics = {}
        for name,kick in (('baseline',None),('local',local),('common',common)):
            value = baseline if kick is None else model.mean(evolution,kick)
            if kick is not None:
                arrays[name+'_kick_delta_ions'] = value['ions']-baseline['ions']
                arrays[name+'_kick_delta_field'] = value['field']-baseline['field']
            diagnostic = model.diagnostics(evolution,value['q'])
            diagnostic['energy_drift_ry'] = diagnostic['energy_ry']-initial_energy[name]
            diagnostic['energy_drift_relative'] = abs(diagnostic['energy_drift_ry'])/max(1.,abs(initial_energy[name]))
            diagnostics[name] = diagnostic
        arrays['spatial_equal_time_ordered'] = np.asarray([model.ordered_pair(evolution,evolution,cell,args.kick_cell,ix) for cell in model.indices])
        arrays['spatial_two_time_ordered'] = np.asarray([model.ordered_pair(evolution,second,cell,args.kick_cell,ix) for cell in model.indices])
        arrays['spatial_two_time_commutator'] = np.asarray([arrays['spatial_two_time_ordered'][k]-model.ordered_pair(second,evolution,args.kick_cell,cell,ix).T for k,cell in enumerate(model.indices)])
        arrays['field_equal_time_ordered'] = model.field_ordered_pair(evolution,evolution)
        arrays['field_two_time_ordered'] = model.field_ordered_pair(evolution,second)
        arrays['ion_field_two_time_ordered'] = model.ion_field_ordered_pair(evolution,second,ix)
        archive = out.arrays(f'time-{number:04d}.npz',arrays)
        receipt = out.json(f'time-{number:04d}.json',dict(status='complete',time_inverse_ry=time,
            reference_time_inverse_ry=args.reference_time,archive=archive,diagnostics=diagnostics,
            arrays={k:dict(shape=list(v.shape),dtype=str(v.dtype)) for k,v in arrays.items()}))
        completed.append(receipt)
        print(json.dumps(dict(event='time-complete',receipt=receipt)),flush=True)
    final = out.json('result.json',dict(schema='P254-clamped-array-material-v1',status='complete',
        source_bindings=sources,inputs=inputs,oracle=oracle_binding,completed_times=completed,
        whole_joint_dynamics_complete=False,electrons_held=True))
    print(json.dumps(dict(event='complete',result=final)),flush=True)


if __name__ == '__main__':
    main()
