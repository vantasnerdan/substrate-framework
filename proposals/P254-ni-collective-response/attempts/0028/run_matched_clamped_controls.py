#!/usr/bin/env python3
"""Actual registered case clamped-electron controls; not coupled electronic E/D."""
from pathlib import Path
import argparse
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
for directory in (ROOT/'src', HERE, HERE.parent/'0025', HERE.parent/'0024'):
    sys.path.insert(0,str(directory))

import numpy as np
from material_run_contract import bound_json, binding
from run_clamped_array_material import BoundedOutputs
from matched_clamped_sources import REGISTRATION, require_oracle
from matched_clamped_case import load_clamped
from matched_clamped_readouts import time_readout


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    registration = bound_json(binding(REGISTRATION))
    parser.add_argument('--case',choices=tuple(registration['cases']),required=True)
    parser.add_argument('--mass',choices=('original','H1','H2'),required=True)
    parser.add_argument('--control',choices=('periodic','independent-regions'),required=True)
    parser.add_argument('--oracle',type=Path,required=True)
    parser.add_argument('--times',type=float,nargs='+',required=True,help='finite explicit inverse-Ry times')
    parser.add_argument('--reference-time',type=float,required=True,help='second ordered-covariance time, inverse Ry')
    parser.add_argument('--kick-cell',type=int,nargs=3,default=(0,0,0))
    parser.add_argument('--kick-atom',type=int,required=True,help='original one-based atom identity')
    parser.add_argument('--kick-axis',choices=('x','y','z'),required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--output-max-bytes',type=int,required=True)
    args = parser.parse_args()
    if not np.all(np.isfinite(args.times+[args.reference_time])):
        raise ValueError('all times must be finite')
    return args, registration


def main():
    args, registration = arguments()
    oracle, sources = require_oracle(args.oracle)
    model, identity = load_clamped(args.case,args.mass,args.control)
    atom_ids = list(identity['atom_indices'])
    if args.kick_atom not in atom_ids:
        raise ValueError('original atom identity absent from this case')
    coordinate = 3*atom_ids.index(args.kick_atom)+'xyz'.index(args.kick_axis)
    local, common = model.kicks(args.kick_cell,coordinate)
    kicks = dict(baseline=None,local=local,common=common)
    origin = model.evolve(0.)
    second = origin if args.reference_time == 0 else model.evolve(args.reference_time)
    initial_energy = {name:model.diagnostics(origin,model.mean(origin,kick)['q'])['energy_ry']
                      for name,kick in kicks.items()}
    out = BoundedOutputs(args.output,args.output_max_bytes)
    inputs = out.json('inputs.json',dict(schema='P254-matched-clamped-inputs-v1',
        source_bindings=sources,oracle=oracle,identity=identity,
        times_inverse_ry=args.times,reference_time_inverse_ry=args.reference_time,
        kick_cell=args.kick_cell,kick_atom_original_one_based=args.kick_atom,kick_axis=args.kick_axis,
        kick_added_energy_ry=model.omega[coordinate],initial_total_energy_ry=initial_energy,
        readout_phase_space_indices=[coordinate,model.d+coordinate],
        coefficient_diagnostics=model.coefficient_diagnostics(),
        field_index='singleton global field' if args.control == 'periodic' else 'original cell order, one independent field per region',
        field_covariance_pair='field[index](t), field at kick region(s); singleton for global field',
        ion_field_covariance_pair='selected ion[cell](t), field at kick region(s); same global field for periodic',
        covariance='connected ordered; identical for baseline and both coherent kicks',
        construction=registration['construction'],limits=registration['limits'],
        interpretation='Mean phase locking is not entanglement. Field multiplicity/preparation changes explicitly in independent-regions control.'))
    completed = []
    for number,time in enumerate(args.times):
        left = origin if time == 0 else model.evolve(time)
        arrays,diagnostics = time_readout(model,left,second,args.kick_cell,coordinate,kicks,initial_energy)
        archive = out.arrays(f'time-{number:04d}.npz',arrays)
        receipt = out.json(f'time-{number:04d}.json',dict(status='complete',time_inverse_ry=time,
            reference_time_inverse_ry=args.reference_time,archive=archive,diagnostics=diagnostics,
            arrays={name:dict(shape=list(value.shape),dtype=str(value.dtype)) for name,value in arrays.items()}))
        completed.append(receipt)
        print(json.dumps(dict(event='time-complete',receipt=receipt)),flush=True)
    final = out.json('result.json',dict(schema='P254-matched-clamped-controls-v1',status='complete',
        source_bindings=sources,oracle=oracle,inputs=inputs,completed_times=completed,
        case=args.case,mass=args.mass,control=args.control,electrons_held=True,
        whole_joint_dynamics_complete=False,controlled_lattice_comparison_complete=False))
    print(json.dumps(dict(event='complete',result=final)),flush=True)


if __name__ == '__main__':
    main()
