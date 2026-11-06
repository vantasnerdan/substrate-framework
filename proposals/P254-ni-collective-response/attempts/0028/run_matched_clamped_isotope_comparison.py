#!/usr/bin/env python3
"""Compare exact H1 and H2 clamped consumers at fixed case/control/potential.

Differences are H2 minus H1. Local/common kicks have equal total energy WITHIN
each isotope, not across isotopes: one local trap quantum changes with mass.
Baseline response/covariance differences do not have that kick-energy confound.
"""
from pathlib import Path
import argparse
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
for directory in (ROOT/'src',HERE,HERE.parent/'0025',HERE.parent/'0024'):
    sys.path.insert(0,str(directory))

from material_run_contract import binding, bound_json
from substrate_framework.material_control_reference import bound_arrays
from run_clamped_array_material import BoundedOutputs
from run_full_field_moments_material import unpack
from matched_clamped_sources import source_bindings


def load_result(path, variant, sources):
    record = binding(path)
    result = bound_json(record)
    if (result.get('schema') != 'P254-matched-clamped-controls-v1'
            or result.get('status') != 'complete' or result.get('mass') != variant
            or result.get('source_bindings') != sources or result.get('electrons_held') is not True):
        raise ValueError('complete exact-source '+variant+' clamped result required')
    inputs = bound_json(result['inputs'])
    if inputs['identity']['mass']['variant'] != variant or inputs['source_bindings'] != sources:
        raise ValueError('isotope input identity differs')
    return record,result,inputs


def paired_inputs(first, second):
    """Only mass-dependent preparation/canonical coupling may differ."""
    for key in ('times_inverse_ry','reference_time_inverse_ry','kick_cell',
                'kick_atom_original_one_based','kick_axis','readout_phase_space_indices','oracle'):
        if first[key] != second[key]:
            raise ValueError('isotope comparison input mismatch: '+key)
    a,b = first['identity'],second['identity']
    for key in ('case','matched_joint_identity','matched_joint_specification','aggregate','archive',
                'control','atom_indices','atom_names','positions_bohr','cell_bohr','valence_charges',
                'cell_indices','extra_delta_curvature'):
        if a[key] != b[key]:
            raise ValueError('fixed physical case identity mismatch: '+key)
    for key in ('nist_source','nist_relative_atomic_masses_amu','original_masses_amu','hydrogen_atom_indices'):
        if a['mass'][key] != b['mass'][key]:
            raise ValueError('isotope source/substitution mismatch: '+key)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--h1',type=Path,required=True)
    parser.add_argument('--h2',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--output-max-bytes',type=int,required=True)
    args = parser.parse_args()
    sources = source_bindings()
    h1,r1,i1 = load_result(args.h1,'H1',sources)
    h2,r2,i2 = load_result(args.h2,'H2',sources)
    paired_inputs(i1,i2)
    if len(r1['completed_times']) != len(r2['completed_times']):
        raise ValueError('complete paired finite times required')
    out = BoundedOutputs(args.output,args.output_max_bytes)
    completed = []
    for number,(a,b) in enumerate(zip(r1['completed_times'],r2['completed_times'])):
        first,second = bound_json(a),bound_json(b)
        for key in ('status','time_inverse_ry','reference_time_inverse_ry','arrays'):
            if first[key] != second[key]:
                raise ValueError('paired time/array semantics differ: '+key)
        x,y = bound_arrays(first['archive']),bound_arrays(second['archive'])
        delta = {name:y[name]-value for name,value in x.items()}
        archive = out.arrays(f'difference-{number:04d}.npz',delta)
        energy = {name:unpack(second['diagnostics'][name]['energy_ry'])-
                  unpack(first['diagnostics'][name]['energy_ry']) for name in ('baseline','local','common')}
        completed.append(out.json(f'difference-{number:04d}.json',dict(
            time_inverse_ry=first['time_inverse_ry'],reference_time_inverse_ry=first['reference_time_inverse_ry'],
            H1=a,H2=b,H2_minus_H1=archive,total_energy_difference_ry=energy)))
    final = out.json('result.json',dict(schema='P254-matched-clamped-isotope-comparison-v1',status='complete',
        source_bindings=sources,H1=h1,H2=h2,case=r1['case'],control=r1['control'],completed_times=completed,
        kick_added_energy_ry=dict(H1=i1['kick_added_energy_ry'],H2=i2['kick_added_energy_ry']),
        scope='Explicit NIST H1 to H2 inertial-mass-only comparison; fixed physical potential, geometry, charge and electronic state. Original standard-weight baseline remains separate.',
        normalization='Same local/common added energy within each isotope, not across isotopes. Baseline and connected covariance differences have no kick-energy change.',
        limits='Not nuclear physics, independent density/resolution convergence, coupled electronic E/D or complete C.'))
    print(json.dumps(dict(event='complete',result=final)),flush=True)


if __name__ == '__main__':
    main()
