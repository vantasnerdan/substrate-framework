#!/usr/bin/env python3
"""Consume ALL64 authenticated actual selected-q outputs; never propagate q.

Manifest schema is documented in joint-spatial-readout-registration.json.
Each selected_q entry binds either one genuine all-mode result/request, or
separate forward and pair result/request files, all using path/bytes/sha256.
The forward may be a native committed readout with its historical request;
its historical producer and readout consumer identities remain distinct.
Source-local pairs bind the native material result and immutable owner instead
of a legacy request; the original context-pair publication remains authoritative.
Both split originals are preserved in input_bindings, never rewritten as a
combined producer receipt. No complete actual material set exists at
implementation; this CLI is not evidence that such a set has been produced.
"""
import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE)]

import numpy as np
from substrate_framework.streamed_material_state import file_binding
from substrate_framework.joint_spatial_readout import reconstruct
from joint_spatial_receipts import load_material


LIMITS = [
    'One selected original atom/axis force per cell, not all-observable or full-state localization.',
    'Connected AB/BA origins refer to original t0 product preparation, not final subsystems.',
    'Coherent displacement does not change connected covariance within the same linearized Gaussian preparation.',
    'All q are explicit; no positive-time symmetry completion, clipping or reality repair.',
    'Finite original BvK periodic volume: minimum-image extent cannot distinguish image overlap or prove infinite-volume localization.',
    'One finite time pair inside the authenticated baseline window: no lifetime, asymptotic decay or late-time persistence inference.',
    'No entanglement inference from coherent mean locking or selected covariance alone.',
    'No local energy density/current/flow inferred from force; optional mean energies concern only the exact unperturbed homogeneous baseline boundary.',
    'Selected q producer accuracy is not certified; reconstruction does not certify temporal, basis or physical model accuracy.',
    'Existing legacy clamped result uses a different scalar/preparation and is not automatically a matched comparison.'
]


def sources():
    paths = [HERE/name for name in ('joint-spatial-readout-registration.json',
        'joint_spatial_receipts.py', 'joint_spatial_provenance.py', 'joint_spatial_source_local.py',
        'joint_spatial_receipt_oracle.py', 'joint_spatial_committed_oracle.py',
        'joint_spatial_source_local_oracle.py', 'joint_spatial_dual_coherent_oracle.py',
        'joint_spatial_dual_coherent.py', 'joint_spatial_dual_authority.py', 'joint_spatial_pair_identity.py',
        'run_source_local_q_oracle.py', 'run_joint_response_oracle.py', 'run_joint_partitioned_q_oracle.py',
        'run_source_local_q_material.py', 'source_local_material_owner.py',
        'source_local_material_runtime.py', 'one-sided-material-original-source-snapshots-001.json',
        'run_joint_spatial_readout.py', 'run_joint_spatial_readout_oracle.py')]
    paths += [ROOT/'src/substrate_framework'/name for name in (
        'joint_spatial_readout.py', 'spatial_readout_metrics.py', 'joint_response_observables.py',
        'material_control_reference.py', 'streamed_material_state.py', 'qe_saved_native.py',
        'partitioned_joint_slots.py')]
    return [file_binding(path) for path in paths]


def encode(value):
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    if isinstance(value, np.ndarray):
        return encode(value.tolist())
    if isinstance(value, (complex, np.complexfloating)):
        return dict(real=float(value.real), imag=float(value.imag))
    if isinstance(value, np.generic):
        return value.item()
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--manifest-sha256', required=True)
    parser.add_argument('--reference-cell', nargs=3, type=int,
                        help='Covariance right-cell; default is the authenticated coherent kick cell')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    binding = file_binding(args.manifest)
    if binding['sha256'] != args.manifest_sha256:
        raise ValueError('explicit spatial manifest digest differs')
    data = load_material(binding)
    common = data['common']
    reference = args.reference_cell if args.reference_cell is not None else common['cell']
    readout = reconstruct(data['q'], data['coherent'], data['ab'], data['ba'],
        mesh=data['mesh'], cell_bohr=data['geometry']['cell_bohr'],
        kick_cell=common['cell'], reference_cell=reference,
        probe_position_bohr=data['geometry']['atom_position_bohr'],
        seed_position_bohr=data['geometry']['seed_atom_position_bohr'])
    bindings = sources()
    result = dict(schema='P254-joint-spatial-readout-v1', status='complete',
        input_manifest=binding, input_bindings=data['manifest'], source_bindings=bindings,
        receipt_assembly_by_q=data['receipt_assembly_by_q'],
        source_file_count=len(bindings), case=common['case'], model_identity=common['model_identity'],
        initialization=data['manifest']['initialization'], original_preparation=data['saved']['preparation'],
        scalar=data['saved']['specification']['scalar'], baseline=data['baseline'],
        baseline_boundary_time_ry_inverse=data['baseline_time'], dynamics=common['dynamics'],
        observable=common['observable'], time_ry_inverse=common['times'][-1],
        second_time_ry_inverse=common['second_times'][-1],
        time_grid_ry_inverse=common['times'], second_time_grid_ry_inverse=common['second_times'],
        coherent_kick=dict(amplitude_canonical_momentum=common['amplitude'],
                           atom_one_based=common['seed_atom'], axis=common['seed_axis'], cell=common['cell']),
        geometry=data['geometry'], q_fractional=data['q'], readout=readout,
        normalization=dict(coherent='sum exp(+2pi i q.R)*Fq/sqrt(N); seed phase/amplitude already included',
            connected='sum exp(+2pi i q.(R-Rref))*Cq/N; AB/BA already opposite-q paired; no extra dagger'),
        mean_diagnostics=data['mean_diagnostics'], producer_diagnostics=data['producer_diagnostics'],
        numerical_accuracy_certified=False, scientific_limits=LIMITS)
    with args.output.open('x') as stream:
        json.dump(encode(result), stream, indent=2, allow_nan=False)
        stream.write('\n')


if __name__ == '__main__':
    main()
