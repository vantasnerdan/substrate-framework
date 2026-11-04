#!/usr/bin/env python3
"""A39 actual A28 nonzero-q full-PW selected-column current application.

Authenticate original saved A28 endpoint probes and restore the original input
background/bases through the admitted loaders. Apply the new physical operator
and its independent adjoint; persist raw actions, overlaps, and all source
bindings. No DFT job, legacy H action, populated response or trajectory is run.
Root must supply its explicit CPU grant and an external fresh output directory.
"""
import argparse
import json
import os
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0028'),
                str(HERE.parent/'0025'), str(HERE.parent/'0024')]

import numpy as np
import run_full_field_moments_material as full
from material_run_contract import bound_json, runtime_context
from oracle_finiteq_current import binding, implementations
from substrate_framework.material_physical_operators import physical_hamiltonian
from substrate_framework.material_finiteq_current import MaterialFiniteQCurrent, MODEL
from substrate_framework.qe_bloch_edge import BlochEdge

A28 = HERE.parent/'0028'


def require_runtime(cpu_list):
    cpus = {int(part) for part in cpu_list.split(',')}
    actual = set(os.sched_getaffinity(0))
    if not cpus or cpus & {0, 1, 2, 3, 4, 5, 11} or not actual <= cpus:
        raise ValueError('actual affinity must be within Root explicit unprotected CPU grant')
    for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        if os.environ.get(name) != '1':
            raise ValueError(name+'=1 required')


def checked(record):
    if binding(record['path']) != record:
        raise ValueError('source bytes changed: '+record['path'])
    return record


def load_endpoint_bases(grid, metadata):
    bases = []
    for label, ik in (('source', 2), ('target', 3)):
        slot = next(s for s in grid.slots if s.ispin == 1 and s.source_ik == ik
                    and s.operation.index == 1 and not np.any(s.reciprocal_shift))
        frame = grid.load_frame(slot)
        if (frame.basis.identity != metadata[label+'_basis_identity']
                or frame.gauge_identity != metadata[label+'_frame_identity']
                or dict(frame.basis.slot.source_binding) != metadata[label+'_binding']):
            raise ValueError('restored original A28 '+label+' basis/frame differs')
        bases.append(frame.basis)
        del frame
    return bases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--oracle', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--authorized-cpus', required=True, help='Root-owned comma-separated CPU IDs')
    parser.add_argument('--columns', type=int, choices=(1, 2), default=1)
    parser.add_argument('--line-order', type=int, default=12)
    parser.add_argument('--atom-block-size', type=int, default=1)
    parser.add_argument('--array-max-bytes', type=int, default=16*1024**2)
    parser.add_argument('--polarization', type=float, nargs=6,
        default=(.23, .11, -.17, .19, 0., 0.),
        metavar=('X_RE', 'X_IM', 'Y_RE', 'Y_IM', 'Z_RE', 'Z_IM'))
    args = parser.parse_args()
    require_runtime(args.authorized_cpus)
    if args.line_order < 1 or args.atom_block_size < 1 or args.array_max_bytes < 1:
        raise ValueError('positive line order, atom block size and output cap required')
    output = args.output.resolve()
    if output == ROOT or ROOT in output.parents or output.exists():
        raise ValueError('fresh external output directory required; no old banks copied')
    sources = implementations()
    oracle_record = binding(args.oracle)
    oracle = bound_json(oracle_record)
    if (oracle.get('status') != 'passed' or oracle.get('model') != MODEL
            or oracle.get('implementation_sources') != sources):
        raise ValueError('passed A39 focused oracle bound to actual current source required')
    a28_record = binding(A28/'physical-spatial-material-001/result.json')
    if a28_record['sha256'] != 'e5818b11f905f5e281a3518602d153b08f07b8483e35a2ab6877c7c34c024b70':
        raise ValueError('original settled A28 material result bytes changed')
    a28 = bound_json(a28_record)
    if a28.get('status') != 'complete-operator-probes-and-scalar-handoff':
        raise ValueError('complete saved A28 material result required')
    archive_record = checked(a28['archives']['finite-q-witness'])
    # These are ACTUAL saved nonzero-q source/target inputs, not freshly drawn
    # random vectors or a substituted Gamma basis. Select columns, never rows.
    with np.load(archive_record['path'], allow_pickle=False) as bank:
        right = bank['source_probe'][:, :args.columns].copy()
        left = bank['target_probe'][:, :args.columns].copy()
        geometry = {name: bank[name].copy() for name in ('source_miller', 'target_miller',
            'source_k', 'target_k', 'q_fractional', 'q_bohr_inv', 'reciprocal_fold')}
    payload_upper = (16*(2*(right.size+left.size)+4*args.columns**2)
                     +sum(v.nbytes for v in geometry.values())+65536)
    if payload_upper > args.array_max_bytes:
        raise ValueError('bounded action/probe payload exceeds output cap before actions')
    started = time.monotonic()
    full.preload_parsers()
    ref_record, manifest_record, _, _, original, grid = full.load_inputs(
        A28/'full-common-reference-004/full-common-reference.json')
    if ref_record != a28['full_reference'] or manifest_record != a28['manifest']:
        raise ValueError('A28 full reference/manifest changed')
    background, background_receipt = full.restore_accepted(original)
    source, target = load_endpoint_bases(grid, a28['finite_q'])
    edge = BlochEdge.create(background.inputs.reference, source, target, geometry['q_fractional'])
    for name, actual in (('source_miller', source.miller), ('target_miller', target.miller),
                        ('source_k', source.k_cart_bohr_inv), ('target_k', target.k_cart_bohr_inv),
                        ('q_bohr_inv', edge.q_cart_bohr_inv), ('reciprocal_fold', edge.reciprocal_fold)):
        if not np.array_equal(actual, geometry[name]):
            raise ValueError('actual saved A28 endpoint geometry changed: '+name)
    if not np.any(edge.q_cart_bohr_inv) or not np.any(edge.reciprocal_fold):
        raise ValueError('actual nonzero-q, nonzero-fold material edge required')
    sh = physical_hamiltonian(background, source, edge.ispin, block_size=2)
    th = physical_hamiltonian(background, target, edge.ispin, block_size=2)
    owner = MaterialFiniteQCurrent(edge, sh, th, line_order=args.line_order,
                                  block_size=2, atom_block_size=args.atom_block_size)
    a = np.asarray(args.polarization).reshape(3, 2) @ np.array([1., 1j])
    action = owner.prepare(a)
    forward = action.apply(right)
    adjoint = action.adjoint_apply(left)
    forward_overlap = left.conj().T @ forward
    adjoint_overlap = adjoint.conj().T @ right
    defect = float(np.linalg.norm(forward_overlap-adjoint_overlap)
                   /max(1., np.linalg.norm(forward_overlap)))
    if not np.isfinite(defect) or defect > 3e-10:
        raise AssertionError(('actual material true-adjoint pair', defect))
    arrays = dict(source_probe=right, target_probe=left, polarization=a,
        forward_L_ry_bohr=forward, adjoint_L_ry_bohr=adjoint,
        forward_overlap_L_ry_bohr=forward_overlap, adjoint_overlap_L_ry_bohr=adjoint_overlap,
        current_overlap_hartree=-forward_overlap/2, **geometry)
    payload_bytes = sum(value.nbytes for value in arrays.values())
    if payload_bytes+4096*len(arrays) > args.array_max_bytes:
        raise ValueError('bounded output cap exceeded')
    output.mkdir(parents=True, exist_ok=False)
    with (output/'operator-actions.npz').open('xb') as stream:
        np.savez(stream, **arrays)
    archive = binding(output/'operator-actions.npz')
    if int(archive['bytes']) > args.array_max_bytes:
        raise ValueError('actual archive exceeds cap')
    report = dict(schema='P254-A39-finiteq-current-material-v1',
        status='complete-selected-column-operator-smoke', model=MODEL,
        units=dict(field_amplitude='bohr^-1', L='Ry bohr', current='electronic Hartree atomic units; J=-L/2'),
        implementation_sources=sources, oracle=oracle_record,
        original_A28_result=a28_record, saved_A28_endpoint_archive=archive_record,
        original_A28_endpoint_metadata=a28['finite_q'], full_reference=ref_record,
        manifest=manifest_record, background_receipt=background_receipt,
        actual_operator_sources=list(background.inputs.sources),
        source_basis_identity=source.identity, target_basis_identity=target.identity,
        runtime_context=runtime_context(), elapsed_seconds=time.monotonic()-started,
        arrays=archive, payload_bytes=payload_bytes, selected_columns=list(range(args.columns)),
        true_adjoint_relative_defect=defect,
        forward_norm=float(np.linalg.norm(forward)), adjoint_norm=float(np.linalg.norm(adjoint)),
        line_order=args.line_order, atom_block_size=args.atom_block_size, column_block_size=2,
        full_endpoint_dimensions=[edge.source_dimension, edge.target_dimension],
        scope='actual saved A28 full-PW selected-column operator application and independent adjoint only',
        no_hidden_weights=True, populated_response=False, full_source_influence=False,
        quadrature_error_certified=False, physical_transverse_uniqueness_claim=False,
        limitations=['declared straight-line ordinary-NC gauged scalar',
            'finite line quadrature; material refinement not certified by this smoke',
            'selected operator columns do not populate the full response',
            'no intrinsic collision/loss, photon modes, boundary model or whole-E closure'])
    with (output/'result.json').open('x') as stream:
        json.dump(full.pack(report), stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(status=report['status'], result=binding(output/'result.json'))))


if __name__ == '__main__':
    main()
