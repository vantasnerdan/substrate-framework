"""Reuse t=0 physical contractions only after exact operator/input comparison.

Scalar identity is deliberately NOT reused. Current C, reference and field
coefficients recompute the current stage from the admitted raw slot moments.
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import run_joint_streamed_material as legacy
import matched_reference_io as io
from matched_joint_case import load_case, slot_identity
from substrate_framework.partitioned_joint_slots import array_identity, digest, publish
from substrate_framework.streamed_joint_state import DiskJointMean, JointBudget
from substrate_framework.streamed_material_state import authenticate, file_binding


def bound_json(binding):
    return json.loads(Path(authenticate(binding)).read_text())


def source_values(value):
    """Location is not content; preserve every non-location source field."""
    if isinstance(value, dict):
        return {k: int(v) if k == 'bytes' else source_values(v) for k, v in value.items()
                if not (k == 'path' and 'sha256' in value)}
    if isinstance(value, (list, tuple)):
        return [source_values(v) for v in value]
    return value


def prove(context, current, old_result_binding, old_model_binding, request_binding, budget):
    old_result, old_model, request = (bound_json(b) for b in
                                    (old_result_binding, old_model_binding, request_binding))
    if (old_result['status'] != 'complete' or old_result['mode'] != 'reduce'
            or request['operation'] != 'reduce' or request['interval'] != 0
            or request['frozen_stage'] is not None or old_model['block_size'] != 8):
        raise ValueError('complete original zero-time physical reduction with known block size required')
    old_state = DiskJointMean(authenticate(old_result['initial_boundary']), verify=False)
    if (old_state.time != 0 or current.time != 0 or request['boundary'] != file_binding(old_state.path)
            or request['model_identity'] != old_state.identity
            or tuple(s.identity for s in current.slots) != tuple(s.identity for s in old_state.slots)):
        raise ValueError('complete same ordered original t0 slots required')
    # Authenticate every unchanged original numerical primitive, not the changed
    # former material CLI. Its exercised source bytes are preserved separately.
    primitives = [b for b in old_model['implementation']
                  if Path(b['path']).parent == legacy.ROOT/'src/substrate_framework']
    for binding in primitives:
        authenticate(binding)
    restored, _, _, _, _, _, _, _ = legacy.load_settled()
    old_background, _ = legacy.restore_accepted(restored[4])
    if source_values(old_background.inputs.sources) != source_values(context.background.inputs.sources):
        raise ValueError('original physical input sources differ')
    for name in ('rho_updown', 'vrs'):
        if not np.array_equal(getattr(old_background, name), getattr(context.background, name)):
            raise ValueError('physical saved background differs: '+name)
    old_ref, new_ref = old_background.inputs.reference, context.background.inputs.reference
    for name in ('cell_bohr', 'positions_bohr', 'fft_smooth', 'atom_names'):
        if not np.array_equal(getattr(old_ref, name), getattr(new_ref, name)):
            raise ValueError('original physical geometry/grid differs: '+name)
    if ({name: (table.upf.sha256, table.omega_bohr3) for name, table in old_background.inputs.tables.items()}
            != {name: (table.upf.sha256, table.omega_bohr3) for name, table in context.background.inputs.tables.items()}):
        raise ValueError('physical radial source/volume differs')
    old_grid, old_order = restored[5], restored[6]
    old_slots = {slot_identity(slot): slot for slot in old_order}
    new_slots = {slot_identity(slot): slot for slot in context.ordered}
    records, bases = [], []
    work = Path(request_binding['path']).parent
    key = digest(request)
    for index, (old, new) in enumerate(zip(old_state.slots, current.slots)):
        budget.check()
        for name in ('identity', 'component', 'weight', 'shape', 'source_identity', 'frame', 'coefficient'):
            if getattr(old, name) != getattr(new, name):
                raise ValueError('original contraction input differs: '+old.identity+' '+name)
        a = old_grid.load_frame(old_slots[old.identity]).basis
        b = context.grid.load_frame(new_slots[new.identity]).basis
        fingerprints = {}
        for name in ('miller', 'reciprocal_bohr_inv', 'k_cart_bohr_inv'):
            x, y = np.asarray(getattr(a, name)), np.asarray(getattr(b, name))
            if not np.array_equal(x, y):
                raise ValueError('actual full-PW basis differs: '+old.identity+' '+name)
            fingerprints[name] = array_identity(x)
        record_binding = file_binding(work/('slot-%04d' % index)/'complete.json')
        record = bound_json(record_binding)
        if (record.get('complete') is not True or record['request'] != key
                or record['index'] != index or record['slot'] != old.identity):
            raise ValueError('incomplete or mismatched original contraction')
        for binding in record['arrays'].values():
            authenticate(binding)
        records.append(record_binding)
        bases.append(dict(slot=old.identity, original_basis=old.basis_identity,
                          current_basis=new.basis_identity, exact_numeric_basis=fingerprints))
    return dict(status='complete', schema='P254-matched-initial-contraction-reuse-v1',
        source=file_binding(__file__), context_identity=context.identity, current_boundary=file_binding(current.path),
        old_result=old_result_binding, old_model=old_model_binding, old_request=request_binding,
        unchanged_primitives=primitives, completed_slot_records=records, numeric_basis_comparisons=bases,
        block_size=8, scope='Physical t0 n/G/ell/H0-energy only; no scalar, trajectory or propagated-state equivalence')


def import_initial_contractions(model, state, proof_binding, *, block_size):
    proof = bound_json(proof_binding)
    authenticate(proof['source'])
    if (proof['status'] != 'complete' or proof['schema'] != 'P254-matched-initial-contraction-reuse-v1'
            or proof['context_identity'] != model.identity or proof['block_size'] != block_size
            or proof['current_boundary'] != file_binding(state.path)):
        raise ValueError('source-exact proof does not bind this current initial model/state')
    for binding in proof['unchanged_primitives']:
        authenticate(binding)
    key, directory, _ = model.executor.request(state)
    for index, (slot, binding) in enumerate(zip(state.slots, proof['completed_slot_records'], strict=True)):
        old = bound_json(binding)
        if old['slot'] != slot.identity:
            raise ValueError('reused original slot order differs')
        for array in old['arrays'].values():
            authenticate(array)
        target = directory/('slot-%04d' % index)
        target.mkdir(exist_ok=True)
        if (target/'complete.json').exists():
            model.executor.receipt(directory, key, index, slot)
            continue
        publish(target/'complete.json', dict(complete=True, request=key, index=index, slot=slot.identity,
            cpu=sorted(os.sched_getaffinity(0))[0], blas_threads=1, arrays=old['arrays'],
            hamiltonian_actions=0, action_columns=0, array_output_bytes=0, elapsed_seconds=0.,
            original_completed_work=binding, exact_input_equivalence=proof_binding,
            origin='authenticated physical contraction reuse, not new H0 work'), model.budget)
    return dict(proof=proof_binding, imported_complete_slots=len(state.slots), new_hamiltonian_actions=0)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case', required=True)
    p.add_argument('--initialization', type=Path, required=True)
    p.add_argument('--old-result', type=Path, required=True)
    p.add_argument('--old-model', type=Path, required=True)
    p.add_argument('--old-request', type=Path, required=True)
    p.add_argument('--cpu', type=int, choices=range(6, 11), required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    os.sched_setaffinity(0, {args.cpu})
    io.runtime()
    budget = JointBudget(max_actions=1, max_seconds=600, max_bytes=4*1024**2)
    context = load_case(args.case)
    initialization = bound_json(file_binding(args.initialization))
    if initialization.get('status') != 'complete' or initialization.get('case') != args.case:
        raise ValueError('completed actual initialization required')
    state = DiskJointMean(authenticate(initialization['saved']['joint_boundary']), expected_identity=context.identity)
    result = dict(status='partial', source=file_binding(__file__))
    try:
        result = prove(context, state, file_binding(args.old_result), file_binding(args.old_model),
                       file_binding(args.old_request), budget)
    except Exception as error:
        result['error'] = type(error).__name__+': '+str(error)
    result['budget'] = budget.snapshot()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    if result['status'] != 'complete':
        raise SystemExit(2)


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()
