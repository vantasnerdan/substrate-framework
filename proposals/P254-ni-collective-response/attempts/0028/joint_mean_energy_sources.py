"""Authenticate saved scalar-only energy inputs; never load frames or apply H0/K.

D/frame identities are preserved through authenticated original boundary/request
metadata, not reread. Only the 256 one_body scalars, two JointStages and admitted
case parameter archives are numerical inputs. The one history-owner execution
change is admitted only by the existing current-source numerical proof.
"""
import json
from pathlib import Path

import numpy as np
from substrate_framework.partitioned_joint_slots import array_identity, digest
from substrate_framework.streamed_joint_reduction import JointStage
from substrate_framework.streamed_material_state import authenticate, file_binding

def bound_json(binding):
    return json.loads(authenticate(binding).read_text())


def same_binding(left, right):
    return (left['path'] == right['path'] and left['sha256'] == right['sha256']
            and int(left['bytes']) == int(right['bytes']))


def parameter_identities(model):
    r = model.reference
    return dict(reference={name: array_identity(getattr(r, name)) for name in
        ('nbar', 'masses', 'curvature', 'ionic_force', 'support_force', 'held_density', 'held_gradient')},
        coefficients={name: array_identity(getattr(model.coefficients, name)) for name in
        ('delta_curvature', 'contact', 'mixed')}, charges=array_identity(model.charges),
        polarizations=array_identity(model.field.polarizations), field_mass=model.field.mass,
        field_volume=model.field.volume_bohr3, measure=model.measure, ncells=r.ncells,
        reality_tolerance=model.reality_tolerance)


def read_stage(reg, model, saved, index):
    """Read every genuine scalar in original weighted order, never a residual."""
    source, record = reg['slot_sources'][index], reg['stages'][index]
    request = bound_json(source['request'])
    boundary = bound_json(source['boundary'])
    original = reg['source_model']
    key = digest(request)
    if (request['schema'] != 'partitioned-joint-slot-request-v1' or request['operation'] != 'reduce'
            or request['interval'] != 0 or request['frozen_stage'] is not None
            or request['model_identity'] != original['identity']
            or request['source_identity'] != original['source_identity']
            or Path(source['request']['path']).parent.name != key
            or not same_binding(request['boundary'], source['boundary'])
            or record['boundary'] != source['boundary']['path']
            or request['time'] != source['time'] or request['time'] != record['time_ry_inverse']
            or request['model_parameters'] != source['model_parameters']
            or request['model_parameters'] != parameter_identities(model)):
        raise ValueError('original request/model/physical parameters or stage metadata differs')
    if (boundary['schema'] != 'streamed-joint-mean-v1' or boundary['complete'] is not True
            or boundary['model_identity'] != original['identity']
            or boundary['time_ry_inverse'] != request['time']
            or not same_binding(boundary['original_electronic_boundary'], saved['original_boundary'])
            or boundary['slots'] != request['slots']
            or len(request['slots']) != 128 or len(source['one_body_sources']) != 128
            or tuple(s['identity'] for s in request['slots']) != model.expected_slots):
        raise ValueError('complete original ordered128-slot boundary required')
    initial_boundary = bound_json(reg['slot_sources'][0]['boundary']) if index else boundary
    if any(s['coefficient'] != initial_boundary['slots'][i]['coefficient']
           for i, s in enumerate(request['slots'])):
        raise ValueError('original D bindings changed')
    stage = JointStage.restore(record)
    r = model.reference
    shapes = dict(density=r.nbar.shape, potential=r.nbar.shape,
        gradient=r.masses.shape, field_gradient_electronic=(3,), u=r.masses.shape,
        p=r.masses.shape, coordinate=(model.field.size,), momentum=(model.field.size,), shift=(3,))
    for name, shape in shapes.items():
        value = getattr(stage, name)
        if value.shape != shape or not np.isfinite(value).all():
            raise ValueError('saved stage dimension/nonfinite value: '+name)
    coords = [stage.u, stage.p, stage.coordinate, stage.momentum]
    boundary_coords = [np.asarray(boundary[name], dtype=float) for name in
                       ('u_bohr', 'p_canonical', 'field_coordinate', 'field_momentum')]
    identities = [array_identity(v) for v in coords]
    if (identities != request['reduction_coordinates'] or identities != request['initial_coordinates']
            or any(not np.array_equal(a, b) for a, b in zip(coords, boundary_coords))
            or not np.array_equal(stage.shift, model.field.cartesian_shift(stage.coordinate))):
        raise ValueError('original saved stage/boundary/reduction coordinates differ')
    one_body, ledger = 0j, []
    for i, (slot, entry) in enumerate(zip(request['slots'], source['one_body_sources'])):
        receipt = bound_json(entry['receipt'])
        if (receipt.get('complete') is not True or receipt['request'] != key or receipt['index'] != i
                or receipt['slot'] != slot['identity'] or entry['slot'] != slot['identity']
                or entry['weight'] != slot['weight'] or slot['weight'] != 1/r.ncells
                or slot['component'] != int(slot['identity'].split(':')[0])-1
                or not same_binding(receipt['arrays']['one_body'], entry['one_body'])
                or Path(entry['receipt']['path']) != Path(source['request']['path']).parent/('slot-%04d'%i)/'complete.json'):
            raise ValueError('complete original indexed slot receipt/weight/scalar binding differs')
        value = np.load(authenticate(entry['one_body']), allow_pickle=False)
        if value.shape != () or value.dtype != np.complex128 or not np.isfinite(value).all():
            raise ValueError('genuine finite original complex128 scalar one_body required')
        scalar = value.item()
        one_body += slot['weight']*scalar
        ledger.append(dict(index=i, slot=slot['identity'], weight=slot['weight'],
            receipt=entry['receipt'], one_body=entry['one_body'], value_ry=[scalar.real, scalar.imag]))
    evidence = dict(time_ry_inverse=request['time'], boundary=source['boundary'], request=source['request'],
        request_identity=key, source_identity=request['source_identity'], model_identity=request['model_identity'],
        model_parameters=request['model_parameters'], stage=record, scalar_sources=ledger,
        weighted_one_body_per_cell_ry=[one_body.real, one_body.imag],
        scalar_count=len(ledger), original_D_and_frame_bindings='unchanged in authenticated request/boundary metadata; not read',
        actual_reads='all128 one_body scalars and complete saved JointStage; no other contraction arrays')
    return stage, one_body, evidence
