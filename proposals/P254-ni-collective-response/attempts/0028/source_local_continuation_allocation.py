"""Bind a wall-only cumulative allocation to the unchanged material owner."""
import argparse
import json
import math
from pathlib import Path

from mixed_material_driver_support import checked_json, source_closure
from run_source_local_q_material import scientific_request
from source_local_material_owner import (
    CAP_NAMES, SCHEMA as OWNER_SCHEMA, execution_caps, isolated_output, material_owner)
from substrate_framework.streamed_material_state import authenticate, file_binding

SCHEMA = 'P254-source-local-material-cumulative-continuation-v1'
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def original_arguments(directory, owner, grant, cpus):
    """Restore selectors and inputs verbatim; execution CPUs are not identity."""
    request = owner['request']
    values = dict(request['selectors'], **request['caps'])
    values.update(output=directory, resume=True, cpus=cpus, frame_admission=[])
    for name, binding in request['explicit_inputs'].items():
        if name == 'frame_admissions':
            values['frame_admission'] = [(item['path'], item['sha256']) for item in binding]
        else:
            values[name] = Path(binding['path'])
            values[name+'_sha256'] = binding['sha256']
    for side in ('left', 'right'):
        binding = grant.get('contexts', {}).get(side+'_context')
        values[side+'_context'] = Path(binding['path']) if binding else None
        values[side+'_context_sha256'] = binding['sha256'] if binding else None
    return argparse.Namespace(**values)


def load_allocation(path, sha256, cpus):
    allocation, grant = checked_json(path, sha256)
    if grant.get('schema') != SCHEMA+'-allocation':
        raise ValueError('explicit cumulative source-local allocation required')
    directory = Path(grant['material_directory']).resolve()
    owner_binding = file_binding(directory/'owner.json')
    owner = json.loads(authenticate(owner_binding).read_text())
    authority = json.loads((directory/'owner-authority.json').read_text())
    previous = json.loads(authenticate(grant['previous_terminal']).read_text())
    if (grant['original_owner'] != owner_binding
            or authority != dict(schema=OWNER_SCHEMA+'-owner-authority', owner=owner_binding)
            or owner.get('schema') != OWNER_SCHEMA+'-owner'
            or previous.get('owner') != owner_binding
            or previous.get('status') != 'failed'
            or previous.get('covariance_executed') is not False
            or grant['original_started_unix'] != owner['started_unix']
            or previous['original_started_unix'] != owner['started_unix']
            or previous['budget']['elapsed_seconds'] < previous['budget']['max_seconds']):
        raise ValueError('unchanged original owner/start and wall-exhausted failed terminal required')
    caps, old = grant['caps'], owner['request']['caps']
    seconds = caps['max_seconds']
    if (set(caps) != set(CAP_NAMES)
            or type(seconds) not in (int, float) or not math.isfinite(seconds)
            or seconds <= max(old['max_seconds'], previous['budget']['max_seconds'])
            or any(caps[name] != old[name] for name in CAP_NAMES if name != 'max_seconds')):
        raise ValueError('wall-only cumulative extension; unchanged action/storage/memory caps required')
    if cpus != grant['cpus']:
        raise ValueError('execution CPUs differ from the bound allocation')
    args = original_arguments(directory, owner, grant, cpus)
    if args.mode != 'pair' or isolated_output(args, ROOT) != directory:
        raise ValueError('same owned original material pair required')
    original_sources = owner['request']['source_bindings']
    for source in original_sources:
        authenticate(source)
    material_owner(directory, scientific_request(args, original_sources, execution_caps(args)),
                   resume=True, started_unix=owner['started_unix'])
    entrypoint = HERE/'run_source_local_material_continuation.py'
    expected = source_closure(original_sources+[file_binding(__file__), file_binding(entrypoint)])
    if grant['source_bindings'] != expected:
        raise ValueError('continuation implementation differs from allocated source bytes')
    return allocation, grant, directory, owner, args, original_sources
