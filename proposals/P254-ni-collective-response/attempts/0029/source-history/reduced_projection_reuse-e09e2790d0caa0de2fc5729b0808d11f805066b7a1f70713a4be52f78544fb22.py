"""Explicit provenance-preserving reuse of unchanged Gamma pilot items.

Never modifies/rebinds the original manifest. Historical producer bytes must be
available as exact archived source snapshots. Only original initial totals,
mechanical coefficients and the basis-independent Gamma kernel are reused.
"""
from pathlib import Path
import json
from substrate_framework.streamed_material_state import authenticate, file_binding


def admit_gamma_pilot(manifest_path, current_spec, source_history):
    binding = file_binding(manifest_path)
    old = json.loads(authenticate(binding).read_text())
    spec = old['specification']
    for name in ('case', 'case_identity', 'context', 'original_boundary', 'neutral_labels',
                 'bands', 'features', 'mesh', 'ndof', 'indexing'):
        if spec.get(name) != current_spec.get(name):
            raise ValueError('Gamma pilot reuse authority differs: '+name)
    if spec['feature_g'][0] != current_spec['feature_g'][0] or spec['feature_partner'][0] != current_spec['feature_partner'][0]:
        raise ValueError('Gamma pilot feature convention differs')
    snapshots = []
    history = Path(source_history)
    for original in spec['implementation']:
        original_path = Path(original['path'])
        snapshot_path = history/(original_path.stem+'-'+original['sha256']+original_path.suffix)
        snapshot = file_binding(snapshot_path)
        if snapshot['sha256'] != original['sha256'] or int(snapshot['bytes']) != int(original['bytes']):
            raise ValueError('historical producer snapshot bytes differ')
        snapshots.append(dict(original_binding=original, archived_snapshot=snapshot))
    for record in spec['physical_implementation']:
        authenticate(record)
    keys = ['initial', 'mechanical', 'kernel/00']
    adopted = {}
    for key in keys:
        if key not in old['items']:
            continue
        item = old['items'][key]
        authenticate(item['archive'])
        adopted[key] = dict(archive=item['archive'], metadata=dict(item['metadata'],
            explicit_reuse=dict(manifest=binding, item=key, original_producer_snapshots=snapshots,
                scope='basis-independent initial/mechanical/Gamma-kernel product only; original manifest unchanged')))
    return binding, adopted, {k: old.get('timings', {}).get(k, {}) for k in adopted}


def adopt_items(cache, admitted):
    """Publish references to original immutable payloads, never copy/overwrite them."""
    _, items, timings = admitted
    changed = False
    for key, item in items.items():
        if key in cache.manifest['items']:
            if cache.manifest['items'][key]['archive'] != item['archive']:
                raise ValueError('explicit reuse collides with a different completed item: '+key)
            continue
        cache.manifest['items'][key] = item
        cache.manifest['timings'][key] = dict(timings[key], adopted_without_operator_recomputation=True)
        changed = True
    if changed:
        cache.publish()
