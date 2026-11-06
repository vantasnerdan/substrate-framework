"""Publish immutable producer output with explicit byte-identical source restoration."""
from pathlib import Path
import json
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent/'0063')]
from history_storage import binding, authenticate, atomic_json
OUT = Path('/var/tmp/p254-a78-fullhistory-common-mean-001')


def main():
    executed = json.loads((OUT/'result.json').read_bytes())
    spec = dict(executed)
    replacements = {
        str(HERE.parent/'0074/fullwindow_force_cache.py'): HERE/'fullwindow_force_cache_original.py.snapshot',
        str(HERE/'produce_fullhistory_mean.py'): HERE/'produce_fullhistory_mean_original.py.snapshot',
        str(HERE/'native_mean.py'): HERE/'native_mean_original.py.snapshot'}
    restored = []
    sources = []
    for original in executed['source_bindings']:
        if original['path'] in replacements:
            archived = binding(replacements[original['path']])
            if any(original[k] != archived[k] for k in ('bytes','sha256')):
                raise ValueError('executed source snapshot differs')
            restored.append(dict(original_binding=original,archived_snapshot=archived))
            sources.append(archived)
        else:
            authenticate(original)
            sources.append(original)
    spec['producer_result'] = binding(OUT/'result.json')
    spec['source_restoration'] = restored
    spec['source_bindings'] = sources
    spec['publication_implementation'] = binding(Path(__file__))
    spec['publication_semantics'] = 'Same actual producer archive/input/material/physical/method/defects; only source restoration paths made explicit after reader cutover. Original producer result remains immutable.'
    atomic_json(OUT/'published-mean.json',spec)
    cache = json.loads(Path(spec['force_cache']['path']).read_bytes())
    records = []
    for original in cache['implementation_bindings']:
        if original['path'] == str(HERE.parent/'0063/saved_mean.py'):
            records.append(dict(original_binding=original,archived_snapshot=spec['producer_reader']['preserved']))
        elif original['path'] == str(HERE.parent/'0074/fullwindow_force_cache.py'):
            records.append(dict(original_binding=original,archived_snapshot=binding(HERE/'fullwindow_force_cache_original.py.snapshot')))
    atomic_json(HERE/'historical-source-snapshots.json',dict(records=records))
    print(json.dumps(dict(published=binding(OUT/'published-mean.json'),snapshots=binding(HERE/'historical-source-snapshots.json'))))


if __name__ == '__main__':
    main()
