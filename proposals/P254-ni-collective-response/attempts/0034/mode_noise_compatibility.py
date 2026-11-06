"""Fixed A34 scientific source profile and exact retrospective paid-result admission."""
import hashlib
import json
from pathlib import Path

MANIFEST = Path(__file__).with_name('scientific-compatibility-001.json')


def _binding(path):
    path = Path(path).resolve(strict=True)
    data = path.read_bytes()
    return dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def load_compatibility():
    """Refuse changed scientific code before either production or consumption."""
    with MANIFEST.open(encoding='utf-8') as stream:
        profile = json.load(stream)
    if profile['schema'] != 'P254-A34-scientific-compatibility-v1':
        raise ValueError('fixed A34 scientific compatibility manifest required')
    sources = profile['scientific_sources']
    paths = [item['path'] for item in sources]
    if not sources or paths != sorted(set(paths)):
        raise ValueError('unique sorted scientific source bindings required')
    for expected in sources:
        if _binding(expected['path']) != expected:
            raise ValueError('incompatible scientific source bytes: ' + expected['path'])
    digest = hashlib.sha256(json.dumps(sources, sort_keys=True,
        separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    identity = dict(schema='P254-A34-scientific-source-identity-v1', sha256=digest)
    if identity != profile['scientific_compatibility_identity']:
        raise ValueError('scientific compatibility identity does not bind source bytes')
    return profile


def admit_result(result, result_binding, profile):
    """No unknown pre-identity result or incompatible producer may enter a sum."""
    identity = result.get('scientific_compatibility_identity')
    if identity == profile['scientific_compatibility_identity']:
        return 'producer-scientific-identity'
    if identity is None:
        for admitted in profile['retrospectively_compatible_results']:
            if (admitted['original_record'] == result['original_record']
                    and admitted['result'] == result_binding):
                return 'exact-retrospective-paid-result'
    raise ValueError('result lacks compatible scientific identity or exact paid-result admission')
