"""Retain producer-bound source bytes independently of mutable live filenames.

The archive is explicitly retrospective recovery of the owner's bound bytes,
not a claim that the producer created an execution snapshot. It never supplies
Python code to import: the consumer separately binds its current implementation.
"""
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from substrate_framework.streamed_material_state import authenticate, file_binding

SCHEMA = 'P254-original-producer-source-archive-v1'


def _matches(data, binding):
    return len(data) == int(binding['bytes']) and sha256(data).hexdigest() == binding['sha256']


def bank_producer_sources(owner_binding, destination):
    """Create a NEW archive only when every original source still matches."""
    owner_bytes = authenticate(owner_binding).read_bytes()
    owner = json.loads(owner_bytes)
    sources = owner['request']['source_bindings']
    if (owner['request'].get('schema') != 'P254-native-mixed-q-material-v1-owner'
            or len({s['path'] for s in sources}) != len(sources)):
        raise ValueError('original Native owner and unique source paths required')
    entries = []
    destination = Path(destination)
    with ZipFile(destination, 'x', compression=ZIP_DEFLATED) as archive:
        archive.writestr('owner.json', owner_bytes)
        for index, binding in enumerate(sources):
            data = Path(binding['path']).read_bytes()
            if not _matches(data, binding):
                raise ValueError('cannot recover changed original source: '+binding['path'])
            member = 'sources/%04d'%index
            archive.writestr(member, data)
            entries.append(dict(original=binding, member=member))
        manifest = dict(schema=SCHEMA, recovered_at=datetime.now(timezone.utc).isoformat(),
            recovery_kind='retrospective exact owner-bound source bytes, not original execution snapshot',
            original_owner=dict(owner_binding), sources=entries)
        archive.writestr('manifest.json', json.dumps(manifest, indent=2)+'\n')
    binding = file_binding(destination)
    read_producer_sources(binding, owner_binding, sources)
    return binding


def read_producer_sources(archive_binding, owner_binding, expected_sources):
    """Authenticate the complete exact source inventory; return original bytes."""
    archive_path = authenticate(archive_binding)
    with ZipFile(archive_path) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)):
            raise ValueError('ambiguous duplicate source archive member')
        manifest = json.loads(archive.read('manifest.json'))
        if (manifest.get('schema') != SCHEMA
                or manifest.get('original_owner') != owner_binding
                or manifest.get('recovery_kind') !=
                    'retrospective exact owner-bound source bytes, not original execution snapshot'):
            raise ValueError('source archive owner/recovery authority differs')
        owner_bytes = archive.read('owner.json')
        if not _matches(owner_bytes, owner_binding):
            raise ValueError('archived original owner bytes differ')
        owner = json.loads(owner_bytes)
        if (owner['request'].get('schema') != 'P254-native-mixed-q-material-v1-owner'
                or owner['request']['source_bindings'] != expected_sources
                or [entry['original'] for entry in manifest['sources']] != expected_sources):
            raise ValueError('source archive differs from complete original typed inventory')
        members = [entry['member'] for entry in manifest['sources']]
        if (len(set(members)) != len(members)
                or set(names) != {'manifest.json', 'owner.json', *members}
                or len({s['path'] for s in expected_sources}) != len(expected_sources)):
            raise ValueError('source archive member coverage differs')
        result = {}
        for entry in manifest['sources']:
            original, member = entry['original'], entry['member']
            if archive.getinfo(member).file_size != int(original['bytes']):
                raise ValueError('archived original source size differs: '+original['path'])
            data = archive.read(member)
            if not _matches(data, original):
                raise ValueError('archived original source bytes differ: '+original['path'])
            result[original['path']] = data
    return result
