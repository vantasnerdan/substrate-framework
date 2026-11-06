"""Import a material irrep-zero component using its original PH/SCF provenance.

This consumes a material XML source bundle, not the two-atom response exporter.
No SCF, PH, spectrum, stationarity re-admission or screening is performed.
"""
from __future__ import annotations
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
from catdft.checkpoints import fingerprint
from catdft.qe_phonon_inputs import phonon_input
from catdft.qe_scf_xml import node, numbers, text
from .qe_frozen_curvature import parse_frozen_curvature, verify_buffer


def _originals(manifest_path):
    raw = Path(manifest_path).read_bytes()
    manifest = json.loads(raw)
    sources, buffers = manifest['files'], {}
    for name, source in sources.items():
        value = Path(source['path']).read_bytes()
        verify_buffer(value, source)
        buffers[name] = value
    import hashlib
    manifest_binding = {'path': str(Path(manifest_path).resolve()), 'bytes': str(len(raw)),
                        'sha256': hashlib.sha256(raw).hexdigest()}
    return manifest, buffers, manifest_binding


def _reference_chain(manifest, buffers):
    sources = manifest['files']
    producer = json.loads(buffers['producer_spec'])
    snapshot = json.loads(buffers['snapshot'])
    ground = json.loads(buffers['ground_state_result'])
    spec = json.loads(buffers['ground_state_spec'])
    if fingerprint(producer) != snapshot['input_sha256']:
        raise ValueError('original PH spec does not own the Native snapshot')
    if ground['status'] != 'complete' or ground['input_sha256'] != fingerprint(spec):
        raise ValueError('ground-state result does not bind its exact specification')
    for name, key in (('ground_state_result', 'ground_state_result'), ('ground_state_spec', 'ground_state_spec')):
        record = sources[name]
        if producer[key] != record['path'] or producer['input_files_sha256'][record['path']] != record['sha256']:
            raise ValueError('PH producer binds a different ground-state source: '+name)
    artifact = ground['data']['artifacts']['data-file-schema.xml']
    if artifact['path'] != sources['ground_state_xml']['path'] or artifact['sha256'] != sources['ground_state_xml']['sha256']:
        raise ValueError('primary SCF artifact differs from its admitted result')
    if (buffers['ground_state_xml'] != buffers['copied_scf_xml']
            or snapshot['files_sha256']['scf.save/data-file-schema.xml'] != sources['copied_scf_xml']['sha256']):
        raise ValueError('PH snapshot is not the own admitted SCF XML')
    native = Path(producer['directory'])/'native'
    if phonon_input(producer, native, False).encode() != buffers['producer_input']:
        raise ValueError('actual original PH input differs from bound q/masses/settings')
    return producer, ground, spec


def _material_identity(manifest, buffers, manifest_binding, producer, ground, spec):
    root = ET.fromstring(buffers['copied_scf_xml'])
    if root.tag != '{http://www.quantum-espresso.org/ns/qes/qes-1.0}espresso' or root.attrib.get('Units') != 'Hartree atomic units':
        raise ValueError('unexpected material SCF XML units/envelope')
    structure = node(root, 'output/atomic_structure')
    atoms = node(structure, 'atomic_positions').findall('atom')
    nat = int(structure.attrib['nat'])
    symbols = [atom.attrib['name'] for atom in atoms]
    if len(atoms) != nat or symbols != ground['data']['native_xml']['geometry']['symbols']:
        raise ValueError('actual Native atom order differs from primary identity')
    positions = np.array([numbers(atom.text or '', 3) for atom in atoms])
    cell = np.array([numbers(text(structure, 'cell/a'+str(i)), 3) for i in (1, 2, 3)])
    species = node(root, 'output/atomic_species').findall('species')
    names = [item.attrib['name'] for item in species]
    if len(set(names)) != len(names) or set(names) != set(symbols):
        raise ValueError('Native species do not resolve every material atom')
    masses = [producer['inputph']['amass('+str(i)+')'] for i in range(1, len(names)+1)]
    if not np.all(np.isfinite(masses)) or np.any(np.asarray(masses) <= 0):
        raise ValueError('positive actual PH masses required')
    for symbol in names:
        record = manifest['files']['pseudo_'+symbol]
        admitted = ground['data']['native_xml']['species'][symbol]
        if (record['sha256'] != admitted['sha256']
                or record['sha256'] != spec['input_files_sha256'][record['path']]
                or snapshot_digest(buffers, admitted['filename']) != record['sha256']):
            raise ValueError('material pseudopotential source changed')
        header = admitted['header']
        if header['pseudo_type'].strip() != 'NC' or any(header[key].strip().lower() not in ('f', 'false') for key in ('is_ultrasoft', 'is_paw', 'core_correction', 'has_so')):
            raise ValueError('ordinary NC/noNLCC/scalar material scope differs')
    flags = {key: text(root, 'input/spin/'+key) for key in ('lsda', 'noncolin', 'spinorbit')}
    if flags != {'lsda': 'true', 'noncolin': 'false', 'spinorbit': 'false'}:
        raise ValueError('material spin scope differs')
    q = producer['q_cartesian_2pi_alat']
    if producer['inputph']['ldisp'] or q != [0, 0, 0]:
        raise ValueError('registered import is the actual single Gamma point')
    identity = {'nat': nat, 'iq': 1, 'label': manifest['case'], 'q_cart_2pi_alat': q,
                'cell': {'units': 'bohr', 'vectors_rows': cell.tolist(), 'positions_atom_rows': positions.tolist()},
                'bundle_manifest': manifest_binding,
                'bundle_kind': manifest['kind'],
                'applicability_flags': {'ordinary_NC': True, 'selected_PP_NLCC': False, **flags},
                'producer_spec': manifest['files']['producer_spec'], 'snapshot': manifest['files']['snapshot'],
                'ground_state_xml': manifest['files']['ground_state_xml'],
                'atom_symbols': symbols, 'species_in_native_order': names,
                'species_index_one_based': [names.index(symbol)+1 for symbol in symbols],
                'PH_species_mass_amu': masses,
                'PH_atom_mass_amu': [masses[names.index(symbol)] for symbol in symbols],
                'SCF_species_mass_amu': [float(text(item, 'mass')) for item in species],
                'reference_contract': ground['data']['native_xml']['model_contract'],
                'stationarity_scope': 'Existing declared-settings primary admission; joint numerical calibration remains unearned.'}
    return identity


def snapshot_digest(buffers, filename):
    return json.loads(buffers['snapshot'])['files_sha256']['scf.save/'+filename]


def load_material_frozen(manifest_path):
    """Return actual unscreened/unmassweighted component with complete Ni/H identity."""
    manifest, buffers, manifest_binding = _originals(manifest_path)
    producer, ground, spec = _reference_chain(manifest, buffers)
    identity = _material_identity(manifest, buffers, manifest_binding, producer, ground, spec)
    value = parse_frozen_curvature(buffers['frozen_xml'], buffers['patterns_xml'],
        dynmat_source=manifest['files']['frozen_xml'], patterns_source=manifest['files']['patterns_xml'], identity=identity)
    return value, {'manifest': manifest_binding, 'originals': manifest['files']}
