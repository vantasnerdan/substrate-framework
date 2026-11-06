"""Bind complete NC exports to native patterns and same-gauge ordinary elph XML."""
from pathlib import Path

import numpy as np

from . import qe_phonon_io as xml
from .qe_io import sha256_file
from .qe_phonon_vertex_decode import open_vertices
from . import qe_phonon_vertex_fermi as fermi
from .qe_phonon_vertex_gauge import screened_hamiltonian_matrix


def validate_request(spec):
    if 'raw_vertices' not in spec:
        return
    request = spec['raw_vertices']
    if not isinstance(request, dict) or set(request) != {'build_receipt', 'cpu_affinity'}:
        raise ValueError('raw_vertices requires a private build_receipt and explicit singleton cpu_affinity')
    cpus = request['cpu_affinity']
    if not isinstance(cpus, list) or len(cpus) != 1 or type(cpus[0]) is not int or cpus[0] < 0:
        raise ValueError('Raw NC vertices require exactly one explicit nonnegative CPU index')
    path = Path(request['build_receipt'])
    if not path.is_absolute():
        raise ValueError('Raw vertex build receipt must be absolute')
    execution = spec['execution']
    if execution['mpi_ranks'] != 1 or execution['kpoint_pools'] != 1:
        raise ValueError('Raw NC vertices support exactly one rank and one pool')
    if (spec['inputph']['ldisp'] or not spec['inputph']['trans']
            or spec['inputph'].get('electron_phonon') != 'epa'):
        raise ValueError('Raw NC vertices require inline single-q trans=true ordinary EPA publication')
    if not spec['inputph'].get('fildvscf'):
        raise ValueError('Raw NC vertices require saved native fildvscf')


def source_paths(spec):
    if 'raw_vertices' not in spec:
        return []
    from . import qe_phonon_vertex_decode, qe_phonon_vertex_gauge, qe_phonon_vertex_recovery, qe_phonon_vertex_launch, qe_runtime
    return [Path(__file__), Path(qe_phonon_vertex_decode.__file__),
            Path(qe_phonon_vertex_gauge.__file__),
            Path(fermi.__file__),
            Path(qe_phonon_vertex_recovery.__file__), Path(qe_phonon_vertex_launch.__file__),
            Path(qe_runtime.__file__),
            Path(qe_phonon_vertex_launch.taskset_path()), Path(spec['raw_vertices']['build_receipt']),
            *qe_runtime.runtime_paths([spec['execution']['ph_executable']])]


def validate_build(spec):
    import json
    from .checkpoints import fingerprint
    path = Path(spec['raw_vertices']['build_receipt'])
    receipt = json.loads(path.read_text())
    if receipt.get('status') != 'complete' or receipt.get('schema_version') != 1:
        raise ValueError('Raw vertex private binary is not a completed bound build')
    if receipt['executable']['path'] != spec['execution']['ph_executable']:
        raise ValueError('Raw vertex executable differs from build receipt')
    roles = receipt['roles']
    if set(roles) != {'native_source', 'native_main', 'native_ph_library', 'native_flags',
                      'patcher', 'writer', 'native_recovery', 'private_source',
                      'native_fermi_source', 'private_fermi_source', 'private_fermi_object'} or not receipt['inputs']:
        raise ValueError('Raw vertex receipt lacks required source/object roles')
    if any(record not in receipt['inputs'] for record in roles.values()):
        raise ValueError('Raw vertex source roles are outside its bound build inputs')
    if receipt['actual_affinity'] != [11]:
        raise ValueError('Raw vertex build was not isolated to reserved CPU11')
    for record in [receipt['executable'], *receipt['inputs']]:
        if sha256_file(record['path']) != record['sha256']:
            raise ValueError(f'Raw vertex build content changed: {record["path"]}')
    return fingerprint(receipt)

def _k_metadata(record, nk, q, reciprocal):
    gamma = np.linalg.norm(q) < 1e-12
    ik = record['k_index']
    ikk = ik if gamma else 2 * ik - 1
    ikq = ikk if gamma else ikk + 1
    spin = 1 if ik <= nk // 2 else 2
    if nk % 2 or (record['ikk'], record['ikq'], record['spin']) != (ikk, ikq, spin):
        raise ValueError('Raw vertex native full-mesh k/spin indexing differs')
    difference = np.array(record['kq_cartesian_2pi_alat']) - record['k_cartesian_2pi_alat'] - q
    fractional = np.linalg.solve(reciprocal, difference)
    if not np.allclose(fractional, np.rint(fractional), rtol=0, atol=1e-10):
        raise ValueError('Raw vertex k+q coordinates are not congruent with k and q')



def _validate_irrep(path, native, iq, irr, count, offset, nat, q):
    if path.is_symlink() or not path.is_file():
        raise ValueError('Raw vertex scientific files must be ordinary bound files')
    phsave = native / 'scratch' / '_ph0' / 'scf.phsave'
    patterns = xml.element(xml.read_xml(phsave / f'patterns.{iq}.xml'), 'IRREPS_INFO')
    representation = xml.element(patterns, f'REPRESENTION.{irr}')
    elph = xml.element(xml.read_xml(phsave / f'elph.{iq}.{irr}.xml'), 'PARTIAL_EL_PHON')
    nk, nb = xml.integer(elph, 'NUMBER_OF_K'), xml.integer(elph, 'NUMBER_OF_BANDS')
    maximum_error, records_count = 0.0, 0
    hamiltonian_maximum_error = 0.0
    diagonal = np.diag_indices(nb)
    geometry = xml.element(xml.read_xml(native / 'dyn.xml'), 'GEOMETRY_INFO')
    reciprocal = xml.values(geometry, 'BG', 9).reshape((3, 3), order='F')
    alat = xml.values(geometry, 'CELL_DIMENSIONS', 6)[0]
    metadata = {}
    scf = xml.read_xml(native / 'scratch' / 'scf.save' / 'data-file-schema.xml')
    metallic = scf.findtext('./input/bands/occupations') == 'smearing'
    expected_shift = metallic and np.array_equal(q, np.zeros(3))
    with open_vertices(path) as (header, records):
        required = {'q_index': iq, 'irrep': irr, 'perturbations': count,
                    'mode_offset': offset, 'bands': nb, 'kpoints': nk, 'atoms': nat}
        if any(header[key] != value for key, value in required.items()):
            raise ValueError('Raw vertex dimensions/indices disagree with native XML')
        if not np.allclose(header['q_cartesian_2pi_alat'], q, rtol=0, atol=1e-12):
            raise ValueError('Raw vertex q differs from requested representative')
        if not np.isclose(header['alat_bohr'], alat, rtol=0, atol=1e-12):
            raise ValueError('Raw/native lattice scale differs')
        if header['fermi_shift_applied'] != expected_shift:
            raise ValueError('Raw native Fermi-shift branch differs from SCF occupations and exact Gamma')
        fermi_source, native_shifts = fermi.validate(native, header)
        for ipert in range(1, count + 1):
            pattern = xml.complex_values(xml.element(representation, f'PERTURBATION.{ipert}'),
                                         'DISPLACEMENT_PATTERN', (3 * nat,))
            if not np.allclose(header['patterns'][:, ipert - 1], pattern, rtol=0, atol=1e-12):
                raise ValueError('Raw vertex native displacement pattern differs')
        for record in records:
            records_count += 1
            _k_metadata(record, nk, q, reciprocal)
            signature = [record['ikk'], record['ikq'], record['spin'],
                         record['k_cartesian_2pi_alat'], record['kq_cartesian_2pi_alat'],
                         record['bands'].tolist()]
            if record['k_index'] in metadata and metadata[record['k_index']] != signature:
                raise ValueError('Raw band/spin/occupation metadata changes across perturbations')
            metadata[record['k_index']] = signature
            knode = xml.element(elph, f'K_POINT.{record["k_index"]}')
            coordinates = xml.values(knode, 'COORDINATES_XK', 3)
            if not np.allclose(record['k_cartesian_2pi_alat'], coordinates, rtol=0, atol=1e-12):
                raise ValueError('Raw/native k coordinates disagree')
            if record['kind'] == 'SCREENED_NATIVE':
                vertices = knode.findall('PARTIAL_ELPH')
                node = vertices[record['perturbation'] - 1]
                if int(node.attrib['perturbation']) != record['perturbation']:
                    raise ValueError('Native perturbation metadata differs')
                pairs = xml.numbers(node.text, 2 * nb * nb).reshape(-1, 2)
                reference = (pairs[:, 0] + 1j * pairs[:, 1]).reshape((nb, nb), order='F')
                error = float(np.max(np.abs(record['matrix'] - reference)))
                maximum_error = max(error, maximum_error)
                if not np.allclose(record['matrix'], reference, rtol=2e-12, atol=2e-13):
                    raise ValueError(f'Same-gauge native screened vertex mismatch: {error}')
                absolute = screened_hamiltonian_matrix(header, record)
                reference[diagonal] += native_shifts[record['perturbation'] - 1]
                error = float(np.max(np.abs(absolute - reference)))
                hamiltonian_maximum_error = max(hamiltonian_maximum_error, error)
                if not np.allclose(absolute, reference, rtol=2e-12, atol=2e-13):
                    raise ValueError(f'Absolute Hamiltonian versus native EPA plus exact native shift mismatch: {error}')
    return {'path': str(path), 'sha256': sha256_file(path), **required,
            'q_cartesian_2pi_alat': header['q_cartesian_2pi_alat'],
            'alat_bohr': header['alat_bohr'], 'complex_matrix_records': records_count,
            'screened_native_max_abs_error_ry_per_bohr': maximum_error,
            'screened_hamiltonian_native_plus_shift_max_abs_error_ry_per_bohr': hamiltonian_maximum_error,
            'fermi_shift_applied': header['fermi_shift_applied'],
            'native_attempt': header['native_attempt'], 'native_fermi_shift_evidence': fermi_source,
            'fermi_shift_ry_per_bohr': [[float(z.real), float(z.imag)]
                                      for z in header['fermi_shift_ry_per_bohr']]}


def validate_exports(spec, native, nat):
    if 'raw_vertices' not in spec:
        return None, []
    build = validate_build(spec)
    phsave = native / 'scratch' / '_ph0' / 'scf.phsave'
    info = xml.element(xml.read_xml(phsave / 'patterns.1.xml'), 'IRREPS_INFO')
    nirr = xml.integer(info, 'NUMBER_IRR_REP')
    files, offset = [], 0
    folder = native / 'raw-vertices'
    wanted = {folder / f'q1-irrep{irr}.raw' for irr in range(1, nirr + 1)}
    if set(folder.glob('*')) != wanted:
        raise ValueError('NC raw vertex publication is incomplete or contains foreign/partial records')
    for irr in range(1, nirr + 1):
        count = xml.integer(xml.element(info, f'REPRESENTION.{irr}'), 'NUMBER_OF_PERTURBATIONS')
        files.append(_validate_irrep(folder / f'q1-irrep{irr}.raw', native, 1, irr, count,
                                     offset, nat, spec['q_cartesian_2pi_alat']))
        offset += count
    if offset != 3 * nat:
        raise ValueError('Raw NC exports do not span all atomic displacement modes')
    return {'format': 'CATDFT_NC_VERTEX_V2', 'build_sha256': build, 'irreps': files,
            'units': 'Ry/bohr', 'matrix_axes': ['bra_band_k_plus_q', 'ket_band_k'],
            'basis_gauge': 'native unsquared complex, no band rotation, no mass weighting',
            'streamed_operators': {
                'BARE': 'native full dvqpsi_us before induced HXC; local and both projector-motion terms',
                'SCREENED_NATIVE': 'same-gauge native ordinary/EPA; metallic Gamma is d(H-mu I)'},
            'on_demand_views': {
                'screened_hamiltonian_matrix': 'absolute dH from native plus matching exact def*identity only in applied branch',
                'screened_hamiltonian_diagonal': 'same absolute diagonal without another dense matrix'},
            'topology': {'mpi_ranks': 1, 'kpoint_pools': 1}}, sorted(wanted) + fermi.files(native)
