"""Exact QE 7.5 collinear density replication into a 2x1x1 cell.

This writes only a derived reciprocal density, never native completion metadata.
The source XML uses Hartree energies, Bohr cell/alat, and dimensionless bg;
charge-density.dat stores physical reciprocal columns (tpiba * bg). The two
complex records are total density and magnetization, NOT spin-up/spin-down.
The original QE tolerance-aware heapsort must be linked into the order driver.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET

import numpy as np
from scipy.io import FortranEOFError, FortranFile

_QE_MODULES = Path('/opt/qebuild/q-e-qe-7.5/Modules')


def required_files(ordering_executable: Path | None = None) -> tuple[Path, ...]:
    """Content-binding closure, including the original ordering algorithm."""
    files = (Path(__file__).resolve(), Path(__file__).parent / 'native/qe_density_order.f90',
             _QE_MODULES / 'kind.f90', _QE_MODULES / 'sort.f90',
             _QE_MODULES / 'recvec_subs.f90')
    return files if ordering_executable is None else files + (Path(ordering_executable).resolve(),)


def _binding(path: Path) -> dict:
    data = path.read_bytes()
    return {'path': str(path.resolve()), 'bytes': len(data),
            'sha256': hashlib.sha256(data).hexdigest()}


def _metadata(xml: Path) -> dict:
    root = ET.parse(xml).getroot()
    for node in root.iter():
        node.tag = node.tag.rsplit('}', 1)[-1]
    if root.attrib.get('Units') != 'Hartree atomic units':
        raise ValueError('source XML must declare Hartree atomic units')

    def node(path: str) -> ET.Element:
        result = root.find(path)
        if result is None:
            raise ValueError(f'source XML missing {path}')
        return result

    def number(path: str) -> float:
        return float(node(path).text or '')

    alat = float(node('output/atomic_structure').attrib['alat'])
    bg = np.column_stack([np.fromstring(node(f'output/basis_set/reciprocal_lattice/b{i}').text or '', sep=' ')
                          for i in (1, 2, 3)])
    cell = np.column_stack([np.fromstring(node(f'output/atomic_structure/cell/a{i}').text or '', sep=' ')
                            for i in (1, 2, 3)])
    fft_node = node('output/basis_set/fft_grid')
    fft = tuple(int(fft_node.attrib[f'nr{i}']) for i in (1, 2, 3))
    ecutrho_ry = 2.0 * number('output/basis_set/ecutrho')
    if bg.shape != (3, 3) or cell.shape != (3, 3):
        raise ValueError('source XML reciprocal/cell columns must have three components')
    if not (np.isfinite(bg).all() and np.isfinite(cell).all() and
            np.isfinite(alat) and alat > 0 and np.isfinite(ecutrho_ry) and ecutrho_ry > 0):
        raise ValueError('source XML has invalid lattice/cutoff values')
    if any(n < 2 for n in fft) or abs(np.linalg.det(cell)) == 0:
        raise ValueError('source XML has invalid FFT/cell')
    if (node('output/basis_set/gamma_only').text or '').strip() != 'false':
        raise ValueError('only full reciprocal, non-gamma density is supported')
    if (node('output/magnetization/lsda').text or '').strip() != 'true':
        raise ValueError('source must be collinear spin-polarized density')
    if (node('output/magnetization/noncolin').text or '').strip() != 'false':
        raise ValueError('noncollinear density is not supported')
    tpiba = 2.0 * np.pi / alat
    if not np.allclose(cell.T @ (tpiba * bg), 2.0 * np.pi * np.eye(3), rtol=1e-13, atol=1e-13):
        raise ValueError('source XML cell and normalized reciprocal lattice disagree')
    return {'alat': alat, 'bg': bg, 'fft': fft, 'tpiba': tpiba,
            'ecutrho_ry': ecutrho_ry, 'ngm': int(number('output/basis_set/ngm')),
            'volume': abs(float(np.linalg.det(cell))),
            'nelec': number('output/band_structure/nelec'),
            'magnetization': number('output/magnetization/total')}


def _read_density(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with FortranFile(path, 'r', header_dtype=np.uint32) as stream:
        header = stream.read_ints(np.int32)
        if header.shape != (3,) or header[0] != 0 or header[1] < 1 or header[2] != 2:
            raise ValueError('expected int32 header (gamma_only=false, ngm>0, nspin=2)')
        ngm = int(header[1])
        reciprocal = stream.read_reals(np.float64)
        miller = stream.read_ints(np.int32)
        if reciprocal.size != 9 or miller.size != 3 * ngm:
            raise ValueError('invalid reciprocal/Miller record dimensions')
        fields = np.column_stack([stream.read_record(np.complex128) for _ in range(2)])
        if fields.shape != (ngm, 2) or not np.isfinite(fields).all():
            raise ValueError('invalid total-density/magnetization records')
        try:
            stream.read_record(np.uint8)
        except FortranEOFError:
            pass
        else:
            raise ValueError('unexpected records after collinear density')
    return header, reciprocal.reshape((3, 3), order='F'), miller.reshape((3, ngm), order='F').T, fields


def _order(executable: Path, fft: tuple[int, int, int], tpiba: float,
           ecutrho_ry: float, bg: np.ndarray) -> np.ndarray:
    lines = [' '.join(map(str, fft)), f'{tpiba:.17g} {ecutrho_ry:.17g}']
    lines.extend(' '.join(f'{value:.17g}' for value in bg[:, col]) for col in range(3))
    result = subprocess.run([str(executable.resolve())], input='\n'.join(lines) + '\n',
                            text=True, capture_output=True, check=True)
    tokens = result.stdout.split()
    if not tokens:
        raise ValueError('ordering driver returned no data')
    ngm = int(tokens[0])
    if ngm < 1 or len(tokens) != 1 + 3 * ngm:
        raise ValueError('ordering driver returned malformed Miller records')
    order = np.array([int(value) for value in tokens[1:]], dtype=np.int32).reshape(ngm, 3)
    if len({tuple(row) for row in order}) != ngm:
        raise ValueError('ordering driver returned duplicate Miller indices')
    return order


def fold_density(source_save: Path, target_file: Path, ordering_executable: Path) -> dict:
    """Replicate actual total/magnetization Fourier fields, preserving amplitudes.

    Source order MUST reproduce index-for-index before target enumeration. The
    native reader skips Miller labels when gamma flags match, so this is not
    replaceable by set equality. Any format/order/invariant error raises.
    """
    source_save, target_file = Path(source_save), Path(target_file)
    ordering_executable = Path(ordering_executable)
    density_file, xml_file = source_save / 'charge-density.dat', source_save / 'data-file-schema.xml'
    protected = {path.resolve() for path in required_files(ordering_executable)}
    if (target_file.resolve().is_relative_to(source_save.resolve()) or
            target_file.resolve() in protected):
        raise ValueError('derived output may not overwrite source scientific inputs or ordering code')
    source_bindings = {'density': _binding(density_file), 'xml': _binding(xml_file)}
    meta = _metadata(xml_file)
    header, reciprocal, source_miller, source_fields = _read_density(density_file)
    if int(header[1]) != meta['ngm']:
        raise ValueError('source density count differs from actual XML ngm')
    if not np.allclose(reciprocal, meta['tpiba'] * meta['bg'], rtol=1e-13, atol=1e-13):
        raise ValueError('source density physical reciprocal header differs from actual XML')
    source_order = _order(ordering_executable, meta['fft'], meta['tpiba'], meta['ecutrho_ry'], meta['bg'])
    if not np.array_equal(source_order, source_miller):
        raise ValueError('original QE G order does not reproduce every primitive source index')
    target_fft = (2 * meta['fft'][0], *meta['fft'][1:])
    target_alat = 2.0 * meta['alat']
    target_tpiba = 2.0 * np.pi / target_alat
    target_reciprocal = reciprocal.copy()
    target_reciprocal[:, 0] /= 2.0
    # Normalized bg changes in columns 2,3 as well, because native alat doubles.
    target_bg = meta['bg'].copy()
    target_bg[:, 1:] *= 2.0
    target_miller = _order(ordering_executable, target_fft, target_tpiba, meta['ecutrho_ry'], target_bg)
    source_index = {tuple(row): i for i, row in enumerate(source_miller)}
    target_fields = np.zeros((len(target_miller), 2), dtype=np.complex128)
    copied = set()
    even_count = 0
    for target_index, (h, k, l) in enumerate(target_miller):
        if h % 2 == 0:
            key = (h // 2, k, l)
            if key not in source_index:
                raise ValueError(f'even target Miller index has no primitive source: {key}')
            index = source_index[key]
            target_fields[target_index] = source_fields[index]
            copied.add(index)
            even_count += 1
    if len(copied) != len(source_miller) or even_count != len(source_miller):
        raise ValueError('target even sector does not bijectively preserve primitive density')
    source_zero = source_index.get((0, 0, 0))
    target_zero = np.flatnonzero(np.all(target_miller == 0, axis=1))
    if source_zero is None or len(target_zero) != 1:
        raise ValueError('density must contain exactly one zero Fourier mode')
    source_integrals = meta['volume'] * source_fields[source_zero]
    target_integrals = 2.0 * meta['volume'] * target_fields[target_zero[0]]
    allowance = 1e-10
    declared = np.array([meta['nelec'], meta['magnetization']])
    if not np.all(np.abs(source_integrals - declared) <= allowance * np.maximum(1.0, np.abs(declared))):
        raise ValueError('source density normalization disagrees with completed native XML')
    if not np.array_equal(target_integrals, 2.0 * source_integrals):
        raise ValueError('replicated charge/magnetization integrals do not double')
    target_file.parent.mkdir(parents=True, exist_ok=True)
    with FortranFile(target_file, 'w', header_dtype=np.uint32) as stream:
        stream.write_record(np.array([0, len(target_miller), 2], dtype=np.int32))
        stream.write_record(target_reciprocal.ravel(order='F'))
        stream.write_record(target_miller.T.ravel(order='F').astype(np.int32))
        for spin in range(2):
            stream.write_record(target_fields[:, spin])
    if source_bindings != {'density': _binding(density_file), 'xml': _binding(xml_file)}:
        raise ValueError('source scientific inputs changed during density transformation')
    return {
        'kind': 'qe75-collinear-density-2x1x1-v1',
        'source': {**source_bindings, 'fft': list(meta['fft']), 'ngm': len(source_miller),
                   'alat_bohr': meta['alat'], 'volume_bohr3': meta['volume'],
                   'reciprocal_physical_columns': reciprocal.T.tolist()},
        'target': {'density': _binding(target_file), 'fft': list(target_fft), 'ngm': len(target_miller),
                   'alat_bohr': target_alat, 'volume_bohr3': 2.0 * meta['volume'],
                   'reciprocal_physical_columns': target_reciprocal.T.tolist()},
        'ordering_executable': _binding(ordering_executable),
        'ordering_sources': [_binding(path) for path in required_files()],
        'source_order_equal': True, 'ecutrho_ry': meta['ecutrho_ry'],
        'density_invariants': {
            'fields': ['total_density', 'magnetization'], 'coefficient_scale': 1.0,
            'even_sector_exact_copy': True, 'odd_sector_zero': True,
            'even_sector_count': even_count, 'odd_sector_count': len(target_miller) - even_count,
            'source_integrals': source_integrals.real.tolist(),
            'target_integrals': target_integrals.real.tolist(),
            'source_integral_imaginary': source_integrals.imag.tolist(),
            'target_integral_imaginary': target_integrals.imag.tolist(),
            'declared_source_integrals': declared.tolist(), 'relative_allowance': allowance,
            'source_inputs_unchanged': True,
        },
        'native_completion_claimed': False,
    }
