"""Run only after Root links the original QE kind/sort with the order driver.

Usage: python oracle_density.py --ordering-executable /absolute/qe_density_order
Optional --source-save and --target-file select real inputs/derived output.
No Native calculation or completed-save XML is created. The consumer below
independently reads actual files and reconstructs real-space fields by IFFT.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import xml.etree.ElementTree as ET

import numpy as np
from scipy.io import FortranEOFError, FortranFile

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / 'src'))
from substrate_framework.qe_density_fold import fold_density

RELATIVE_ALLOWANCE = 1e-10


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _consume(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    # Intentionally independent of the production module's binary reader.
    with FortranFile(path, 'r', header_dtype=np.uint32) as stream:
        header = stream.read_record(np.int32)
        _require(header.shape == (3,) and header[0] == 0 and header[2] == 2,
                 'oracle requires full-G, two-field int32 density header')
        count = int(header[1])
        reciprocal_values = stream.read_record(np.float64)
        miller_values = stream.read_record(np.int32)
        _require(reciprocal_values.size == 9 and miller_values.size == count * 3,
                 'oracle reciprocal/Miller record size mismatch')
        reciprocal = reciprocal_values.reshape(3, 3, order='F')
        miller = miller_values.reshape(count, 3)
        columns = [stream.read_record(np.complex128) for _ in range(2)]
        _require(all(column.shape == (count,) for column in columns), 'oracle field count mismatch')
        fields = np.stack(columns, axis=1)
        try:
            stream.read_record(np.uint8)
        except FortranEOFError:
            pass
        else:
            raise ValueError('oracle found unexpected density records')
    _require(np.isfinite(fields).all(), 'oracle found non-finite density')
    _require(len({tuple(row) for row in miller}) == count, 'oracle found duplicate Fourier indices')
    return reciprocal, miller, fields


def _realspace(miller: np.ndarray, fields: np.ndarray, fft: tuple[int, int, int]) -> np.ndarray:
    indices = np.mod(miller, np.array(fft))
    flat = np.ravel_multi_index(indices.T, fft)
    _require(len(np.unique(flat)) == len(miller), 'Fourier labels alias on declared FFT grid')
    spectrum = np.zeros((*fft, 2), dtype=np.complex128)
    spectrum[indices[:, 0], indices[:, 1], indices[:, 2], :] = fields
    # QE rho(G) is a density, so inverse synthesis is sum_G rho(G)e^(iGr),
    # without NumPy's default 1/N. norm='forward' provides that convention.
    return np.fft.ifftn(spectrum, axes=(0, 1, 2), norm='forward')


def _binding(path: Path) -> dict:
    data = path.read_bytes()
    return {'path': str(path.resolve()), 'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def run_oracle(source_save: Path, target_file: Path) -> dict:
    """Independent file/IFFT consumer, with frozen acceptance and bad-map probes."""
    source_path = source_save / 'charge-density.dat'
    xml_path = source_save / 'data-file-schema.xml'
    root = ET.parse(xml_path).getroot()
    for element in root.iter():
        element.tag = element.tag.rsplit('}', 1)[-1]
    _require(root.attrib.get('Units') == 'Hartree atomic units', 'unexpected XML units')
    output = root.find('output')
    _require(output is not None, 'source XML missing output')
    cell = np.column_stack([np.fromstring(output.findtext(f'atomic_structure/cell/a{i}', ''), sep=' ')
                            for i in (1, 2, 3)])
    source_volume = abs(float(np.linalg.det(cell)))
    fft_node = output.find('basis_set/fft_grid')
    _require(fft_node is not None, 'source XML missing actual FFT grid')
    source_fft = tuple(int(fft_node.attrib[f'nr{i}']) for i in (1, 2, 3))
    target_fft = (2 * source_fft[0], source_fft[1], source_fft[2])
    declared = np.array([float(output.findtext('band_structure/nelec', 'nan')),
                         float(output.findtext('magnetization/total', 'nan'))])
    source_b, source_m, source_f = _consume(source_path)
    target_b, target_m, target_f = _consume(target_file)
    _require(source_fft == (24, 24, 24) and target_fft == (48, 24, 24),
             'actual FFT grids differ from frozen preregistration')
    _require(len(source_m) == 4645 and len(target_m) == 9271,
             'actual G counts differ from frozen preregistration')
    expected_b = source_b.copy()
    expected_b[:, 0] *= 0.5
    _require(np.array_equal(target_b, expected_b), 'target physical reciprocal header is not exact replication')
    source_reciprocal_volume = (2 * np.pi)**3 / abs(float(np.linalg.det(source_b)))
    target_volume = (2 * np.pi)**3 / abs(float(np.linalg.det(target_b)))
    _require(abs(source_reciprocal_volume - source_volume) <= RELATIVE_ALLOWANCE * source_volume,
             'source reciprocal determinant disagrees with XML cell volume')
    _require(abs(target_volume - 2 * source_volume) <= RELATIVE_ALLOWANCE * 2 * source_volume,
             'target reciprocal header does not double real-space volume')
    source_grid = _realspace(source_m, source_f, source_fft)
    target_grid = _realspace(target_m, target_f, target_fft)
    expected_grid = np.concatenate((source_grid, source_grid), axis=0)
    norms = np.max(np.abs(source_grid), axis=(0, 1, 2))
    limits = RELATIVE_ALLOWANCE * np.maximum(1.0, norms)

    def consume_errors(fields: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        grid = _realspace(target_m, fields, target_fft)
        errors = np.max(np.abs(grid - expected_grid), axis=(0, 1, 2))
        integrals = np.mean(grid, axis=(0, 1, 2)) * target_volume
        return errors, integrals

    errors = np.max(np.abs(target_grid - expected_grid), axis=(0, 1, 2))
    _require(np.all(errors <= limits), 'actual target IFFT does not spatially replicate both source fields')
    source_imaginary = np.max(np.abs(source_grid.imag), axis=(0, 1, 2))
    target_imaginary = np.max(np.abs(target_grid.imag), axis=(0, 1, 2))
    _require(np.all(source_imaginary <= limits) and np.all(target_imaginary <= limits),
             'source or target IFFT is not real within frozen allowance')
    source_integrals = np.mean(source_grid, axis=(0, 1, 2)) * source_volume
    target_integrals = np.mean(target_grid, axis=(0, 1, 2)) * target_volume
    integral_limits = RELATIVE_ALLOWANCE * np.maximum(1.0, np.abs(declared))
    _require(np.all(np.abs(source_integrals - declared) <= integral_limits),
             'source grid charge/magnetization is not normalized to completed native XML')
    _require(np.all(np.abs(target_integrals - 2 * source_integrals) <= 2 * integral_limits),
             'target grid charge/magnetization does not double with volume')
    # Both negative controls use precisely the same IFFT and integral consumer.
    factor_errors, factor_integrals = consume_errors(2.0 * target_f)
    _require(np.all(factor_errors > limits), 'oracle fails to reject factor-two error in each field')
    _require(np.all(np.abs(factor_integrals - 2 * source_integrals) > 2 * integral_limits),
             'charge/magnetization consumer fails to reject doubled coefficients')
    wrong_odd = target_f.copy()
    source_by_miller = {tuple(row): source_f[i] for i, row in enumerate(source_m)}
    changed = 0
    for index, (h, k, l) in enumerate(target_m):
        if h % 2:
            value = source_by_miller.get(((h - 1) // 2, k, l))
            if value is not None:
                wrong_odd[index] = value
                changed += 1
    _require(changed > 0, 'odd-sector sensitivity probe has no populated Fourier indices')
    odd_errors, odd_integrals = consume_errors(wrong_odd)
    _require(np.all(odd_errors > limits), 'oracle fails to reject wrong odd-sector copying in each field')
    return {
        'schema_version': 1, 'relative_allowance': RELATIVE_ALLOWANCE,
        'source': {'density': _binding(source_path), 'xml': _binding(xml_path),
                   'fft': list(source_fft), 'ngm': len(source_m), 'volume_bohr3': source_volume},
        'target': {'density': _binding(target_file), 'fft': list(target_fft),
                   'ngm': len(target_m), 'volume_bohr3': target_volume},
        'fields': ['total_density', 'magnetization'],
        'spatial_max_abs_errors': errors.tolist(), 'spatial_abs_limits': limits.tolist(),
        'source_field_max_abs': norms.tolist(),
        'source_imaginary_max_abs': source_imaginary.tolist(),
        'target_imaginary_max_abs': target_imaginary.tolist(),
        'source_integrals': source_integrals.real.tolist(),
        'target_integrals': target_integrals.real.tolist(),
        'declared_source_integrals': declared.tolist(),
        'factor_two_negative_control': {'rejected_both_fields': True,
                                        'spatial_max_abs_errors': factor_errors.tolist(),
                                        'integrals': factor_integrals.real.tolist()},
        'odd_copy_negative_control': {'rejected_both_fields': True, 'copied_odd_count': changed,
                                      'spatial_max_abs_errors': odd_errors.tolist(),
                                      'integrals': odd_integrals.real.tolist()},
        'passed': True, 'native_completion_claimed': False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    inventory_path = Path(__file__).with_name('source-inventory.json')
    inventory = json.loads(inventory_path.read_text())
    source_default = next(Path(item['path']).parent for item in inventory['files']
                          if Path(item['path']).name == 'charge-density.dat')
    parser.add_argument('--source-save', type=Path, default=source_default)
    parser.add_argument('--ordering-executable', type=Path, required=True)
    parser.add_argument('--target-file', type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='qe-density-oracle-') as temporary:
        target = args.target_file or Path(temporary) / 'charge-density.dat'
        receipt = fold_density(args.source_save, target, args.ordering_executable)
        evidence = run_oracle(args.source_save, target)
        print(json.dumps({'transformation': receipt, 'oracle': evidence}, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
