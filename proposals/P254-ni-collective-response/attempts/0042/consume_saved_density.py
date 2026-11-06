"""Read paid density spectral rows and consume original saved CAR outputs; no replay."""
from pathlib import Path
import argparse
import json
import resource
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[3] / 'src')] + [
    str(HERE.parent / name) for name in ('0038', '0037', '0031', '0030', '0028')]
import numpy as np
from export_frozen_bath_comparator import car_moments, projected_factors
from expanded_influence_readout import save_json
from integrate_saved_spectral import load_saved_spectral
from substrate_framework.streamed_material_state import authenticate, file_binding


def bound_json(path, size, digest):
    binding = dict(path=str(path), bytes=size, sha256=digest)
    return json.loads(Path(authenticate(binding)).read_bytes()), binding


def consume(result, manifest):
    spec = manifest['specification']
    for field, key in (('basis', 'basis'), ('current_bank', 'current_bank'),
                       ('spectral_archive', 'spectral_reference'), ('filter_archive', 'filter')):
        if spec[field] != result[key]:
            raise ValueError('historical row authority differs: ' + field)
    for binding in spec['implementation']:
        authenticate(binding)
    saved = load_saved_spectral(result['spectral_reference'])
    v, vi = saved['right_eigenvectors'], saved['inverse_right_eigenvectors']
    size = len(v)
    forward = np.empty((3, size, size), complex)
    backward = np.empty_like(forward)
    cursor = 0
    for key, item in sorted(manifest['items'].items()):
        start, stop = item['metadata']['row_start'], item['metadata']['row_stop']
        if start != cursor or not start < stop <= size:
            raise ValueError('missing or overlapping paid spectral rows: ' + key)
        with np.load(authenticate(item['archive']), allow_pickle=False) as chunk:
            f, b = chunk['forward'], chunk['backward']
        if f.shape != (3, stop-start, size) or b.shape != (3, size, stop-start):
            raise ValueError('actual three-probe spectral row shapes differ')
        forward[:, start:stop], backward[:, :, start:stop] = f, b
        cursor = stop
    if cursor != size or manifest['status'] != 'complete':
        raise ValueError('original spectral convolution is incomplete')
    reconstructed = v @ forward @ vi
    reconstructed_adjoint = vi.conj().T @ backward @ v.conj().T
    with np.load(authenticate(result['operators']), allow_pickle=False) as archive:
        matrices = dict(density=archive['density_residual'], matched=archive['matched'])
        adjoints = dict(density=archive['density_residual_adjoint'],
                        matched=archive['matched_adjoint'])
    defects = dict(density=float(np.linalg.norm(reconstructed[2]-matrices['density'])),
                   density_adjoint=float(np.linalg.norm(reconstructed_adjoint[2]-adjoints['density'])))
    exact_rows = (np.array_equal(reconstructed[2], matrices['density'])
                  and np.array_equal(reconstructed_adjoint[2], adjoints['density']))
    if not exact_rows:
        raise ValueError('paid row reconstruction is not byte-identical to original density output')
    basis_packet = json.loads(Path(authenticate(result['inputs']['a37_basis_complete'])).read_bytes())
    with np.load(authenticate(result['basis']), allow_pickle=False) as archive:
        frame, gram, density = (archive[key] for key in
                                ('projected_frame', 'gram', 'density_coefficients'))
    factors = projected_factors(matrices, adjoints, frame, basis_packet['basis_identity'],
                                result['original_state_authority'], 'bohr')
    moments = car_moments(frame, gram, density, factors, result['original_factor_identity'],
                          basis_packet['basis_identity'], result['original_state_authority'])
    def decode(value):
        return complex(value['real'], value['imag'])
    mean_errors = {name: abs(value-decode(result['moments']['means'][name]))
                   for name, value in moments['means'].items()}
    original = {(row['left'], row['right']): decode(row['value'])
                for row in result['moments']['covariance']}
    errors = {row['left']+'->'+row['right']:
              abs(row['value']-original[row['left'], row['right']])
              for row in moments['covariance']}
    if any(mean_errors.values()) or any(errors.values()):
        raise ValueError('saved actual density/matched CAR outputs differ from producer')
    return dict(status='complete', paid_row_blocks=len(manifest['items']),
                all_original_rows=size, probes=3, exact_density_reconstruction=exact_rows,
                matrix_defects=defects, mean_errors=mean_errors, covariance_errors=errors,
                actual_saved_moments=moments,
                scope='original source0/Gamma density and matched saved-output consumer only',
                work=dict(mean_solves=0, filter_integrations=0, physical_vertex_evaluations=0,
                          physical_H_actions=0, spectral_eigensolves=0, moment_integrations=0))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('result', 'rows', 'output'):
        parser.add_argument('--'+name, required=True, type=Path)
    for name in ('result', 'rows'):
        parser.add_argument('--'+name+'-bytes', required=True, type=int)
        parser.add_argument('--'+name+'-sha256', required=True)
    args = parser.parse_args()
    started = perf_counter()
    result, rb = bound_json(args.result, args.result_bytes, args.result_sha256)
    manifest, mb = bound_json(args.rows, args.rows_bytes, args.rows_sha256)
    output = consume(result, manifest)
    output.update(original_result=rb, original_rows=mb, implementation=file_binding(__file__),
                  elapsed_seconds=perf_counter()-started,
                  maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    binding = save_json(args.output, output)
    print(json.dumps(dict(result=binding, status=output['status'],
                          exact_density_reconstruction=output['exact_density_reconstruction'],
                          covariance_errors=output['covariance_errors'],
                          elapsed_seconds=output['elapsed_seconds']), allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
