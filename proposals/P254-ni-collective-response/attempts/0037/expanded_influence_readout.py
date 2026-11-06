"""Persist complex A37 readouts and label the old8/new692 model contrast."""
from __future__ import annotations

import json
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[3] / 'src'), str(HERE.parent / '0028')]

import numpy as np
from clamped_scientific_readout import encode
from substrate_framework.streamed_material_state import file_binding


VERTICES = ('full_vertex', 'retained_vertex', 'omitted_vertex')
SCOPE = (
    'Original source0/Gamma/Ni1-x final output at t=30 inverse Ry; two selected '
    'ionic force modes; all346 original raw state columns and692 extra empty '
    'directions; frozen literal initial H; NONSTATIONARY original state. '
    'Initial Pactive force-vertex split propagated within the same finite bath; '
    'all ordered cross terms retained. Conditional forced-output electronic '
    'contribution, NOT independent additive noise or a corrected A32 total. '
    'Not full moving/backreaction/RPA closure, all128 sources/64q, field/collision '
    'loss, convergence, an error bound, physical damping, lifetime, heat, or '
    'whole physical E/C/N completion.'
)


def _native(value):
    """Remove NumPy scalar leaves BEFORE the existing complex-pair encoder."""
    if isinstance(value, dict):
        return {key: _native(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_native(item) for item in value]
    if isinstance(value, np.ndarray):
        return _native(value.tolist())
    if isinstance(value, np.generic):
        return _native(value.item())
    return value


def _sync_directory(directory):
    descriptor = os.open(directory, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def save_json(path, value):
    """Publish only complete, fsynced strict JSON, never replace an old result."""
    path = Path(path)
    payload = json.dumps(encode(_native(value)), indent=2, sort_keys=True,
                         allow_nan=False) + '\n'
    with path.open('x') as stream:
        stream.write(payload)
        stream.flush()
        os.fsync(stream.fileno())
    _sync_directory(path.parent)
    return file_binding(path)


def save_npz(path, arrays):
    path = Path(path)
    with path.open('xb') as stream:
        np.savez(stream, **arrays)
        stream.flush()
        os.fsync(stream.fileno())
    _sync_directory(path.parent)
    _sync_directory(path.parent.parent)
    return file_binding(path)


def covariance_map(moments, *, encoded=False):
    rows = moments['covariance']
    values = {}
    for row in rows:
        value = row['value']
        values[row['left'], row['right']] = (
            complex(value['real'], value['imag']) if encoded else value)
    if len(rows) != 9 or set(values) != {(a, b) for a in VERTICES for b in VERTICES}:
        raise ValueError('require all9 ordered full/retained/omitted covariances')
    return values


def reconstruction(moments, matrices, adjoints):
    covariance = covariance_map(moments)
    parts = VERTICES[1:]
    return dict(
        raw_full_minus_retained_omitted_sum=covariance['full_vertex', 'full_vertex']
            - sum(covariance[a, b] for a in parts for b in parts),
        raw_mean_reconstruction= moments['means']['full_vertex']
            - sum(moments['means'][name] for name in parts),
        raw_operator_reconstruction_fro=float(np.linalg.norm(
            matrices['full_vertex'] - matrices['retained_vertex'] - matrices['omitted_vertex'])),
        raw_adjoint_operator_reconstruction_fro=float(np.linalg.norm(
            adjoints['full_vertex'] - adjoints['retained_vertex'] - adjoints['omitted_vertex'])),
        raw_supplied_adjoint_defect_fro={name: float(np.linalg.norm(
            adjoints[name] - matrix.conj().T)) for name, matrix in matrices.items()},
        raw_ordered_covariance_conjugacy_defect=[dict(left=a, right=b,
            value=covariance[a, b] - covariance[b, a].conjugate())
            for a in VERTICES for b in VERTICES])


def model_contrast(reference, moments, reference_binding):
    old = covariance_map(reference['filtered_covariance'], encoded=True)
    new = covariance_map(moments)
    old_means = reference['filtered_covariance']['means']
    return dict(
        reference_result=reference_binding,
        interpretation='Same observable/state/source/normalization; old8 versus new692 '
            'extra-direction model-enrichment contrast only. Differences are NOT '
            'error bounds, convergence, independent noise or corrected A32 totals.',
        extra_rank_old=8, extra_rank_new=692,
        ordered_covariance=[dict(left=a, right=b, old8=old[a, b], new692=new[a, b],
            new692_minus_old8=new[a, b] - old[a, b]) for a in VERTICES for b in VERTICES],
        means=[dict(vertex=name, old8=complex(old_means[name]['real'], old_means[name]['imag']),
            new692=moments['means'][name], new692_minus_old8=moments['means'][name]
            - complex(old_means[name]['real'], old_means[name]['imag'])) for name in VERTICES])


def raw_bank_defects(bank):
    h, ha = bank['hamiltonian'], bank['hamiltonian_adjoint']
    rho, projector = bank['projected_density'], bank['active_projector']
    frame = bank['projected_frame']
    return dict(
        raw_state_commutator_frobenius=float(np.linalg.norm(h @ rho - rho @ h)),
        raw_hamiltonian_supplied_adjoint_defect_fro=float(np.linalg.norm(ha - h.conj().T)),
        raw_projected_density_hermiticity_fro=float(np.linalg.norm(rho - rho.conj().T)),
        raw_projected_density_coordinate_defect_fro=float(np.linalg.norm(
            rho - frame @ bank['density_coefficients'] @ frame.conj().T)),
        raw_frame_gram_defect_fro=float(np.linalg.norm(frame.conj().T @ frame - bank['gram'])),
        raw_active_projector_hermiticity_fro=float(np.linalg.norm(projector - projector.conj().T)),
        raw_active_projector_idempotence_fro=float(np.linalg.norm(projector @ projector - projector)),
        raw_defect_policy='Raw defects only; no clipping, symmetrization, repair or rank/quality gate.')
