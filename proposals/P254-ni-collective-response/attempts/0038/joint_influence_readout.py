"""Assemble the ionic/current operator sum and its complete ordered CAR readout."""
from __future__ import annotations

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[3] / 'src')] + [str(HERE.parent / name)
    for name in ('0037', '0031', '0030', '0029', '0028')]

import numpy as np
from export_frozen_bath_comparator import car_moments, projected_factors
from integrate_saved_spectral import integrate_saved_spectral

LABELS = ('ionic', 'field', 'joint', 'retained', 'omitted')
PARTS = ('full_vertex', 'retained_vertex', 'omitted_vertex')
SCOPE = ('Original source0/Gamma/Ni1-x t=30 inverse Ry; unchanged nonstationary '
    'S/D/R, all346 original columns plus692 A37 extra directions, frozen raw H. '
    'Joint two-mode ionic-force and uniform-current electronic contribution with '
    'both ordered cross terms in the same finite bath. Conditional forced output, '
    'NOT independently additive noise or corrected A32 total; NOT full moving, '
    'all128-source/64q, density/backreaction/loss closure, convergence, error bound, '
    'physical damping, lifetime, nuclear rate or whole physical E/C/N completion.')


def field_source_amplitude(ncells, weight, matrix):
    """A36 source0 covariance N*w*M*C_J*M.T implies amplitude sqrt(N*w)*M."""
    if type(ncells) is not int or ncells < 1 or not np.isfinite(weight) or weight <= 0:
        raise ValueError('positive exact cell count and source weight required')
    matrix = np.asarray(matrix)
    if np.iscomplexobj(matrix) or matrix.shape != (3, 3) or not np.all(np.isfinite(matrix)):
        raise ValueError('real finite original three-polarization map required')
    return np.sqrt(ncells * weight) * matrix


def integrate_field(saved, currents, current_adjoints, projector, bounds, field, amplitude):
    """Route raw J and Pactive J Pactive through six channels and two probes.

    The field filter has canonical Pi-force coefficients, not Pi state noise.
    Ionic F/sqrt(N) operators are neither consumed nor renormalized here.
    """
    field = np.asarray(field, dtype=np.complex128)
    if field.ndim != 4 or field.shape[1:] != (8, 3, 1):
        raise ValueError('require original field[segment,degree7,Pi-force,first probe]')
    current = np.asarray(currents, dtype=np.complex128)
    adjoint = np.asarray(current_adjoints, dtype=np.complex128)
    if current.shape[0] != 3 or adjoint.shape != current.shape:
        raise ValueError('three Cartesian current and true-adjoint channels required')
    retained = np.stack([projector @ value @ projector for value in current])
    retained_adjoint = np.stack([projector @ value @ projector for value in adjoint])
    weights = np.einsum('sdap,ac->sdcp', field, amplitude)
    batched = np.zeros(weights.shape[:2] + (6, 2), dtype=np.complex128)
    batched[:, :, :3, 0] = weights[:, :, :, 0]
    batched[:, :, 3:, 1] = weights[:, :, :, 0]
    return integrate_saved_spectral(saved, np.concatenate((current, retained)),
        np.concatenate((adjoint, retained_adjoint)), bounds, batched)


def joint_operators(ionic, ionic_adjoints, filtered):
    """Keep paid ionic legs untouched; add new field legs before contraction."""
    field = dict(zip(PARTS, (filtered['operators'][0], filtered['operators'][1],
        filtered['operators'][0] - filtered['operators'][1])))
    field_adjoints = dict(zip(PARTS, (filtered['adjoint_operators'][0],
        filtered['adjoint_operators'][1],
        filtered['adjoint_operators'][0] - filtered['adjoint_operators'][1])))
    joint = {name: ionic[name] + field[name] for name in PARTS[:2]}
    joint['omitted_vertex'] = joint['full_vertex'] - joint['retained_vertex']
    joint_adjoints = {name: ionic_adjoints[name] + field_adjoints[name] for name in PARTS[:2]}
    joint_adjoints['omitted_vertex'] = (joint_adjoints['full_vertex']
                                     - joint_adjoints['retained_vertex'])
    matrices = dict(zip(LABELS, (ionic[PARTS[0]], field[PARTS[0]],
        joint[PARTS[0]], joint[PARTS[1]], joint[PARTS[2]])))
    adjoints = dict(zip(LABELS, (ionic_adjoints[PARTS[0]], field_adjoints[PARTS[0]],
        joint_adjoints[PARTS[0]], joint_adjoints[PARTS[1]], joint_adjoints[PARTS[2]])))
    arrays = {'field_operator_' + name: value for name, value in field.items()}
    arrays.update({'field_adjoint_operator_' + name: value for name, value in field_adjoints.items()})
    arrays.update({'joint_operator_' + name: value for name, value in joint.items()})
    arrays.update({'joint_adjoint_operator_' + name: value for name, value in joint_adjoints.items()})
    return matrices, adjoints, arrays


def joint_moments(bank, matrices, adjoints, basis_identity, factor_identity, authority):
    """Reuse unchanged A31 CAR, including its small state-admissibility diagnostics."""
    frame = bank['projected_frame']
    return car_moments(frame, bank['gram'], bank['density_coefficients'],
        projected_factors(matrices, adjoints, frame, basis_identity, authority, 'bohr'),
        factor_identity, basis_identity, authority)


def joint_reconstruction(moments, matrices, adjoints):
    values = {(row['left'], row['right']): row['value'] for row in moments['covariance']}
    if len(moments['covariance']) != 25 or set(values) != {(a, b) for a in LABELS for b in LABELS}:
        raise ValueError('all25 ordered complex ionic/field/joint/retained/omitted pairs required')
    return dict(ionic_field_cross_order=values['ionic', 'field'],
        field_ionic_cross_order=values['field', 'ionic'],
        raw_joint_minus_ionic_field_both_cross_orders=values['joint', 'joint']
            - values['ionic', 'ionic'] - values['field', 'field']
            - values['ionic', 'field'] - values['field', 'ionic'],
        raw_joint_minus_retained_omitted_sum=values['joint', 'joint']
            - sum(values[a, b] for a in ('retained', 'omitted') for b in ('retained', 'omitted')),
        raw_joint_mean_reconstruction=moments['means']['joint']
            - moments['means']['ionic'] - moments['means']['field'],
        raw_joint_operator_reconstruction_fro=float(np.linalg.norm(
            matrices['joint'] - matrices['ionic'] - matrices['field'])),
        raw_joint_adjoint_reconstruction_fro=float(np.linalg.norm(
            adjoints['joint'] - adjoints['ionic'] - adjoints['field'])),
        raw_joint_retained_omitted_operator_fro=float(np.linalg.norm(
            matrices['joint'] - matrices['retained'] - matrices['omitted'])),
        raw_joint_retained_omitted_adjoint_fro=float(np.linalg.norm(
            adjoints['joint'] - adjoints['retained'] - adjoints['omitted'])),
        raw_supplied_adjoint_defects_fro={name: float(np.linalg.norm(
            adjoints[name] - value.conj().T)) for name, value in matrices.items()},
        raw_ordered_covariance_conjugacy_defects=[dict(left=a, right=b,
            value=values[a, b] - values[b, a].conjugate()) for a in LABELS for b in LABELS])
