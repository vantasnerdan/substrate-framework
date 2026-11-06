#!/usr/bin/env python3
"""Independent finite scalar/preparation/identity checks, NOT material results.

No material H0, coefficients, trajectories or original frame storage are run.
The separate production initialize/read-back commands prove the actual data path.
"""
import argparse
from copy import deepcopy
import json
import os
from pathlib import Path
from types import SimpleNamespace
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
import matched_reference_io as io
import numpy as np
from matched_joint_case import joint_blocks, model_identity
from matched_joint_sources import registration, validate_aggregate
from substrate_framework.material_matched_reference import slot_identity
from substrate_framework.material_joint_field_scalar import coordinate_correction, JointFieldCoefficients
from substrate_framework.material_modes import AMU_RY


TOL = 1e-10


def residual(a, b):
    return float(np.max(np.abs(np.asarray(a)-b)))/(1+float(np.max(np.abs(b))))


def scalar_checks(coeff, masses, charges):
    # Entrywise polynomial independent of coordinate_correction's matrix code.
    u = np.array([.2, -.3, .13, -.1, .07, .19])
    p = np.array([.5, -.7, .11, .31, -.23, .09])
    a = np.array([.13, -.21, .09])
    ell = np.array([-.12, .34, -.17])
    common = np.diag([2., 3., 5., 7., 11., 13.])
    common[0, 4] = common[4, 0] = .17
    # Physical aggregate C already includes every spatial onsite correction.
    energy = sum(common[i,j]*u[i]*u[j]/2 for i in range(6) for j in range(6))
    energy += sum((p[i]-charges[i//3]*a[i%3])**2/(2*masses[i]) for i in range(6))
    energy += sum(a[i]*ell[i] for i in range(3))
    energy += sum(a[i]*coeff.contact[i,j]*a[j]/2 for i in range(3) for j in range(3))
    energy += sum(u[i]*coeff.mixed[i,j]*a[j] for i in range(6) for j in range(3))
    force = np.array([-sum(common[i,j]*u[j] for j in range(6))
                      -sum(coeff.mixed[i,j]*a[j] for j in range(3)) for i in range(6)])
    gradient = np.array([ell[j]+sum(coeff.contact[j,k]*a[k] for k in range(3))
        +sum(coeff.mixed[i,j]*u[i] for i in range(6))
        -sum(charges[n]*(p[3*n+j]-charges[n]*a[j])/masses[3*n+j] for n in range(2)) for j in range(3)])
    result = coordinate_correction(coeff, u, p, a, masses, charges, np.zeros(6), ell)
    actual_energy = u@common@u/2+sum(p*p/masses)/2+result['energy']
    return dict(scalar=residual(actual_energy, energy),
        force=residual(-common@u+result['force_correction'], force),
        field_gradient=residual(result['field_gradient'], gradient),
        canonical_velocity=residual(result['velocity'], [(p[i]-charges[i//3]*a[i%3])/masses[i] for i in range(6)]))


def preparation_checks(field, ionic, photon, masses, charges, contact, diagonal):
    d = len(masses)
    # Independently construct oscillator widths and ordered CCR blocks.
    w = np.zeros((2*d, 2*d), complex)
    for i in range(d):
        width = (masses[i]*diagonal[i])**.5
        w[i,i], w[d+i,d+i] = 1/(2*width), width/2
        w[i,d+i], w[d+i,i] = .5j, -.5j
    stiffness = np.array([[64*(contact[i,j]+(sum(charges[n]**2/masses[3*n+i] for n in range(2)) if i == j else 0))
                           for j in range(3)] for i in range(3)])
    q, p = photon['covariance_qq'], photon['covariance_pp']
    symplectic = np.block([[np.zeros((3,3)), np.eye(3)], [-np.eye(3), np.zeros((3,3))]])
    return dict(ionic_ordered=residual(ionic.local.ordered, w),
        ionic_mean=residual(ionic.mean_ordered, w/64),
        q_minus_pair=residual(ionic.q_pair((1,2,3),(3,2,1)), w),
        wrong_q_pair=residual(ionic.q_pair((1,0,0),(1,0,0)), np.zeros_like(w)),
        field_stiffness=residual(photon['clamped_stiffness_raw'], stiffness),
        field_stationary=residual(stiffness@q, p/field.mass),
        field_purity=residual(q@p, np.eye(3)/4),
        field_ccr=residual(photon['ordered_covariance']-photon['ordered_covariance'].T, 1j*symplectic),
        maxwell_mass=residual(field.mass, 64*24/(8*np.pi)))


def identity_checks():
    rows = []
    for spin in (1,2):
        for index in np.ndindex(4,4,4):
            s = SimpleNamespace(ispin=spin, mesh_index=index, source_ik=1,
                source_binding=dict(path='/finite/source', bytes='17', sha256='finite-source-'+str(spin)),
                weight=1/64, operation=SimpleNamespace(index=1), orbital_conjugation=False,
                reciprocal_shift=np.zeros(3, int))
            rows.append((len(rows),s,np.array([.3,.7])))
    context = dict(case='finite-case-A', population_ledger=dict(sha256='finite-ledger-A'))
    lineage = [dict(slot=slot_identity(s), receipt=dict(sha256='finite-slot-'+str(i))) for i,s,_ in rows]
    value = dict(status='complete', complete_case=True, context=context, slot_count=128,
        model='matched-initial-source-orbit-declared-geometric-v1',
        positive_time_symmetry=False, source_coefficients_recomputed=False,
        coefficient_defect_upper_bound=None, coefficient_lineage=lineage,
        source_receipts=[x['receipt'] for x in lineage], original_source_count=2)
    validate_aggregate(value, context, rows)
    mutations = {}
    for name in ('wrong_case_same_shape', 'missing_slot', 'duplicate_slot', 'wrong_source', 'wrong_weight'):
        wrong = deepcopy(value)
        if name == 'wrong_case_same_shape': wrong['context']['case'] = 'finite-case-B'
        elif name == 'missing_slot': wrong['coefficient_lineage'].pop()
        elif name == 'duplicate_slot': wrong['coefficient_lineage'][1] = wrong['coefficient_lineage'][0]
        elif name == 'wrong_source': wrong['coefficient_lineage'][0]['slot']['source']['sha256'] = 'other'
        else: wrong['coefficient_lineage'][0]['slot']['weight'] *= 2
        try:
            validate_aggregate(wrong, context, rows)
        except ValueError:
            mutations[name] = True
        else:
            mutations[name] = False
    changed = deepcopy(context)
    changed['population_ledger']['sha256'] = 'finite-ledger-B'
    mutations['model_binds_original_ledger'] = model_identity(changed) != model_identity(context)
    return mutations


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    reg, _, sources = registration()
    diagonal = np.array([2.,3.,5.,7.,11.,13.])
    lattice = np.zeros((4,4,4,2,3,2,3), complex)
    lattice[0,0,0] = np.diag(diagonal).reshape(2,3,2,3)
    arrays = dict(atomic_masses_amu=np.array([1.7,5.3]), valence_charges=np.array([2.,7.]),
        cell_bohr=np.diag([2.,3.,4.]), common_force_lattice_ry_bohr2=lattice,
        contact=np.array([[2.,.11,-.07],[.11,3.,.13],[-.07,.13,4.]]),
        mixed_ionic_field=np.arange(18).reshape(2,3,3)/73-.09)
    coeff, field, ionic, photon, masses, charges = joint_blocks(arrays, authority='explicit finite oracle fixture')
    errors = scalar_checks(coeff,masses,charges)
    errors.update(preparation_checks(field,ionic,photon,masses,charges,arrays['contact'],diagonal))
    errors['original_masses'] = residual(masses, np.repeat(arrays['atomic_masses_amu']*AMU_RY,3))
    errors['no_second_delta'] = residual(coeff.delta_curvature, np.zeros((6,6)))
    mutations = identity_checks()
    doubled = JointFieldCoefficients(np.diag(diagonal), coeff.contact, coeff.mixed, 'intentional duplicate curvature')
    mutations['duplicate_spatial_correction'] = max(scalar_checks(doubled,masses,charges).values()) > TOL
    wrong = dict(photon, ordered_covariance=photon['ordered_covariance'].T)
    mutations['field_ordered_sign'] = preparation_checks(field,ionic,wrong,masses,charges,arrays['contact'],diagonal)['field_ccr'] > TOL
    wrong = dict(photon, clamped_stiffness_raw=photon['clamped_stiffness_raw']/64)
    mutations['global_field_normalization'] = preparation_checks(field,ionic,wrong,masses,charges,arrays['contact'],diagonal)['field_stiffness'] > TOL
    status = 'passed' if max(errors.values()) <= TOL and all(mutations.values()) else 'failed'
    result = dict(schema='P254-matched-joint-oracle-v1', status=status, registration=reg,
        implementation_sources=sources, scope='independent finite joins only; not material evidence',
        errors=errors, tolerance=TOL, mutations=mutations,
        actual_case_proof='run_matched_joint_case.py initialize followed by read-back', material_trajectory=False)
    io.save_json(args.output,result)
    if status != 'passed': raise SystemExit(2)


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()
