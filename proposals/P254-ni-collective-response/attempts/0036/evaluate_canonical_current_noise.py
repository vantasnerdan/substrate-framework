#!/usr/bin/env python3
"""Admit actual A35 saved sectors and map their raw initial canonical force noise.

Only JSON/source authentication, archive cell extraction and small matrix maps
are performed. No preparation, vertex, Hamiltonian or trajectory work is run.
Exit 0 means complete128, 2 valid unrenormalized partial, and 1 invalid.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(HERE), str(HERE.parent / '0035'), str(ROOT / 'src')]
from aggregate_uniform_current_noise import (
    PROBES, PARTS, aggregate, authenticate, binding, complex_record, complex_sum,
    complex_value, load_case, load_context, load_json, parse_covariance,
    same_binding, save_json,
)
from canonical_current_noise import canonical_force_covariance, current_force_matrix
from substrate_framework.physical_field_environment import UniformMaxwellField
from substrate_framework.prepared_reference import array

PLAN = HERE / 'implementation-contract-001.json'
PLAN_BYTES = 7113
PLAN_SHA256 = 'f6f8a3e9b182db6fcf5a502fdf4d762f71d2c382efd98faf61a89c12136fafdd'
SCHEMA = 'P254-A36-initial-canonical-uniform-field-force-noise-v1'
SCOPE = (
    'Proposed initial external-Gamma electronic generalized-force covariance entering '
    'the ORIGINAL global canonical Pi-dot. Preserves raw complex ordered pairs and '
    'RR/RE/ER/EE, original state and full/partial source coverage. Not initial Pi-state '
    'covariance, independent added noise, total force-current joint bank, ionic noise, '
    'moving influence kernel, finite-q density, stationarity/FDT, temperature, loss, '
    'heat, linewidth, lifetime, rates or whole E/C/N closure. Independent map oracle '
    'is separate; this consumer does not assert that oracle has passed.'
)


def fixed_plan():
    receipt = binding(PLAN)
    if receipt['bytes'] != PLAN_BYTES or receipt['sha256'] != PLAN_SHA256:
        raise ValueError('fixed settled adjacent A36 plan bytes required')
    plan = load_json(PLAN)
    if plan['schema'] != 'P254-A36-canonical-uniform-field-current-noise-plan-v1':
        raise ValueError('wrong canonical current-noise plan')
    sources = [authenticate(item) for item in plan['settled_sources']]
    originals = {key: authenticate(plan[key]) for key in (
        'original_initial', 'original_boundary', 'original_preparation',
        'matched_aggregate', 'matched_archive')}
    return receipt, plan, sources, originals


def admit_current(path, plan, originals):
    """Rebuild through the owning consumer, not a second summation convention."""
    receipt = binding(path)
    current = load_json(path)
    contract = authenticate(current['contract'], HERE.parent / '0035' / PLAN.name)
    pinned = next(item for item in plan['settled_sources']
                  if Path(item['path']).resolve() == Path(contract['path']))
    if not same_binding(contract, pinned):
        raise ValueError('current aggregate A35 contract differs from settled pin')
    inputs = current['inputs']
    if not same_binding(inputs['initial_result'], originals['original_initial']):
        raise ValueError('current aggregate initial identity differs from original')
    # The existing owner authenticates coefficients, all128 source/weight/selection
    # identities, the A34 scientific profile, and original boundary/preparation.
    context = load_context(inputs['initial_result']['path'],
                           inputs['coefficient_manifest']['path'], contract['path'])
    for key, expected in (('original_boundary', originals['original_boundary']),
                          ('preparation', originals['original_preparation'])):
        if not same_binding(context[key], expected):
            raise ValueError('current context differs from original: ' + key)
    if not isinstance(current['sectors'], list) or not current['sectors']:
        raise ValueError('at least one explicitly committed original sector required')
    cases = []
    for sector in current['sectors']:
        record = sector['original_record']
        if type(record) is not int or record not in range(128):
            raise ValueError('literal original sector record in 0..127 required')
        done = authenticate(sector['done'])
        cases.append(load_case(done['path'], context, expected_record=record))
    rebuilt = aggregate(cases, context)
    # Compare every key, ordered list, numeric representation and covariance, not
    # merely metadata labels. Object key order is immaterial; no tolerances apply.
    encode = lambda value: json.dumps(value, sort_keys=True, separators=(',', ':'),
                                      allow_nan=False)
    if encode(current) != encode(rebuilt):
        raise ValueError('saved current aggregate differs from actual A35 reconstruction')
    authenticate(receipt)
    return receipt, rebuilt, context, cases


def canonical_field(plan, originals, context):
    saved = context['initial']['saved']
    if not same_binding(saved['specification']['aggregate'], originals['matched_aggregate']):
        raise ValueError('original scalar matched aggregate differs from fixed plan')
    matched = load_json(originals['matched_aggregate']['path'])
    if matched['status'] != 'complete' or not same_binding(
            matched['archive'], originals['matched_archive']):
        raise ValueError('original matched archive linkage differs')
    construction = plan['canonical_construction']
    ncells, polarizations = construction['ncells'], construction['polarizations']
    if type(ncells) is not int or ncells != 64 or polarizations != np.eye(3).tolist():
        raise ValueError('original N64 and three identity polarizations required')
    # Extract only the existing physical cell; do not instantiate joint_blocks,
    # which would replay ionic and conditional-field preparation.
    with np.load(originals['matched_archive']['path'], allow_pickle=False) as archive:
        cell = array(archive['cell_bohr'], (3, 3), 'original physical cell', real=True)
    cell_volume = float(abs(np.linalg.det(cell)))
    field = UniformMaxwellField(ncells * cell_volume, polarizations,
        saved['model_identity'] + '; original BvK4^3, not cavity')
    return field, cell, ncells


def mapped_covariance(current, ncells, polarizations):
    full, blocks, defects = parse_covariance(current['covariance'])

    def matrix(values):
        return np.asarray([[values[(a, b)] for b in PROBES] for a in PROBES], complex)

    def mapped(values):
        return canonical_force_covariance(matrix(values), ncells, polarizations)

    force_full = mapped(full)
    force_blocks = {(pa, pb): mapped({(a, b): blocks[(a, b, pa, pb)]
                    for a in PROBES for b in PROBES}) for pa in PARTS for pb in PARTS}
    force_defect = mapped(defects)
    weighted_defects = {(row['left'], row['right']):
                        complex_value(row['weighted_raw_reconstruction_defect'])
                        for row in current['covariance']}
    force_weighted_defect = mapped(weighted_defects)
    current_reconstruction = {(row['left'], row['right']):
                              complex_value(row['block_reconstruction'])
                              for row in current['covariance']}
    force_reconstruction = mapped(current_reconstruction)
    names = ['canonical_field_force_' + axis for axis in 'xyz']
    rows = []
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            parts = [dict(part_left=pa, part_right=pb,
                          value=complex_record(force_blocks[pa, pb][i, j].item()))
                     for pa in PARTS for pb in PARTS]
            reconstructed = complex_sum(complex_value(part['value']) for part in parts)
            rows.append(dict(left=a, right=b, full=complex_record(force_full[i, j].item()),
                blocks=parts, block_reconstruction=complex_record(reconstructed),
                reconstruction_defect=complex_record((force_full[i, j] - reconstructed).item()),
                mapped_current_block_reconstruction=complex_record(force_reconstruction[i, j].item()),
                mapped_current_reconstruction_defect=complex_record(force_defect[i, j].item()),
                mapped_weighted_raw_reconstruction_defect=complex_record(
                    force_weighted_defect[i, j].item())))
    return names, rows


def evaluate(path):
    plan_receipt, plan, sources, originals = fixed_plan()
    current_receipt, current, context, cases = admit_current(path, plan, originals)
    field, cell, ncells = canonical_field(plan, originals, context)
    names, covariance = mapped_covariance(current, ncells, field.polarizations)
    saved = context['initial']['saved']
    return dict(schema=SCHEMA, status=current['status'], contract=plan_receipt,
        current_aggregate=current_receipt, current_contract=current['contract'],
        scientific_compatibility_identity=current['scientific_compatibility_identity'],
        settled_sources=sources, implementation_sources=[binding(HERE / name) for name in
            ('canonical_current_noise.py', 'evaluate_canonical_current_noise.py')],
        originals=originals, original_inputs=current['inputs'],
        original_model_identity=saved['original_model_identity'],
        joint_model_identity=saved['model_identity'],
        original_state=dict(case=context['initial']['case'], time_ry_inverse=0,
            conditional_preparation=context['initial']['conditional_preparation'],
            source_states=saved['specification']['source_states'],
            preparation=originals['original_preparation'],
            boundary=originals['original_boundary']),
        original_records=current['original_records'],
        missing_original_records=current['missing_original_records'],
        required_original_records=current['required_original_records'],
        sector_count=current['sector_count'],
        sectors=[dict(original_record=case['original_record'], done=case['done'],
                      result=case['result'], source=case['payload']['source'],
                      state_authority=case['payload']['state_authority']) for case in cases],
        probe_names=names, units='electronic_canonical_field_generalized_force_ry_bohr',
        covariance_units='(Ry*bohr)^2', field_normalization=field.normalization(),
        normalization=dict(ncells=ncells, physical_cell_bohr=cell.tolist(),
            primitive_volume_bohr3=float(abs(np.linalg.det(cell))),
            source_current_normalization=current['normalization'],
            current_force_matrix=current_force_matrix(field.polarizations).tolist(),
            rule='C_Pidot_el = N * M * [sum_explicit_spin_k(w_k*C_J,k)] * M.T; M=2e',
            mean_rule='Pi-dot_el=-N*e*<L_Ry>=2*N*e*<J_Ha>',
            fluctuation_rule='global source amplitude sqrt(N)*M; covariance factor N, not N^2',
            scalar_reason='E_total=N*E_matter_per_cell+Pi^2/(2*M_a); field Poisson bracket 1',
            conversion_reason='J_Ha=-L_Ry/2; -e*L_Ry=+2*e*J_Ha; each covariance leg gets 2',
            normalization_sources=['physical_field_environment.py:3-16,67-84',
                'material_joint_field_scalar.py:94-95,123-148',
                'matched_joint_case.py:38-39,118',
                'joint_bloch_response.py:224-240', 'joint_ordered_sources.py:121-124'],
            quadrature_applied_once=True, quadrature_weight=1 / 64,
            extra_spin_factor=1, additional_quadrature_factor=1,
            spatial_covariance_factor=ncells, partial_renormalization=False,
            mass_division=False, oscillator_normalization=False,
            momentum_state_covariance=False, electronic_generalized_force_only=True),
        raw_current_covariance=current['covariance'], covariance=covariance,
        raw_defect_policy='Report original, mapped original and newly reconstructed defects separately; never repair.',
        numerical_repairs_applied=False, independent_oracle_status='not_executed_by_consumer',
        scope=SCOPE)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--current-aggregate', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True,
                        help='new exclusive JSON path; parent directory must exist')
    args = parser.parse_args()
    try:
        result = evaluate(args.current_aggregate)
        code = 0 if result['status'] == 'complete' else 2
    except Exception as error:
        result = dict(schema=SCHEMA, status='invalid', error=repr(error), scope=SCOPE)
        code = 1
    try:
        receipt = save_json(args.output, result)
    except Exception as error:
        print(json.dumps(dict(event='canonical-current-noise-exit', status='invalid',
                              error=repr(error), exit_code=1), allow_nan=False),
              file=sys.stderr, flush=True)
        return 1
    print(json.dumps(dict(event='canonical-current-noise-exit', status=result['status'],
                          output=receipt, exit_code=code), allow_nan=False), flush=True)
    return code


if __name__ == '__main__':
    raise SystemExit(main())
