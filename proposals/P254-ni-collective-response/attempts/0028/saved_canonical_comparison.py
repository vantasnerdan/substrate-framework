"""Same-time/seed Gamma coupled versus frozen-electron canonical response.

Uses each model's own baseline and the existing additive kinetic component.
Differences are raw conditional-model comparisons, not resolved screening,
finite-dose redistribution, heat, spatial transport or a lattice ranking.
"""
import numpy as np

import saved_canonical_response as algebra
from saved_canonical_response_readout import bound_json, geometry, same_binding


def decoded(value):
    return np.asarray(value['real']) + 1j*np.asarray(value['imag'])


def compare_response(coupled, control):
    """Consume actual API report and bound finite-control report, without evolution."""
    if (coupled['schema'] != 'P254-saved-canonical-Gamma-atom-response-v1'
            or coupled.get('complete') is not True
            or control['schema'] != 'P254-saved-clamped-Gamma-response-v1'
            or control['status'] != 'complete'
            or coupled['model_identity'] != control['model_identity']
            or coupled['seed'] != control['original_seed']
            or coupled['actual_read_times_inverse_ry'][-1] != control['time_ry_inverse']
            or any(coupled['q_fractional']) or any(control['q'])):
        raise ValueError('genuine same-model/time/Gamma/original-seed reports required')
    for key in ('mixed_initial', 'mean_initial', 'mean_endpoint', 'metadata'):
        if not same_binding(coupled['source_inputs'][key], control['source_inputs'][key]):
            raise ValueError('different original source input: '+key)
    g = geometry(bound_json(coupled['source_inputs']['metadata']))
    if g['ncells'] != control['ncells']:
        raise ValueError('different original cell count')
    small = control['endpoint']
    stages = dict(mean={k: decoded(v) for k, v in control['own_unexcited_baseline'].items()})
    stages.update({direction: {k: decoded(v) for k, v in small[direction].items()}
                   for direction in algebra.DIRECTIONS})
    mechanical = algebra.mechanical_stages(stages, masses=g['masses'],
        charges=g['charges'], polarizations=g['polarizations'])
    kinetic = algebra.kinetic_atomwise(mechanical, masses=g['masses'], ncells=g['ncells'])
    clamped = algebra.atom_table(mechanical, kinetic, atom_ids=g['atom_ids'],
        species=g['species'], positions=g['positions'])
    endpoint = coupled['response']['endpoint']
    rows = []
    for left, right in zip(endpoint['atoms'], clamped['atoms'], strict=True):
        if (left['atom_original_one_based'] != right['atom_original_one_based']
                or left['species'] != right['species']):
            raise ValueError('original atom/species ordering differs')
        directions = {}
        for direction in algebra.DIRECTIONS:
            a, b = left['response'][direction], right['response'][direction]
            directions[direction] = {
                key: dict(coupled=a[key], clamped=b[key], difference=a[key]-b[key])
                for key in ('u_bohr_xyz', 'velocity_bohr_ry_xyz')}
        rows.append(dict(atom_original_one_based=left['atom_original_one_based'],
            species=left['species'], response=directions,
            kinetic_mixed_global_ry=dict(coupled=left['kinetic_mixed_global_ry'],
                clamped=right['kinetic_mixed_global_ry'],
                difference=left['kinetic_mixed_global_ry']-right['kinetic_mixed_global_ry'])))
    totals = {}
    for species in ('all', 'Ni', 'H'):
        left, right = endpoint['totals'][species], clamped['totals'][species]
        totals[species] = dict(
            kinetic_mixed_global_ry=dict(coupled=left['kinetic_mixed_global_ry'],
                clamped=right['kinetic_mixed_global_ry'],
                difference=left['kinetic_mixed_global_ry']-right['kinetic_mixed_global_ry']),
            displacement_weight_bohr2={direction: dict(
                coupled=left['response'][direction]['displacement_weight_bohr2'],
                clamped=right['response'][direction]['displacement_weight_bohr2'])
                for direction in algebra.DIRECTIONS})
    return dict(schema='P254-matched-Gamma-canonical-comparison-v1', complete=True,
        model_identity=g['model'], time_ry_inverse=control['time_ry_inverse'],
        seed=coupled['seed'], q=[0, 0, 0], atoms=rows, totals=totals,
        clamped_endpoint=clamped, coupled_baseline=coupled['baseline_unexcited_mean']['endpoint'],
        clamped_baseline=mechanical['mean'],
        source_inputs=coupled['source_inputs'],
        interpretation='Raw coupled moving-electron versus original frozen-electron Gamma response under the SAME canonical seed and endpoint time. Excitation derivatives and each own baseline are separate. Differences are not stationary screening, finite-dose energy redistribution, resolved material effects, spatial transport, heat or lattice ranking; relevant numerical/model error remains unresolved.')
