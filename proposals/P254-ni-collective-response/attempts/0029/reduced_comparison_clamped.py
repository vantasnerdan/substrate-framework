"""Canonical-dose clamped control from saved observables or the SAME scalar.

Saved one-trap-quantum deltas are linear means and can be scaled by J/Jold.
Finite kinetic and Maxwell energies are then evaluated with their OWN saved
moving baseline; they are never multiplied by (J/Jold)^2. Covariances are not
scaled at all. Only genuinely missing times use ClampedArray.evolve (quadratic
mechanical control, no new electronic solver, Native or DFT initialization).
"""
from pathlib import Path

import numpy as np
import clamped_scientific_readout as clamped
from substrate_framework.material_clamped_array import ClampedArray
from substrate_framework.streamed_material_state import file_binding, authenticate
from reduced_material_readout import spatial_pulse, finite_kinetic_exchange
from reduced_comparison_inputs import (
    array_archive, baseline_readout, bound_json, compact_spatial, decode, physical_context, require,
)


def discover_saved(index_path):
    """Expand the actual A28 compact-summary index; retain isotope/control scope."""
    index_binding = file_binding(Path(index_path))
    index = bound_json(index_binding)
    require(index['schema'] == 'P254-clamped-scientific-summary-index-v1'
            and index['status'] == 'complete', 'original clamped publication index required')
    publications = []
    for binding in index['completed'].values():
        summary = bound_json(binding)
        require(summary['schema'] == 'P254-clamped-scientific-compact-summary-v1'
                and summary['status'] == 'complete', 'completed original compact summary required')
        authenticate(summary['original_scientific_output'])
        for readout in summary['readouts']:
            inputs = bound_json(readout['original_input'])
            result = bound_json(readout['original_result'])
            require(result['inputs'] == readout['original_input'] and result['status'] == 'complete'
                    and result['schema'] == 'P254-matched-clamped-controls-v1', 'saved clamped producer/input mismatch')
            require(result['source_bindings'] == inputs['source_bindings']
                    and result['oracle'] == inputs['oracle'], 'saved clamped source/oracle differs')
            for source in inputs['source_bindings']:
                authenticate(source)
            oracle = bound_json(inputs['oracle'])
            require(oracle['status'] == 'passed' and oracle['source_bindings'] == inputs['source_bindings'],
                    'saved clamped original oracle/source differs')
            geometry = clamped.validate_identity(inputs, result)
            require(readout['case'] == result['case'] and readout['mass'] == result['mass']
                    and readout['control'] == result['control'], 'compact/source case or control differs')
            receipts = {}
            require(len(inputs['times_inverse_ry']) == len(result['completed_times']), 'clamped physical-time coverage differs')
            for time, receipt_binding in zip(inputs['times_inverse_ry'], result['completed_times']):
                receipt = bound_json(receipt_binding)
                require(receipt['status'] == 'complete' and receipt['time_inverse_ry'] == time
                        and receipt['reference_time_inverse_ry'] == inputs['reference_time_inverse_ry'],
                        'saved clamped receipt physical-time identity differs')
                require(time not in receipts, 'duplicate clamped physical time')
                receipts[time] = dict(receipt, receipt_binding=receipt_binding)
            publications.append(dict(index=index_binding, compact_summary=binding,
                original_publication=summary['original_scientific_output'], readout=readout,
                inputs=inputs, result=result, geometry=geometry, receipts=receipts))
    return publications


def matching_publications(run, publications):
    """Original/H1/H2 stay distinct. Independent-region controls are not coupled controls."""
    out = []
    spec, seed = run['specification'], run['seed']
    for publication in publications:
        item, inputs = publication['readout'], publication['inputs']
        identity, geometry = inputs['identity'], publication['geometry']
        if item['case'] != spec['case'] or item['control'] != 'periodic':
            continue
        require(physical_context(identity['matched_joint_specification']) == physical_context(spec['context']),
                'clamped/reduced physical source or original electronic preparation differs')
        require(inputs['kick_atom_original_one_based'] == seed['atom'] and inputs['kick_axis'] == seed['axis']
                and inputs['kick_cell'] == seed['cell'], 'clamped/reduced original seed differs')
        for key, actual in (('atom_indices', run['geometry']['ids']), ('atom_names', run['geometry']['names']),
                            ('cell_bohr', run['geometry']['cell']), ('valence_charges', run['geometry_arrays']['valence_charges'])):
            require(np.array_equal(clamped.decoded(identity[key]), actual), 'clamped/reduced geometry differs: '+key)
        require(np.array_equal(geometry['mesh'], run['geometry']['mesh']), 'clamped/reduced q coverage differs')
        if item['mass'] == 'original':
            prep = identity['preparation']
            for key, source in (('mass_ry', 'masses'), ('ionic_ordered', 'ionic_ordered'),
                                ('field_ordered', 'field_ordered'), ('trap_omega_ry', 'trap_omega')):
                require(np.array_equal(clamped.decoded(prep[key]), run['mechanical_arrays'][source]),
                        'clamped/reduced original mechanical preparation differs: '+key)
            require(clamped.real_metadata(prep['field_mass'], 'field mass').item()
                    == float(run['mechanical_arrays']['field_mass']), 'original Maxwell preparation mass differs')
        out.append(publication)
    return out


def q_from_ions_field(ions, field, geometry):
    d = geometry['d']
    ions = np.asarray(ions, complex)
    require(ions.shape == (*geometry['mesh'], 2*d) and np.shape(field) == (1, 6), 'periodic ions/global field layout required')
    q = np.fft.fftn(ions, axes=(0, 1, 2), norm='ortho').reshape(geometry['n'], 2*d)
    return (np.r_[q[0], np.asarray(field)[0]], *tuple(q[1:]))


def scalar_difference(model, baseline, delta):
    """Exact quadratic finite difference, retaining both raw transpose cross terms."""
    value = model.linear@delta[0]
    for q, h in enumerate(model.h):
        other = int(model.minus[q])
        value += (baseline[other]@h@delta[q] + delta[other]@h@baseline[q]
                  + delta[other]@h@delta[q])/2
    return value


def model_from_saved(publication):
    identity = publication['inputs']['identity']
    aggregate = bound_json(identity['aggregate'])
    require(aggregate['archive'] == identity['archive'], 'original clamped aggregate archive differs')
    a = array_archive(identity['archive'])
    for key in ('atom_indices', 'atom_names', 'cell_bohr', 'positions_bohr', 'valence_charges'):
        require(np.array_equal(clamped.decoded(identity[key]), a[key]), 'original aggregate geometry differs: '+key)
    require(np.array_equal(clamped.real_metadata(identity['mass']['original_masses_amu'], 'original masses'),
                           a['atomic_masses_amu']), 'original isotope authority differs')
    prep = decode(identity['preparation'])
    d = len(np.asarray(prep['mass_ry']))
    return ClampedArray(a['common_force_lattice_ry_bohr2'], np.asarray(prep['mass_ry']),
        a['valence_charges'], np.zeros((d, d)), a['contact'], np.asarray(a['mixed_ionic_field']).reshape(d, 3),
        np.asarray(prep['force_per_cell_ry_bohr']), np.asarray(prep['field_gradient_per_cell']), np.eye(3),
        float(np.real(prep['field_mass'])), np.asarray(prep['ionic_ordered']),
        np.asarray(prep['field_ordered']), np.asarray(prep['trap_omega_ry']))


def saved_correlation(arrays, receipt, geometry):
    n = geometry['n']
    minus = np.ravel_multi_index(tuple((-geometry['cells'] % geometry['mesh']).T), geometry['mesh'])
    equal_ab = arrays['spatial_equal_time_ordered'][:, 0, 0]
    equal_ba = arrays['spatial_equal_time_ordered'][minus, 0, 0]
    ab = arrays['spatial_two_time_ordered'][:, 0, 0]
    comm = arrays['spatial_two_time_commutator'][:, 0, 0]
    def components(a, b):
        require(np.shape(a) == np.shape(b) == (n,), 'clamped covariance cell coverage differs')
        return dict(ordered_ab=a, ordered_ba=b, symmetrized=(a+b)/2, commutator=a-b)
    return dict(available=True, reference_time_inverse_ry=receipt['reference_time_inverse_ry'],
        noise_scope='original clamped ionic/global-field product preparation; no electronic noise',
        equal_time=components(equal_ab, equal_ba), two_time=components(ab, ab-comm),
        units='bohr^2; selected original kick atom/axis displacement')


def evolved_correlation(model, evolution, reference, coordinate, geometry):
    ix = [coordinate]
    def read(left, right):
        ab = np.array([model.ordered_pair(left, right, cell, geometry['anchor'], ix)[0, 0]
                       for cell in geometry['cells']])
        ba = np.array([model.ordered_pair(right, left, geometry['anchor'], cell, ix)[0, 0]
                       for cell in geometry['cells']])
        return dict(ordered_ab=ab, ordered_ba=ba, symmetrized=(ab+ba)/2, commutator=ab-ba)
    return dict(available=True, reference_time_inverse_ry=reference.time_inverse_ry,
        noise_scope='original clamped ionic/global-field product preparation; no electronic noise',
        equal_time=read(evolution, evolution), two_time=read(evolution, reference),
        units='bohr^2; selected original kick atom/axis displacement')


def canonical_control(run, publication=None):
    """Return same-time canonical control; saved isotope variants are explicit counterfactuals."""
    geometry = run['geometry']
    if publication is None:
        model = ClampedArray(**run['mechanical_arrays'])
        variant, receipts = 'original', {}
        origin = dict(mechanical=run['manifest']['items']['mechanical']['archive'])
        scope = 'exact original matched preparation, electrons held fixed'
    else:
        model = model_from_saved(publication)
        variant, receipts = publication['readout']['mass'], publication['receipts']
        origin = {key: publication[key] for key in ('index', 'compact_summary', 'original_publication')}
        origin.update(original_result=publication['readout']['original_result'],
                      original_input=publication['readout']['original_input'])
        scope = ('exact original matched preparation, electrons held fixed' if variant == 'original' else
                 'isotope mass AND conditional vacuum counterfactual; no matching reduced isotope trajectory')
    seed = run['seed']
    old_amplitude = np.sqrt(2*model.masses[run['coordinate']]*model.omega[run['coordinate']])
    scale = seed['amplitude']/old_amplitude
    observations = []
    reference_cache = {}
    for actual in run['observations']:
        t = actual['time_inverse_ry']
        desired_reference = min(t, run['configuration']['correlation_reference'])
        receipt = receipts.get(t)
        if receipt is not None:
            arrays = array_archive(receipt['archive'])
            clamped.validate_observation(arrays, publication['geometry'])
            metadata = {key: dict(shape=list(v.shape), dtype=str(v.dtype)) for key, v in arrays.items()}
            require(metadata == receipt['arrays'], 'saved clamped NPZ/receipt metadata differs')
            baseline = q_from_ions_field(arrays['baseline_ions'], arrays['baseline_field'], geometry)
            raw_delta = q_from_ions_field(arrays['local_kick_delta_ions'], arrays['local_kick_delta_field'], geometry)
            delta = tuple(scale*x for x in raw_delta)
            correlation = saved_correlation(arrays, receipt, publication['geometry'])
            field_covariance = arrays['field_equal_time_ordered'][0]
            state_origin = dict(mode='reused original saved baseline, canonical-normalized linear mean and unscaled covariance',
                                receipt=receipt['receipt_binding'], archive=receipt['archive'])
        else:
            evolution = model.evolve(t)
            baseline = model.mean(evolution)['q']
            initial, _ = model.kicks(np.array(seed['cell'], int), run['coordinate'])
            delta = tuple(scale*(s@x) for s, x in zip(evolution.transitions, initial))
            if desired_reference not in reference_cache:
                reference_cache[desired_reference] = model.evolve(desired_reference)
            correlation = evolved_correlation(model, evolution, reference_cache[desired_reference], run['coordinate'], geometry)
            field_covariance = model.field_ordered_pair(evolution, evolution)
            state_origin = dict(mode='new canonical-unit response through existing ClampedArray mechanical scalar',
                                time_inverse_ry=t, reference_time_inverse_ry=desired_reference)
        # Saved controls have their own reference time. Do not relabel a saved
        # t,0 covariance as t,s when the reduced consumer requests another s.
        kinetic = finite_kinetic_exchange(baseline[0], delta, geometry, model.masses,
                                          run['geometry_arrays']['valence_charges'],
                                          float(run['mechanical_arrays']['field_mass']) if publication is None else
                                          float(np.real(clamped.decoded(publication['inputs']['identity']['preparation']['field_mass']))))
        field_mass = (float(run['mechanical_arrays']['field_mass']) if publication is None else
                      float(np.real(clamped.decoded(publication['inputs']['identity']['preparation']['field_mass']))))
        observations.append(dict(time_inverse_ry=t, source=state_origin,
            spatial=compact_spatial(spatial_pulse(delta, geometry), run['coordinate']),
            mechanical_energy_change_ry=scalar_difference(model, baseline, delta),
            finite_kinetic_exchange=kinetic,
            baseline=baseline_readout(baseline[0], geometry, model.masses, run['geometry_arrays']['valence_charges']),
            connected_correlation=correlation,
            maxwell_covariance=dict(equal_time_ordered=field_covariance,
                free_energy_ry=np.trace(field_covariance[3:, 3:])/(2*field_mass),
                coherent_impulse_change_ry=0j,
                scope='unscaled original clamped covariance; no electronic noise'),
            requested_correlation_reference_matches=correlation['reference_time_inverse_ry'] == desired_reference))
    return dict(case=run['specification']['case'], mass=variant, control='periodic', source=origin,
        physical_scope=scope, seed=seed, original_saved_impulse_bohr_inverse=old_amplitude,
        mean_amplitude_conversion=scale, finite_energy_rescaling_applied=False,
        covariance_amplitude_rescaling_applied=False, observations=observations)
