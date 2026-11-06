"""Source-resolved lattice panels and observed nested-model ranges.

Structures keep their own registered preparations. These are descriptive
conditional-model contrasts, never matched state subtractions, physical error
bounds, lifetime estimates or a calibrated lattice ranking.
"""
from itertools import combinations

import numpy as np

from reduced_comparison_inputs import CASES, digest

REPRESENTATIONS = ((8, 7), (8, 27), (12, 7), (12, 27))
REQUIRED_TIMES = (0., 1., 10., 30.)


def _signature(run):
    coordinate = run['coordinate']
    return dict(seed=run['seed'],
                selected_species=str(run['geometry_arrays']['atom_names'][coordinate // 3]),
                selected_mass=float(run['mechanical_arrays']['masses'][coordinate]),
                tolerances={key: run['configuration'][key] for key in ('rtol', 'atol')})


def _samples(curve):
    rows = []
    for sample in curve['samples']:
        row = dict(time_inverse_ry=sample['time_inverse_ry'], eligible=False)
        normalized = sample.get('ordered_coherence', {})
        if not sample['connected_available'] or not normalized.get('defined', False):
            row['reason'] = normalized.get('reason', 'connected coherence unavailable')
        elif sample['reference_time_inverse_ry'] != 0.:
            row['reason'] = 'different reference time; no replacement or interpolation'
        else:
            degree = np.asarray(normalized['complex_degree'])
            if degree.shape != (64, 1, 1):
                row['reason'] = 'requires the actual single selected probe on the complete64-cell grid'
            else:
                value = complex(degree[0, 0, 0])
                row.update(eligible=True, reference_time_inverse_ry=0.,
                           ordered_degree=value, magnitude=float(abs(value)),
                           phase_radians=float(np.angle(value)),
                           variance_left=normalized['variance_left'],
                           variance_right=normalized['variance_right'],
                           variance_real_parts_positive=normalized['variance_real_parts_positive'],
                           variance_imaginary_parts=normalized['variance_imaginary_parts'],
                           noise_scope=sample['noise_scope'])
        rows.append(row)
    return rows


def _envelope(case, time, slots):
    values, missing = [], []
    for bands, features in REPRESENTATIONS:
        records = slots.get((case, bands, features), [])
        if len(records) != 1:
            missing.append(dict(bands=bands, features=features,
                                state='absent' if not records else 'ambiguous'))
            continue
        sample = next((row for row in records[0]['samples']
                       if row['time_inverse_ry'] == time), None)
        if sample is None or not sample['eligible']:
            missing.append(dict(bands=bands, features=features,
                                state='missing_sample' if sample is None else 'undefined_sample',
                                reason='no exact-time sample' if sample is None else sample['reason']))
            continue
        values.append(sample['magnitude'])
    out = dict(case=case, time_inverse_ry=time, complete_observed_nested_set=not missing,
               observed_representations=len(values), missing=missing)
    if not missing:
        low, high = min(values), max(values)
        out.update(minimum_magnitude=low, maximum_magnitude=high,
                   observed_spread=high-low)
    return out


def _group_report(signature, slots):
    panels = []
    complete = True
    all_times = set(REQUIRED_TIMES)
    for bands, features in REPRESENTATIONS:
        cases = []
        for case in CASES:
            records = slots.get((case, bands, features), [])
            unique = len(records) == 1
            cases.append(dict(case=case, unique_actual_run=unique, actual_runs=records))
            if not unique:
                complete = False
                continue
            record = records[0]
            eligible = {row['time_inverse_ry'] for row in record['samples'] if row['eligible']}
            all_times.update(row['time_inverse_ry'] for row in record['samples'])
            complete = complete and record['status'] == 'complete' and set(REQUIRED_TIMES) <= eligible
        panels.append(dict(bands=bands, features=features, cases=cases))
    envelopes, contrasts = [], []
    for time in sorted(all_times):
        by_case = {case: _envelope(case, time, slots) for case in CASES}
        envelopes.extend(by_case.values())
        for left, right in combinations(CASES, 2):
            a, b = by_case[left], by_case[right]
            available = a['complete_observed_nested_set'] and b['complete_observed_nested_set']
            row = dict(left_case=left, right_case=right, time_inverse_ry=time,
                       complete_observed_nested_sets=available)
            if available:
                if a['minimum_magnitude'] > b['maximum_magnitude']:
                    relation = 'left_greater'
                elif b['minimum_magnitude'] > a['maximum_magnitude']:
                    relation = 'right_greater'
                else:
                    relation = 'overlap_or_touch'
                row['conditional_model_relation'] = relation
            else:
                row['reason'] = 'missing or ambiguous nested-model samples; no best-resolution substitution'
            contrasts.append(row)
    return dict(signature=signature, numerical_comparison_complete=complete,
                required_times_inverse_ry=list(REQUIRED_TIMES),
                case_panels=panels, observed_nested_ranges=envelopes,
                descriptive_pairwise_contrasts=contrasts)


def compare_lattices(runs, curves):
    """Compare like numerical representations without equating preparations."""
    groups = {}
    for run, curve in zip(runs, curves, strict=True):
        signature = _signature(run)
        group = groups.setdefault(digest(signature), dict(signature=signature, slots={}))
        case = run['specification']['case']
        representation = (case, run['configuration']['bands'], run['configuration']['features'])
        record = dict(source=run['binding'], status=run['status'],
                      original_case_identity=run['specification']['case_identity'],
                      original_model_source_identity=run['model_source_identity'],
                      own_preparation=True, samples=_samples(curve['coherence']))
        group['slots'].setdefault(representation, []).append(record)
    reports = [_group_report(group['signature'], group['slots']) for group in groups.values()]
    return dict(
        numerical_comparison_complete=any(group['numerical_comparison_complete'] for group in reports),
        selection='same registered canonical probe/species/mass/seed and solver tolerances; same representation and exact physical/reference times',
        observable='origin-cell selected ordered displacement coherence, normalized within each own moving unexcited baseline',
        comparison_groups=reports,
        physical_lattice_ranking=None, material_accuracy_certified=False,
        ranking_status='unresolved: different registered preparations and uncontrolled omitted full-PW/passive noise and geometry errors',
        scope='Descriptive conditional reduced-model comparisons only. Observed nested min/max are sensitivities, not error bounds or confidence intervals. No matched subtraction across structures, physical lifetime, entanglement or irreversible heat.')
