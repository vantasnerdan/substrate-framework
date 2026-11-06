"""Compare ACTUAL A29 finite results with source-matched saved mechanical controls.

Required: --finite RESULT [RESULT ...] --output NEW_JSON. Each RESULT is an
actual run_reduced_material_dynamics.py result.json, with its manifest snapshot,
geometry binding, all selected baseline/finite NPZ boundaries and exact source
bytes (live or authenticated historical snapshots) available. Supply the four
structures and all (8,7),(8,27),(12,7),(12,27)
results to fill the coverage table; absent inputs/times stay explicitly missing.
No result pathname is guessed or asserted to exist.

By default --clamped-index and --full-pw consume the ACTUAL original A28 index
and canonical Gamma publication linked from P254/SOLUTION.md. Override their
paths to authenticated copies retaining their nested original bindings.
--early-tangent requests like-for-like FIRST Gamma response on the existing
active q tangent at the original full-PW sample times, not finite-dose matching
or original trajectory replay. Without it, the analytic/finite mismatch is
reported rather than replaced by a guessed comparison.
The inter-structure panel retains each original preparation. It compares only
the same canonical probe/species/mass/seed, representation and tolerances.
Observed nested ranges are descriptive sensitivities, not physical error bars.


Example (shell variables MUST name already-produced actual result files):
  python proposals/P254-ni-collective-response/attempts/0029/run_reduced_material_comparison.py \\
    --finite "$LOADED_8_7_RESULT" "$LOADED_8_27_RESULT" "$LOADED_12_7_RESULT" \\
             "$LOADED_12_27_RESULT" "$PRISTINE_RESULT" "$VACANCY_RESULT" "$DISPERSED_RESULT" \\
    --early-tangent --output "$NEW_COMPARISON_JSON"

Repeat the remaining structure's nested inputs in the same --finite list for
full coverage. Output creation is exclusive. No input/producer file is changed.
"""
import argparse
from itertools import combinations
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
A28 = HERE.parent/'0028'
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(A28)]

import numpy as np
from clamped_scientific_readout import encode
from substrate_framework.streamed_material_state import file_binding
from reduced_comparison_inputs import CASES, digest, read_finite_runs, require
from reduced_comparison_clamped import discover_saved, matching_publications, canonical_control
from reduced_comparison_early import assess_early, difference, read_full_pw
from reduced_coherence_readout import coherence_curve
from reduced_comparison_lattices import compare_lattices

LIMITS = [
    'Controlled dynamical active8/12 saved-band,7/27 density-feature approximation; not full-PW material accuracy.',
    'Every model retains its own nonstationary unexcited baseline. Baseline-subtracted finite means and analytic plus/minus/mixed derivatives are different observables.',
    'Nested and solver-tolerance differences are observed sensitivities, not rigorous error bars, convergence proofs or an experimental calibration.',
    'Fixed4^3 geometry, finite q sampling, omitted passive/complement response/noise and active-passive initial coherences have unresolved physical errors.',
    'Raw complex scalar energy and ordered/symmetric/commutator covariances are retained. No real-part projection, conjugate replacement or covariance repair.',
    'Finite energy ledgers are mean scalar differences. Saved active Gaussian correlations are linearized about the unexcited baseline, not the finite kicked state; no missing covariance energy is inferred.',
    'Participation and cell-centre extent measure displacement weight, not quantum coherence, entanglement or a physical lifetime.',
    'Connected active Gaussian fluctuations are not passive full-material noise or an entanglement witness.',
    'Reversible electronic/ionic/Maxwell exchange is not heat, dissipation, loss or an energy current.',
    'Different structures have different registered preparations; case-resolved rows are not matched subtractions across unlike preparations and imply no lattice ranking.',
    'H1/H2 clamped variants change inertial masses and conditional vacua. Without a reduced isotope trajectory they are explicit mechanical counterfactuals, not isotope validation of coupled dynamics.',
    'Independent-region publications change intercell coupling, field multiplicity and field vacuum; they are inventoried but not relabeled as the one-global-field matched control.',
    'No nuclear preparation, transition matrix element, branching fraction, energy-disposal mechanism or nuclear rate follows from these responses.'
]


def spatial_difference(left, right):
    require(left['selected_original_cartesian_coordinate'] == right['selected_original_cartesian_coordinate'],
            'selected original spatial probe differs')
    out = dict(selected_displacement_by_cell_complex_bohr=difference(
        left['selected_displacement_by_cell_complex_bohr'], right['selected_displacement_by_cell_complex_bohr']))
    for species in ('all_ions', 'Ni', 'H'):
        a, b = left[species], right[species]
        if not a['applicable'] or not b['applicable']:
            out[species] = dict(supported=False, reason='species absent in original structure')
            continue
        require(np.array_equal(a['original_atom_indices'], b['original_atom_indices']), 'species/original atom ordering differs')
        ma, mb = a['metric'], b['metric']
        out[species] = {key: (dict(supported=False, reason='zero displacement has undefined normalized metric')
                            if ma[key] is None or mb[key] is None else difference(ma[key], mb[key]))
                        for key in ('total_squared_weight', 'participation_cells', 'rms_distance_bohr',
                                    'squared_weight_by_cell')}
    return out


def covariance_difference(left, right):
    if not left['available'] or not right['available']:
        return dict(supported=False, reason='connected correlations not present on both sides')
    out = dict(supported=True, scope='raw difference between declared noise models, not a material-noise accuracy estimate',
               left_noise_scope=left['noise_scope'], right_noise_scope=right['noise_scope'])
    out['equal_time'] = {key: difference(left['equal_time'][key], right['equal_time'][key])
                         for key in ('ordered_ab', 'ordered_ba', 'symmetrized', 'commutator')}
    if left['reference_time_inverse_ry'] == right['reference_time_inverse_ry']:
        out['reference_time_inverse_ry'] = left['reference_time_inverse_ry']
        out['two_time'] = {key: difference(left['two_time'][key], right['two_time'][key])
                           for key in ('ordered_ab', 'ordered_ba', 'symmetrized', 'commutator')}
    else:
        out['two_time'] = dict(supported=False, reason='different actual reference physical times',
            left_reference_time_inverse_ry=left['reference_time_inverse_ry'],
            right_reference_time_inverse_ry=right['reference_time_inverse_ry'])
    return out


def kinetic_difference(left, right):
    species = sorted(set(left['ionic_by_species_ry']) & set(right['ionic_by_species_ry']))
    return dict(ionic_by_species_ry={key: difference(left['ionic_by_species_ry'][key], right['ionic_by_species_ry'][key])
                                     for key in species},
        ionic_by_cell_atom_ry=difference(left['ionic_by_cell_atom_ry'], right['ionic_by_cell_atom_ry']),
        ionic_total_ry=difference(left['ionic_total_ry'], right['ionic_total_ry']),
        maxwell_free_mean_change_ry=difference(left['maxwell_free_mean_change_ry'], right['maxwell_free_mean_change_ry']))


def compare_control(run, control):
    # canonical_control constructs this exact seed/time grid, not nearest times.
    require(control['seed'] == run['seed'], 'control canonical impulse differs')
    rows = []
    for a, b in zip(run['observations'], control['observations']):
        require(a['time_inverse_ry'] == b['time_inverse_ry'], 'control physical time differs')
        rows.append(dict(time_inverse_ry=a['time_inverse_ry'],
            spatial=spatial_difference(b['spatial'], a['spatial']),
            kinetic=kinetic_difference(b['finite_kinetic_exchange'], a['finite_kinetic_exchange']),
            shifted_total_mean_scalar_ry=difference(b['mechanical_energy_change_ry'], a['energy_change_ry']['total']),
            connected_correlation=covariance_difference(b['connected_correlation'], a['connected_correlation']),
            own_baseline_displacement_difference=difference(b['baseline']['displacement_by_atom_bohr'], a['baseline']['displacement_by_atom_bohr']),
            own_baseline_velocity_difference=difference(b['baseline']['velocity_by_atom_bohr_ry'], a['baseline']['velocity_by_atom_bohr_ry'])))
    return dict(finite_result=run['binding'], clamped_mass=control['mass'],
                supported_matched_preparation=control['mass'] == 'original',
                scope=control['physical_scope'], comparison_sign='reduced minus canonical clamped',
                observations=rows)


def sensitivity_pair(left, right):
    lc, rc = left['configuration'], right['configuration']
    if left['model_source_identity'] != right['model_source_identity'] or left['seed'] != right['seed']:
        return None
    changed = [key for key in ('bands', 'features', 'rtol', 'atol') if lc[key] != rc[key]]
    if not changed:
        kind = 'repeat-result sensitivity'
    elif changed == ['bands']:
        kind = 'nested active-band sensitivity at fixed density-feature family/tolerances'
    elif changed == ['features']:
        kind = 'nested density-feature sensitivity at fixed bands/tolerances'
    elif set(changed) <= {'rtol', 'atol'}:
        kind = 'solver-tolerance sensitivity at fixed electronic representation'
    else:
        kind = 'combined representation/tolerance change; causes not separately resolved'
    ltime = {row['time_inverse_ry']: row for row in left['observations']}
    rtime = {row['time_inverse_ry']: row for row in right['observations']}
    rows = []
    for time in sorted(set(ltime) & set(rtime)):
        a, b = ltime[time], rtime[time]
        rows.append(dict(time_inverse_ry=time, spatial=spatial_difference(a['spatial'], b['spatial']),
            energy_change_ry={key: difference(a['energy_change_ry'][key], b['energy_change_ry'][key])
                              for key in ('mechanical', 'band', 'interaction', 'screening', 'total')},
            kinetic=kinetic_difference(a['finite_kinetic_exchange'], b['finite_kinetic_exchange']),
            connected_correlation=covariance_difference(a['connected_correlation'], b['connected_correlation']),
            own_baseline_displacement_difference=difference(a['baseline']['displacement_by_atom_bohr'], b['baseline']['displacement_by_atom_bohr']),
            own_baseline_velocity_difference=difference(a['baseline']['velocity_by_atom_bohr_ry'], b['baseline']['velocity_by_atom_bohr_ry'])))
    return dict(kind=kind, changed=changed, left=left['binding'], right=right['binding'],
        comparison_sign='right minus left', observations=rows,
        left_only_times_inverse_ry=sorted(set(ltime)-set(rtime)),
        right_only_times_inverse_ry=sorted(set(rtime)-set(ltime)),
        scope='observed numerical sensitivity only; no tolerance gate or bound on omitted full-PW physics')


def energy_diagnostics(run):
    initial = next((row for row in run['observations'] if row['time_inverse_ry'] == 0.), None)
    expected = run['seed']['amplitude']**2/(2*run['mechanical_arrays']['masses'][run['coordinate']])
    rows = []
    for row in run['observations']:
        energy = row['energy_change_ry']
        rows.append(dict(time_inverse_ry=row['time_inverse_ry'],
            total_mean_difference_drift_from_saved_zero_ry=None if initial is None else
                energy['total']-initial['energy_change_ry']['total'],
            scalar_component_sum_defect_ry=energy['total']-sum(energy[k] for k in
                ('mechanical', 'band', 'interaction', 'screening')),
            imaginary_components_ry={key: np.imag(value) for key, value in energy.items()}))
    return dict(expected_centered_canonical_input_ry=expected,
        initial_input_scope='analytic J^2/(2 M_selected), centered original initial means; NOT a rescaling of finite kinetic energies',
        saved_time_zero_available=initial is not None,
        saved_zero_minus_expected_input_ry=None if initial is None else initial['energy_change_ry']['total']-expected,
        observations=rows,
        scope='energy consistency diagnostics, not a material accuracy bound or irreversible heat')




def coverage(runs):
    rows = []
    for case in CASES:
        for bands in (8, 12):
            for features in (7, 27):
                found = [run for run in runs if run['specification']['case'] == case
                         and run['configuration']['bands'] == bands and run['configuration']['features'] == features]
                rows.append(dict(case=case, bands=bands, features=features, available=bool(found),
                    results=[dict(source=run['binding'], status=run['status'], seed=run['seed'],
                                  times_inverse_ry=[x['time_inverse_ry'] for x in run['observations']],
                                  missing_times_inverse_ry=run['missing_times_inverse_ry']) for run in found]))
    return rows


def compare(args):
    runs = read_finite_runs(args.finite)
    require(len({(r['binding']['path'], r['binding']['sha256']) for r in runs}) == len(runs), 'duplicate finite result path')
    publications = discover_saved(args.clamped_index)
    evidence = [read_full_pw(path) for path in args.full_pw]
    controls, control_comparisons, early = [], [], []
    control_cache = {}
    for run in runs:
        matched = matching_publications(run, publications)
        variants = {p['readout']['mass']: p for p in matched}
        if not variants:
            variants = {'original': None}
        require('original' in variants, 'original-mass clamped control is required alongside isotope controls')
        for variant, publication in variants.items():
            key = digest(dict(case=run['specification']['case'], context=run['specification']['context'],
                mechanical=run['manifest']['items']['mechanical']['archive'], variant=variant,
                seed=run['seed'], times=[x['time_inverse_ry'] for x in run['observations']],
                reference=run['configuration']['correlation_reference']))
            if key not in control_cache:
                control_cache[key] = canonical_control(run, publication)
                controls.append(dict(identity=key, **control_cache[key]))
            if variant == 'original':
                item = compare_control(run, control_cache[key])
                item['clamped_control_identity'] = key
                control_comparisons.append(item)
            else:
                control_comparisons.append(dict(finite_result=run['binding'],
                    clamped_control_identity=key, clamped_mass=variant, supported=False,
                    reason='no reduced trajectory with this isotope mass AND conditional preparation; canonical isotope control is reported separately, not subtracted from a different preparation'))
        for original in evidence:
            if original['raw']['original_metadata_identity']['case'] == run['specification']['case']:
                early.append(assess_early(run, original, calculate_tangent=args.early_tangent))
        if not any(original['raw']['original_metadata_identity']['case'] == run['specification']['case'] for original in evidence):
            early.append(dict(finite_result=run['binding'], supported=False,
                              reason='no supplied original full-PW early canonical response for this structure'))
    sensitivities = [pair for a, b in combinations(runs, 2) if (pair := sensitivity_pair(a, b)) is not None]
    curves = [dict(case=run['specification']['case'], source=run['binding'],
        bands=run['configuration']['bands'], features=run['configuration']['features'],
        tolerances={key: run['configuration'][key] for key in ('rtol', 'atol')}, seed=run['seed'],
        original_matched_case_identity=run['specification']['case_identity'],
        original_coefficient_manifest_snapshot=run['manifest_binding'],
        status=run['status'], missing_times_inverse_ry=run['missing_times_inverse_ry'],
        observations=run['observations'], energy_diagnostics=energy_diagnostics(run),
        coherence=coherence_curve(run['observations'])) for run in runs]
    all_coverage = coverage(runs)
    lattices = compare_lattices(runs, curves)
    return dict(schema='P254-activeband-actual-matched-comparison-v1', status='complete',
        whole_campaign_complete=False, material_accuracy_certified=False, lattice_ranking=None,
        consumer_sources=[file_binding(HERE/name) for name in ('run_reduced_material_comparison.py',
            'reduced_coherence_readout.py', 'reduced_comparison_lattices.py',
            'reduced_comparison_inputs.py', 'reduced_comparison_clamped.py', 'reduced_comparison_early.py', 'reduced_ode_initial_step.py')],
        scope='actual source-bound finite data; exact canonical clamped response; supported early first-response assessment or explicit mismatch',
        authentication_scope='result/configuration digest, exact manifest snapshot, all recorded coefficient item archives once per CLI process, original source/context/code receipts, geometry, every consumed baseline/finite state, original clamped receipts/arrays and full-PW small-coordinate boundaries; no raw PW trajectory replay',
        source_recoveries=[dict(finite_result=run['binding'], **receipt)
                           for run in runs for receipt in run['source_recoveries']],
        source_recovery_scope='retrospective source recovery, not original execution binding; original saved result bindings are unchanged',
        units=dict(time='inverse Ry; exact saved physical times, no interpolation', energy='Ry', displacement='bohr', momentum='bohr^-1'),
        coverage=all_coverage,
        all_nested_four_structure_inputs_present=all(row['available'] for row in all_coverage),
        numerical_comparison_complete=lattices['numerical_comparison_complete'],
        inter_structure_comparison=lattices,
        actual_curves=curves, canonical_clamped_controls=controls, matched_clamped_comparisons=control_comparisons,
        numerical_sensitivities=sensitivities,
        sensitivity_gaps=['No matched sensitivity pair supplied for '+case for case in CASES
                          if not any(run['specification']['case'] == case and any(
                              p['left'] == run['binding'] or p['right'] == run['binding'] for p in sensitivities) for run in runs)],
        original_full_PW_early_assessments=early,
        saved_independent_region_controls=[dict(case=p['readout']['case'], mass=p['readout']['mass'],
            source=p['readout']['original_result'], used_as_matched_one_global_field_control=False,
            reason='changed intercell coupling, field multiplicity and conditional field vacuum')
            for p in publications if p['readout']['control'] == 'independent-regions'],
        physical_and_calibration_limits=LIMITS)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--finite', type=Path, nargs='+', required=True, help='existing actual A29 finite result.json paths; any missing case/representation is reported')
    parser.add_argument('--clamped-index', type=Path, default=A28/'clamped-scientific-summaries-001'/'result.json',
                        help='actual original compact-summary index, not a new run plan')
    parser.add_argument('--full-pw', type=Path, nargs='+', default=[A28/'saved-canonical-response-original-001.json'],
                        help='actual original canonical plus/minus/mixed Gamma publications')
    parser.add_argument('--early-tangent', action='store_true', help='compute existing reduced Gamma FIRST tangent at actual full-PW times; mixed remains unsupported')
    parser.add_argument('--output', type=Path, required=True, help='new exclusive output JSON; never overwrite prior science')
    args = parser.parse_args()
    require(not args.output.exists(), 'comparison output already exists')
    report = compare(args)
    with args.output.open('x') as stream:
        json.dump(encode(report), stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
    print(json.dumps(file_binding(args.output), sort_keys=True))


if __name__ == '__main__':
    main()
