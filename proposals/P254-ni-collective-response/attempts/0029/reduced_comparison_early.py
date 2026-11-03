"""Assess original full-PW analytic responses, never guess a finite-dose state.

Optional early response uses the existing active-model q tangent and existing
SciPy DOP853 integration, on its own nonstationary active baseline. It is NOT a
second electronic solution, full-PW replay, or a second-order mixed solver.
Original mixed coordinates/kinetic response remain explicit unsupported targets.
"""
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from clamped_scientific_readout import encode
from substrate_framework.streamed_material_state import file_binding, authenticate
from reduced_comparison_inputs import (
    baseline_readout, bound_json, decode, digest, physical_context, require, same_payload,
)
from reduced_ode_initial_step import bounded_first_step


def difference(a, b):
    """b-a with a dimensionful absolute diagnostic; no acceptance threshold."""
    a, b = np.asarray(a), np.asarray(b)
    require(a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all(),
            'comparison observable layout/nonfinite differs')
    delta = b-a
    scale = float(np.max(np.abs(a))) if a.size else 0.
    return dict(difference=delta, max_abs_difference=float(np.max(np.abs(delta))) if delta.size else 0.,
                max_abs_left=scale, max_abs_right=float(np.max(np.abs(b))) if b.size else 0.,
                relative_to_left_max=None if scale == 0 else float(np.max(np.abs(delta)))/scale)


def read_full_pw(path):
    binding = file_binding(Path(path))
    report = bound_json(binding)
    require(report.get('schema') == 'P254-saved-canonical-Gamma-atom-response-v1'
            and report['complete'] is True, 'actual original saved canonical Gamma response required')
    require(report['q_fractional'] == [0., 0., 0.] and report['seed'] ==
            dict(atom=1, axis='x', amplitude=1., cell=[0, 0, 0]), 'original Gamma local impulse seed required')
    for source in report['source_bindings']:
        authenticate(source)
    for item in report['source_inputs'].values():
        if item is not None:
            authenticate(item)
    metadata = bound_json(report['source_inputs']['metadata'])
    require(report['original_metadata_identity'] == metadata['identity'], 'early report original metadata differs')
    require(digest(metadata['identity']['matched_joint_specification']) == report['model_identity'],
            'full-PW model specification identity differs')
    first = bound_json(report['source_inputs']['mean_initial'])
    last = bound_json(report['source_inputs']['mean_endpoint'])
    require(first['time_ry_inverse'] == report['actual_read_times_inverse_ry'][0] == 0.
            and last['time_ry_inverse'] == report['actual_read_times_inverse_ry'][-1]
            and first['original_electronic_boundary'] == last['original_electronic_boundary'],
            'original early physical-time/preparation identity differs')
    # Reuse the admitted saved-state small-coordinate reader. It authenticates
    # checkpoint/mean/source/D/seed bindings and never touches electronic legs.
    from saved_canonical_response_readout import read_saved_canonical_response
    reread = read_saved_canonical_response(**report['source_inputs'])
    require(reread['model_identity'] == report['model_identity'] and
            digest(encode(reread['response'])) == digest(report['response']),
            'published full-PW response differs from authenticated original small-coordinate readout')
    return dict(binding=binding, report=decode(report), raw=report, original_mean=first)


def early_compatibility(run, evidence):
    report, first = evidence['raw'], evidence['original_mean']
    if report['original_metadata_identity']['case'] != run['specification']['case']:
        return 'no original full-PW evidence for this structure'
    if physical_context(report['original_metadata_identity']['matched_joint_specification']) != physical_context(run['specification']['context']):
        return 'original full-PW and active-model physical context/preparation differ'
    if report['seed'] != run['seed']:
        return 'original full-PW and finite-run canonical impulse seeds differ'
    original = run['specification']['original_boundary']
    if not same_payload(first['original_electronic_boundary'], original):
        return 'original full-PW and active-model raw electronic preparation boundaries differ'
    raw = bound_json(original)
    slots = lambda x: [{k: v for k, v in row.items() if k != 'frame'} for row in x['slots']]
    if slots(raw) != slots(first):
        return 'original full-PW and active-model raw source/D/spin/k slots differ'
    from clamped_scientific_readout import decoded
    prep = report['original_metadata_identity']['preparation']
    for key, other in (('mass_ry', 'masses'), ('ionic_ordered', 'ionic_ordered'),
                       ('field_ordered', 'field_ordered'), ('trap_omega_ry', 'trap_omega')):
        if not np.array_equal(decoded(prep[key]), run['mechanical_arrays'][other]):
            return 'original full-PW and reduced mechanical preparation differs: '+key
    return None


def tangent_rows(run, times):
    """Two INDEPENDENT Gamma first directions; q0 J/sqrt(N), physical J/N.

    Neither original minus nor mixed is replaced by a conjugate or inferred
    zero. This evolves two explicitly seeded first directions; the mixed one
    requires second variation and is deliberately not fabricated.
    """
    from reduced_projection_cache import load_projected_coefficients
    from reduced_electronic_dynamics import ActiveBandDynamics
    from reduced_electronic_correlations import ActiveBandCorrelations
    config = run['configuration']
    c = load_projected_coefficients(run['manifest_path'], bands=config['bands'],
                                    features=config['features'], allow_partial=True)
    model = ActiveBandDynamics(c)
    rtol, atol = config['rtol'], config['atol']
    endpoint = max(times)
    require(endpoint > 0, 'positive actual early endpoint required')
    baseline = model.integrate_baseline((0., endpoint), rtol=rtol, atol=atol)
    tangent = ActiveBandCorrelations(c, baseline, rtol=rtol, atol=atol)
    initial = np.zeros((tangent.state_size(0), 2), complex)
    initial[c.d+run['coordinate'], 0] = run['seed']['amplitude']/np.sqrt(c.n)
    initial[c.d+run['coordinate'], 1] = run['seed']['amplitude']/np.sqrt(c.n)
    shape = initial.shape
    # The existing tangent evaluates every actual intermediate baseline, not a
    # stationary initial rho, and preserves all128 explicit slots.
    rhs = lambda time, value: tangent.tangent(0, time, value.reshape(shape)).ravel()
    vector = initial.ravel()
    first = bounded_first_step(rhs, (0., endpoint), vector, rtol=rtol, atol=atol)
    result = solve_ivp(rhs, (0., endpoint), vector, method='DOP853', rtol=rtol, atol=atol,
                       t_eval=times, first_step=first)
    require(result.success, 'existing early q-tangent integration failed: '+result.message)
    rows = []
    masses = c.mechanical.masses.reshape(-1, 3)
    charges = run['geometry_arrays']['valence_charges'][:, None]
    for index, time in enumerate(result.t):
        state = result.y[:, index].reshape(shape)
        require(np.isfinite(state).all(), 'nonfinite early tangent; no repair')
        directions = {}
        for column, label in enumerate(('plus', 'minus')):
            x = state[:c.mechanical_sizes[0], column]
            u = x[:c.d].reshape(-1, 3)/np.sqrt(c.n)
            p = x[c.d:2*c.d].reshape(-1, 3)/np.sqrt(c.n)
            a, pi = x[2*c.d:2*c.d+3], x[2*c.d+3:2*c.d+6]
            directions[label] = dict(u_bohr_xyz=u, p_bohr_inverse_xyz=p,
                velocity_bohr_ry_xyz=(p-charges*a[None, :])/masses,
                global_field_coordinate=a, global_field_momentum=pi)
        rows.append(dict(time_inverse_ry=float(time), directions=directions,
            baseline=baseline_readout(baseline.state(float(time)).x, run['geometry'],
                                      c.mechanical.masses, run['geometry_arrays']['valence_charges'])))
    return rows


def assess_early(run, evidence, *, calculate_tangent=False):
    reason = early_compatibility(run, evidence)
    report = evidence['report']
    base = dict(full_PW_source=evidence['binding'], finite_result=run['binding'],
        original_full_PW_model_identity=report['model_identity'],
        reduced_input_case_identity=run['specification']['case_identity'],
        definitions=report['conventions'], actual_full_PW_times_inverse_ry=report['actual_read_times_inverse_ry'],
        finite_dose_comparison=dict(supported=False,
            reason='plus/minus/mixed are independent analytic derivatives, not finite-dose endpoint states'),
        mixed_comparison=dict(supported=False,
            reason='first q tangent has no second mixed variation; mixed coordinates and N/2 kinetic product rule require mean*mixed and both independent first cross terms'),
        full_PW_mixed_response={label: dict(
            atoms=[dict(atom_original_one_based=row['atom_original_one_based'],
                        response=row['response']['mixed'], kinetic_mixed_global_ry=row['kinetic_mixed_global_ry'])
                   for row in report['response'][label]['atoms']],
            totals=report['response'][label]['totals']) for label in ('initial', 'endpoint')},
        material_calibration='two early Gamma samples of one structure do not calibrate finite64-q active-band dynamics, passive noise, late times or other structures')
    if reason is not None:
        return dict(base, supported=False, reason=reason)
    if not calculate_tangent:
        return dict(base, supported=False, reason='no like-for-like first-response calculation requested; use --early-tangent for existing q-tangent consumer')
    times = report['actual_read_times_inverse_ry']
    calculated = tangent_rows(run, times)
    rows = []
    for label, reduced in zip(('initial', 'endpoint'), calculated):
        original = report['response'][label]
        require(original['time_ry_inverse'] == reduced['time_inverse_ry'], 'early response time differs')
        comparisons = {}
        for direction in ('plus', 'minus'):
            comparisons[direction] = {key: difference(
                np.asarray([row['response'][direction][key] for row in original['atoms']]),
                reduced['directions'][direction][key]) for key in
                ('u_bohr_xyz', 'p_bohr_inverse_xyz', 'velocity_bohr_ry_xyz')}
            for key, target in (('coordinate', 'global_field_coordinate'), ('momentum', 'global_field_momentum')):
                comparisons[direction][target] = difference(original['global_field'][direction][key], reduced['directions'][direction][target])
        baseline_full = report['baseline_unexcited_mean'][label]
        baseline_difference = {
            'displacement': difference(baseline_full['u'], reduced['baseline']['displacement_by_atom_bohr']),
            'canonical_momentum': difference(baseline_full['p'], reduced['baseline']['canonical_momentum_by_atom_bohr_inverse']),
            'velocity': difference(baseline_full['velocity'], reduced['baseline']['velocity_by_atom_bohr_ry'])}
        rows.append(dict(time_inverse_ry=reduced['time_inverse_ry'], first_response=comparisons,
                         own_baseline_difference=baseline_difference))
    return dict(base, supported=True, supported_scope='Gamma first canonical response only; raw difference, no accuracy or finite-dose inference',
                comparison_sign='reduced minus original full-PW', comparisons=rows,
                method='existing128-slot q tangent on freshly integrated nonstationary reduced baseline; no original PW trajectory replay')
