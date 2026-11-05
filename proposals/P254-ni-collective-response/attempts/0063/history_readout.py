"""Original A53 observables plus an optional, separate A69 endpoint comparison."""
import json
import sys
import numpy as np
from history_storage import HERE, authenticate, binding
sys.path.append(str(HERE.parent/'0053'))
from consume_material_step import gamma_readout, json_value


def raw_difference(first, second):
    if first.shape != second.shape:
        raise ValueError('iterate dimensions differ')
    return max(float(np.max(np.abs(first[index]-second[index])))
               for index in np.ndindex(first.shape[:-2]))


def represented_control_readout(history):
    if history.control is None:
        return dict(status='unavailable',
            reason='No actual represented-control history supplied; not zero or a physical comparison.')
    return gamma_readout(history.control['control_ret'], history.control['control_less'],
                         int(history.aux['n']))


def projected_control_readout(history, path=None):
    """Read the actual original one-probe control; never construct a full bank."""
    if path is None:
        return dict(status='unavailable',
            reason='No actual projected-control path supplied; not zero or a full114 control.')
    sys.path.append(str(HERE.parent/'0069'))
    from projected_gamma_control import consume as consume_projected
    from substrate_framework.material_clamped_array import ClampedArray
    from saved_mean import POSITIONS

    control_binding = binding(path)
    record = json.loads(authenticate(control_binding).read_bytes())
    saved = consume_projected(control_binding['path'])
    reference_binding = getattr(history, 'reference_identity', None)
    if getattr(history, 'provider', None) is not None:
        reference_binding = history.provider.spec['input_identity']['numerical_reference']
    elif hasattr(history, 'spec'):
        reference_binding = history.spec['input_identity']['numerical_reference']
    if reference_binding != record['reference']:
        raise ValueError('projected control and history original numerical reference differ')
    reference = json.loads(authenticate(reference_binding).read_bytes())
    manifest_binding = record['input_bindings']['manifest']
    if (reference['coefficient_manifest'] != manifest_binding or
            record['baseline'] not in reference['physical_history_bindings']):
        raise ValueError('projected control original manifest/baseline differs')
    baseline = json.loads(authenticate(record['baseline']).read_bytes())
    manifest = json.loads(authenticate(manifest_binding).read_bytes())
    compatibility = json.loads(authenticate(record['source_compatibility']).read_bytes())
    if (compatibility['reference'] != reference_binding or
            compatibility['baseline'] != record['baseline'] or not compatibility['compatible']):
        raise ValueError('projected control original source compatibility differs')
    indexing = manifest['specification']['indexing']
    for key in ('k_indices', 'q_indices', 'slot_spin', 'slot_k', 'target', 'minus'):
        if indexing[key] != baseline['coefficient_specification']['indexing'][key]:
            raise ValueError('original projected-control source layout differs: '+key)
    preparation_bindings = []
    for key, item in record['input_bindings']['consumed_items'].items():
        if item['archive'] != manifest['items'][key]['archive']:
            raise ValueError('projected control consumed source differs: '+key)
        authenticate(item['archive'])
        preparation_bindings.append(item['archive'])
    for item in record['source_bindings']:
        authenticate(item)
    n = int(history.aux['n'])
    probe = dict(name='canonicalGammaNi1/x', index=0, normalization='1sqrtN',
                 coefficient=float(1/np.sqrt(n)), q=0)
    if (n != 64 or record['n'] != n or saved['probe'] != probe or
            record['initial_preparation']['mean_added'] is not False or
            record['dimensions']['propagated_probe_columns'] != 1):
        raise ValueError('original canonical probe/normalization/preparation differs')
    times = np.asarray(history.times)
    if times[0] != 0. or times[-1] != 30. or not np.all(np.diff(times) > 0):
        raise ValueError('literal genuine first/final saved nodes0/30 required; no interpolation')
    consumed = record['input_bindings']['consumed_items']
    for slot in range(128):
        with np.load(consumed['slot/%03d' % slot]['archive']['path'], allow_pickle=False) as raw:
            if not np.array_equal(1j*raw['rho0'], history.aux['source_electron_less0'][slot]):
                raise ValueError('original full128 coherent CAR preparation differs')
    # Reuse the owning original product-preparation assembly, not a new W0 rule.
    with np.load(consumed['mechanical']['archive']['path'], allow_pickle=False) as raw:
        mechanical = ClampedArray(**{key: raw[key] for key in raw.files})
    expected = -1j*mechanical.w0[0][np.ix_(POSITIONS, POSITIONS)].T
    if not np.array_equal(expected, history.aux['source_gamma_less'][0, 0]):
        raise ValueError('original full114 Gamma product preparation differs')
    ret, less = history.gamma_slot()
    if not np.array_equal(less[0, 0], expected):
        raise ValueError('history no longer has the original Gamma preparation')
    current = gamma_readout(ret, less, n)
    first, last = 1j*less[0, 0, 0, 0]/n, 1j*less[-1, -1, 0, 0]/n
    covariance = saved['raw_complex_covariance']
    if not np.isclose(covariance[0, 0], first, rtol=4*np.finfo(float).eps, atol=0.):
        raise ValueError('projected control original initial variance differs')
    return dict(status='available', eligible=True, control=control_binding,
        kind=record['kind'], probe=probe, times_inverse_ry=[0., 30.],
        history_saved_times_inverse_ry=times, history_endpoint_indices=[0, len(times)-1],
        n=n, history=dict(**current, initial_variance=first, final_variance=last),
        projected_control=dict(raw_complex_covariance=covariance,
            connected_cross=covariance[1, 0], initial_variance=covariance[0, 0],
            final_variance=covariance[1, 1], normalized=saved['normalized_ordered']),
        raw_connected_cross_difference=current['connected_cross']-covariance[1, 0],
        compatibility=dict(original_reference=reference_binding, baseline=record['baseline'],
            manifest=manifest_binding, source_compatibility=record['source_compatibility'],
            full128_original_rho0_equal=True, full114_original_position_W0_equal=True,
            consumed_archive_bindings=preparation_bindings),
        source_bindings=[control_binding, record['raw_archive'], *record['source_bindings']],
        interpretation='Original represented coupled Gaussian baseline at literal0/30 only; separate from legacy full114 controls. Not matched evolved-common-mean control, proof of intrinsic-only loss, temporal accuracy or physical closure. No control solve, full114 history or interpolation.')


def consume(history, previous, *, projected_control=None):
    nt = len(history.times)
    trace = np.zeros(nt, dtype=complex)
    minreal, maxreal, maximag, herm, prepared = np.inf, -np.inf, 0., 0., 0.
    delta_ret, delta_less = 0., 0.
    for slot in range(128):
        ret, less = history.electron_slot(slot)
        rho = -1j*less[np.arange(nt), np.arange(nt)]
        ev = np.linalg.eigvals(rho)
        trace += np.trace(rho, axis1=-2, axis2=-1)
        minreal, maxreal = min(minreal, ev.real.min()), max(maxreal, ev.real.max())
        maximag = max(maximag, float(np.abs(ev.imag).max()))
        herm = max(herm, float(np.abs(rho-rho.swapaxes(-1, -2).conj()).max()))
        prepared = max(prepared, float(np.abs(less[0, 0]-history.aux['source_electron_less0'][slot]).max()))
        old_ret, old_less = previous.electron_slot(slot)
        delta_ret = max(delta_ret, raw_difference(ret, old_ret))
        delta_less = max(delta_less, raw_difference(less, old_less))
    ret, less = history.gamma_slot()
    old_ret, old_less = previous.gamma_slot()
    n = int(history.aux['n'])
    return dict(schema='p254-common-history-readout-v1', times_inverse_ry=history.times,
        intrinsic=gamma_readout(ret, less, n),
        represented_control=represented_control_readout(history),
        projected_comparison=projected_control_readout(history, projected_control),
        original_bare=gamma_readout(history.aux['source_gamma_ret'], history.aux['source_gamma_less'], n),
        canonical_fermion_trace=trace,
        raw_fermion_eigenvalue_real_min=float(minreal), raw_fermion_eigenvalue_real_max=float(maxreal),
        raw_fermion_eigenvalue_imaginary_max=maximag, raw_fermion_hermiticity_max=herm,
        fermion_initial_source_maxabs=prepared,
        gamma_initial_source_maxabs=float(np.abs(less[0, 0]-history.aux['source_gamma_less'][0, 0]).max()),
        all128_canonical_fermions_consumed=128, all_Gamma_position_coordinates_consumed=114,
        delta_from_previous=dict(electron_ret_maxabs=delta_ret, electron_less_maxabs=delta_less,
            gamma_ret_raw_coordinate_maxabs=raw_difference(ret, old_ret),
            gamma_less_raw_coordinate_maxabs=raw_difference(less, old_less),
            interpretation='Raw iterate differences only; mixed Gamma units are not an error norm, conservation closure or convergence criterion.'),
        interpretation='Original unprojected A53 readout formulas. No refill, Hermitization or state repair; not temporal accuracy, calibrated loss, heat, rate or closed all-q energy conservation.')
