"""Gamma ionic111 refinement: five actual new pairs, four literal old pairs.

Only the launcher runs this module. Partial archives contain genuine running
internal-q sums, with per-pair coverage; they are never complete kernel tables.
The qualified A43/A47/A51 scientific proofs are inherited, not replayed.
"""
from pathlib import Path
import json
import sys
import time

import numpy as np

from refinement_common import HERE, TIMES, NEW_PAIRS, OLD_MAP, binding, save, catalog, reference

sys.path[:0] = [str(HERE.parent / name) for name in ('0041', '0051')]
from material_electron_ion_inputs import load_material
from cubic_contour import CubicContour, ordered_matrix, ordered_action

SOURCE = Path('/var/tmp/p254-a47-full37-source-001/source.json')
KEYS = ('greater', 'lesser', 'retarded')
PAYLOAD_CAP = 8 * 1024**2
NUMERIC_BYTES = 3 * 3 * 3 * 111 * 111 * np.dtype(complex).itemsize
FINALIZATION_RESERVE = 15.0


def _payload(output):
    return sum(p.stat().st_size for p in output.rglob('*') if p.is_file())


def _checkpoint(output, arrays, completed, timings, inputs):
    counts = np.full((3, 3), len(completed), dtype=np.int64)
    for a in OLD_MAP:
        for b in OLD_MAP:
            counts[a, b] = 64
    path = output / 'Gamma-Pi.npz'
    temporary = output / 'Gamma-Pi.tmp.npz'
    np.savez(temporary, **arrays, times_inverse_ry=TIMES,
             external_q=np.asarray([0]), external_rows=np.arange(111),
             completed_internal_q=np.asarray(completed, dtype=np.int64),
             new_time_pairs=np.asarray(NEW_PAIRS, dtype=np.int64),
             pair_internal_q_counts=counts, pair_complete=counts == 64)
    temporary.replace(path)
    scope = dict(
        schema='p254-a55-cubic-refinement-scope-v1',
        status='complete' if len(completed) == 64 else 'partial',
        archive=binding(path), times_inverse_ry=np.asarray(TIMES).tolist(),
        new_time_pairs=[list(p) for p in NEW_PAIRS],
        literal_old_pair_mapping={str(k): int(v) for k, v in OLD_MAP.items()},
        completed_internal_q=list(completed), required_internal_q=list(range(64)),
        pair_internal_q_counts=counts.tolist(), pair_complete=(counts == 64).tolist(),
        external_q=[0], external_coordinates=111, internal_coordinates=111,
        new_pair_semantics='exact running sum over completed_internal_q; complete only at all64',
        old_pair_semantics='literal archived full64 sums, never recomputed per q',
        per_q_new_pair_count=5, old_pairs_computed_per_q=0,
        unfinished_pairs_filled_with_zero=False, timings=timings, inputs=inputs,
        numeric_array_payload_bytes=int(NUMERIC_BYTES),
        component_persistent_payload_cap_bytes=PAYLOAD_CAP,
        atomic_checkpoint_peak_bound_bytes=2 * int(NUMERIC_BYTES) + 1024**2)
    save(output / 'scope.json', scope)
    if _payload(output) > PAYLOAD_CAP:
        raise RuntimeError('cubic component persistent payload cap exceeded')
    return scope


def _affected_orientation(left, right, dg1, dg2, matrix):
    """One finite new (0.5,0) q=0 boundary, not the inherited proof campaign."""
    x = np.linspace(.1, .9, 111) * (1 + .37j)
    y = x[::-1] * (.3 - .8j)
    action = ordered_action(left, right, dg1, dg2, x)
    dual = ordered_action(left, right, dg1, dg2, y, transpose=True)
    record = dict(
        scope='one actual greater contribution at external Gamma, internal q0, times (0.5,0)',
        new_time_pair=[1, 0], internal_q=0,
        matrix_action_error=float(np.max(np.abs(action - matrix @ x))),
        true_bilinear_transpose_error=float(np.max(np.abs(dual - matrix.T @ y))),
        wrong_Hilbert_adjoint_sensitivity=float(np.linalg.norm(dual - matrix.conj().T @ y)),
        actual_matrix_norm=float(np.linalg.norm(matrix)),
        inherited_oracles_replayed=False)
    numbers = [v for v in record.values() if isinstance(v, float)]
    if (not np.isfinite(numbers).all()
            or max(record['matrix_action_error'], record['true_bilinear_transpose_error']) > 2e-10
            or record['wrong_Hilbert_adjoint_sensitivity'] <= 0):
        raise AssertionError(record)
    return record


def _consume(output, old_binding):
    scope = json.loads((output / 'scope.json').read_bytes())
    if binding(output / 'Gamma-Pi.npz') != scope['archive']:
        raise ValueError('saved cubic archive binding changed')
    with np.load(output / 'Gamma-Pi.npz', allow_pickle=False) as saved:
        arrays = {key: saved[key] for key in KEYS}
        times = saved['times_inverse_ry']
        counts = saved['pair_internal_q_counts']
        complete = saved['pair_complete']
        completed = saved['completed_internal_q']
    if (not np.array_equal(times, TIMES) or counts.shape != (3, 3)
            or not np.array_equal(complete, counts == 64)
            or any(a.shape != (3, 3, 111, 111) or not np.isfinite(a).all()
                   for a in arrays.values())):
        raise ValueError('saved full3x3 ionic111 arrays or explicit pair coverage differ')
    if binding(old_binding['path']) != old_binding:
        raise ValueError('primary old cubic archive binding changed')
    with np.load(old_binding['path'], allow_pickle=False) as old:
        unchanged = all(arrays[key][a, b].tobytes() == old[key][oa, ob].tobytes()
                        for key in KEYS for a, oa in OLD_MAP.items()
                        for b, ob in OLD_MAP.items())
    theta = (times[:, None] >= times[None, :])[:, :, None, None]
    causal_error = float(np.max(np.abs(arrays['retarded'] - theta *
                                       (arrays['greater'] - arrays['lesser']))))
    # Original Gamma coordinate ordering is ionic111 followed by field3.
    # Structural zero field blocks are an embedding, not invented cubic values.
    embedded = np.zeros((114, 114), dtype=complex)
    embedded[:111, :111] = arrays['retarded'][1, 0]
    x = np.linspace(.1, .9, 114) * (1 + .37j)
    y = x[::-1] * (.3 - .8j)
    action = embedded @ x
    dual = embedded.T @ y
    record = dict(
        status='passed', scope=scope['status'], archive=scope['archive'],
        primary_old_archive=old_binding, literal_oldfour_byte_identical=unchanged,
        times_inverse_ry=times.tolist(), completed_internal_q=completed.tolist(),
        pair_internal_q_counts=counts.tolist(), pair_complete=complete.tolist(),
        causal_saved_array_error=causal_error,
        actual_greater_norm=float(np.linalg.norm(arrays['greater'])),
        actual_lesser_norm=float(np.linalg.norm(arrays['lesser'])),
        actual_retarded_norm=float(np.linalg.norm(arrays['retarded'])),
        actual_new_pair_retarded_action_norm=float(np.linalg.norm(action)),
        actual_new_pair_noise_norm=float(np.linalg.norm(.5j *
            (arrays['greater'][1, 0] + arrays['lesser'][1, 0]))),
        original114_embedding_bilinear_error=float(abs(y @ action - dual @ x)),
        field_output_norm=float(np.linalg.norm(action[111:])),
        actual_new_pair_complex_sample=[float(arrays['greater'][1, 0, 0, 93].real),
                                        float(arrays['greater'][1, 0, 0, 93].imag)],
        consumer_pair=[1, 0], consumer_pair_complete=bool(complete[1, 0]),
        action_scope='full64 sum' if complete[1, 0] else 'genuine partial internal-q sum only',
        ionic_block_only=True, field_cubic_terms_added=False, coupled_evolution=False)
    if (not unchanged or causal_error > 2e-12
            or record['original114_embedding_bilinear_error'] > 2e-10
            or record['field_output_norm'] != 0):
        raise AssertionError(record)
    save(output / 'saved-consumer.json', record)
    return record


def run(output, wall_cap=388.713):
    """Compute bounded genuine work; return explicit complete/partial evidence."""
    started = time.monotonic()
    wall_cap = float(wall_cap)
    if not np.isfinite(wall_cap) or wall_cap <= 0:
        raise ValueError('finite positive cubic component wall cap required')
    output = Path(output).resolve()
    output.relative_to((HERE / 'material').resolve())
    output.mkdir(parents=True, exist_ok=False)
    if (not np.array_equal(TIMES, [0., .5, 1.])
            or list(map(tuple, NEW_PAIRS)) != [(0, 1), (1, 0), (1, 1), (1, 2), (2, 1)]
            or OLD_MAP != {0: 0, 2: 1}):
        raise ValueError('three-node cubic refinement contract differs')
    old_binding = catalog()['archives']['cubic']
    if binding(old_binding['path']) != old_binding:
        raise ValueError('primary old cubic archive failed byte authentication')
    arrays = {key: np.full((3, 3, 111, 111), np.nan + 1j * np.nan,
                           dtype=complex) for key in KEYS}
    with np.load(old_binding['path'], allow_pickle=False) as old:
        if (not np.array_equal(old['times_inverse_ry'][:2], [0., 1.])
                or not np.array_equal(old['completed_internal_q'], np.arange(64))
                or not np.array_equal(old['external_rows'], np.arange(111))
                or not np.array_equal(old['external_q'], [0])):
            raise ValueError('old cubic archive is not full Gamma111/all64 at original0/1')
        for key in KEYS:
            if old[key].shape != (4, 4, 111, 111):
                raise ValueError('old cubic archive axes differ')
            for a, oa in OLD_MAP.items():
                for b, ob in OLD_MAP.items():
                    arrays[key][a, b] = old[key][oa, ob]
    c, manifest, coefficient_binding = load_material()
    provider = CubicContour.from_source(SOURCE)
    source_record = provider.source.source_record
    if (provider.d, provider.n, c.d, c.f) != (111, 64, 111, 3):
        raise ValueError('full original cubic and Gamma axes required')
    if source_record['coefficient_manifest'] != coefficient_binding:
        raise ValueError('compact cubic source and actual material bindings differ')
    ref = reference(c, coefficient_binding)
    timings, decisions, completed = [], [], []
    inputs = dict(
        primary_old_archive=old_binding,
        original_catalog=binding(HERE.parent / '0053' / 'actual-input-catalog-001.json'),
        compact_source=binding(SOURCE), compact_source_arrays=source_record['arrays'],
        compact_source_implementation_bindings=source_record['source_bindings'],
        coefficient_manifest=coefficient_binding, actual_reference_bindings=ref.bindings,
        implementation_bindings=[binding(p) for p in (
            Path(__file__), HERE / 'refinement_common.py',
            HERE.parent / '0051' / 'cubic_contour.py',
            HERE.parent / '0047' / 'cubic_feedback.py',
            HERE.parent / '0052' / 'general_moving_reference.py',
            HERE.parent / '0041' / 'material_electron_ion_inputs.py')])
    affected = None
    try:
        if not np.array_equal(ref.times, TIMES):
            raise ValueError('exact requested-time lazy reference differs')
        for k in range(64):
            unit = time.monotonic()
            l, left, right = provider.pair(0, k)
            vertex_seconds = time.monotonic() - unit
            for a, b in NEW_PAIRS:
                dg1, dl1 = ref.boson(k, a, b)
                dg2, dl2 = (dg1, dl1) if l == k else ref.boson(l, a, b)
                for key, first, second in (('greater', dg1, dg2), ('lesser', dl1, dl2)):
                    contribution = ordered_matrix(left, right, first, second)
                    if k == 0:
                        arrays[key][a, b] = contribution
                    else:
                        arrays[key][a, b] += contribution
                    if k == 0 and (a, b) == (1, 0) and key == 'greater':
                        affected = _affected_orientation(left, right, first, second, contribution)
                arrays['retarded'][a, b] = ((arrays['greater'][a, b] - arrays['lesser'][a, b])
                                           if TIMES[a] >= TIMES[b] else 0.)
            completed.append(k)
            seconds = time.monotonic() - unit
            timings.append(dict(internal_q=k, routed_partner_q=int(l), seconds=seconds,
                                vertex_seconds=vertex_seconds, new_time_pair_count=5,
                                old_time_pair_count=0))
            scope = _checkpoint(output, arrays, completed, timings, inputs)
            elapsed = time.monotonic() - started
            peak = max(t['seconds'] for t in timings)
            forecast = elapsed + (64 - len(completed)) * peak + FINALIZATION_RESERVE
            if k == 0:
                first = dict(first_actual_internal_q=timings[0], elapsed_seconds=elapsed,
                             measured_q_seconds=seconds, remaining_internal_q=63,
                             forecast_seconds=forecast, wall_cap_seconds=wall_cap,
                             fits=forecast <= wall_cap, forecast_basis='elapsed + remaining64q * measured_peak_q + 15s finalization',
                             action='continue' if forecast <= wall_cap else 'preserve genuine partial; full64 forecast exceeds cubic grant')
                decisions.append(first)
                save(output / 'first-unit-decision.json', first)
                save(output / 'affected-orientation.json', affected)
                if not first['fits']:
                    break
            if len(completed) < 64 and elapsed + max(1.5 * peak, 10.) + FINALIZATION_RESERVE > wall_cap:
                decisions.append(dict(action='preserve genuine partial before next indivisible q',
                                      elapsed_seconds=elapsed, completed_internal_q=list(completed),
                                      next_q_reserve_seconds=max(1.5 * peak, 10.),
                                      forecast_seconds=forecast, wall_cap_seconds=wall_cap))
                break
    finally:
        ref.close()
    consumer = _consume(output, old_binding)
    result = dict(
        schema='p254-a55-gamma-cubic-three-node-refinement-v1', status=scope['status'],
        archive=scope['archive'], scope=binding(output / 'scope.json'),
        actual_saved_consumer=binding(output / 'saved-consumer.json'), consumer=consumer,
        affected_comparison=binding(output / 'affected-orientation.json'),
        inputs=inputs, completed_internal_q=completed, required_internal_q=list(range(64)),
        new_time_pairs=[list(p) for p in NEW_PAIRS], times_inverse_ry=np.asarray(TIMES).tolist(),
        old_pairs_computed_per_q=0, internal_coordinate_count=111, external_q=[0],
        decisions=decisions, timings=timings, elapsed_seconds=time.monotonic() - started,
        measured_q_seconds=[t['seconds'] for t in timings],
        latest_full64_forecast_seconds=forecast,
        numeric_array_payload_bytes=int(NUMERIC_BYTES),
        component_persistent_payload_cap_bytes=PAYLOAD_CAP,
        atomic_checkpoint_peak_bound_bytes=2 * int(NUMERIC_BYTES) + 1024**2,
        payload_bytes_before_result=_payload(output),
        source_proofs='A43/A47/A51 scientific proofs inherited unchanged; no old proof oracle replay',
        conventions=dict(D='-i ordered connected Cartesian covariance', Pi='-i connected cubic force covariance',
                         theta_at_zero=1, canonical_N='1/sqrt64 in each full111 vertex; no extra weight'),
        limits=['Only Gamma external ionic111 cubic block; field3 embedded structurally with no invented field cubic',
                'Five midpoint-affected pairs only; oldfour literal archived full64 sums',
                'Partial new pairs are only the explicitly listed internal-q sums, not completed kernels',
                'No coupled evolution, new scientific proof, physical accuracy, damping or lifetime claim'])
    save(output / 'result.json', result)
    if _payload(output) > PAYLOAD_CAP:
        raise RuntimeError('cubic component persistent payload cap exceeded')
    return result
