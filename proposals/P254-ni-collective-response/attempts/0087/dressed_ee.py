"""Generic A60 dressed EE contractions, selecting only uncommitted cells.

No Gaussian U separation: both ordered electron callbacks retain all128 legs;
every selected exchange output still sums the original 64 inner h. Only the
last, external-output part of A46 direct_q is sliced on a resumed selection.
"""
import time
import numpy as np
from streamed_ee import _operands, _pair_tiles
from all_slot_coulomb_shared import static_superoperators, exchange_superoperator
from all_slot_coulomb_blas import direct_q
from all_slot_coulomb_contractions import validate_physical


def direct_selected(factors, q, forward, reverse, external):
    """A46 direct_q verbatim algebra, full bubble, selected final outputs.

The bubble's slot contraction is NOT truncated to external ownership. Slicing
only the final external legs prevents recomputing saved direct output cells.
    """
    c = factors.c
    external = np.asarray(external, dtype=int)
    if np.array_equal(external, np.arange(c.ns)):
        return direct_q(factors, q, forward, reverse)
    qm = int(c.minus[q])
    p = c.feature_partner[q][:len(factors.weights[q])]
    r, m = len(p), c.m
    a = c.target[qm, np.arange(c.ns)]
    y = c.density[qm][:, p].conj() @ forward[:, None]
    z = c.density[q][a].conj() @ reverse[a, None]
    yy = y.transpose(1, 0, 2, 3).reshape(r, -1)
    zz = z.swapaxes(-1, -2).transpose(1, 0, 2, 3).reshape(r, -1)
    w = factors.weights[q] / c.n
    kernel = (yy @ zz.T) * w[:, None] * w[None, :]
    b = c.target[qm, external]
    x = c.density[q][b].conj() @ forward[b, None]
    right = c.density[qm][external][:, p].conj().reshape(len(external), r, m*m)
    dressed = (kernel[None] @ right).reshape(len(external), r, m, m)
    return x.transpose(0, 2, 1, 3).reshape(len(external), m, r*m) @ dressed.reshape(len(external), r*m, m)


def stream_missing_ee(factors, electron, times, q, store, *, pair_batch=2,
                      external_batch=2, static_slots=64):
    """Persist each missing pair/external tile; committed outputs are excluded.

Static groups, full-slot operands and direct bubbles may be rebuilt on resume.
These are shared input contractions, not a license to replay saved outputs.
At most pair_batch*static_slots direct output cells and one external exchange
subtile can be uncommitted at a boundary; each delivered tile is durable.
    """
    c = factors.c
    validate_physical(factors)
    if (c.ns, c.m, c.nq, c.n) != (128, 12, 64, 64):
        raise ValueError('original all128/full12/inner64 inputs required')
    times = np.asarray(times, dtype=float)
    if not np.array_equal(times, np.arange(31, dtype=float)):
        raise ValueError('exact complete saved31 time grid required')
    if any(type(v) is not int or v <= 0 for v in (pair_batch, external_batch, static_slots)):
        raise ValueError('positive integer tile sizes required')
    run = store.record['runs'][-1]
    manifest = dict(schema='p254-a87-missing-cell-dressed-ee-v1', owned_outer_q=[q],
        exchange_inner_h=64, input_slots=128, external_slots=128, bands=12,
        reverse_legs='independent electron(a,b) and electron(b,a) callbacks',
        pair_batch=pair_batch, external_batch=external_batch, static_slots=static_slots,
        static_groups=0, operand_pair_evaluations=0, electron_callback_calls=0,
        direct_output_cells=0, exchange_output_cells=0, blocks=0,
        static_seconds=0., operands_seconds=0., direct_seconds=0., exchange_seconds=0.,
        static_payload_peak=0, operands_payload_peak=0, output_payload_peak=0,
        scientific_output_recomputed_for_committed_cells=0,
        shared_work_rebuild='static group, full128 operands and full128 direct bubble are rebuilt only for missing output work')
    run['ee_stream'] = manifest

    def counted_electron(a, b):
        manifest['electron_callback_calls'] += 1
        return electron(a, b)

    for begin in range(0, c.ns, static_slots):
        group = np.arange(begin, min(begin + static_slots, c.ns))
        if store.seen[:, :, group].all():
            continue
        tick = time.monotonic()
        matrix, b, cc = static_superoperators(factors, q, group)
        manifest['static_groups'] += 1
        manifest['static_seconds'] += time.monotonic() - tick
        manifest['static_payload_peak'] = max(manifest['static_payload_peak'], matrix.nbytes)
        missing_pairs = ((a, z) for a in range(31) for z in range(31)
                         if not store.seen[a, z, group].all())
        for tile in _pair_tiles(31, pair_batch, missing_pairs):
            # Different resume masks must never force a committed cell back into
            # a rectangular contraction. Group only identical missing patterns.
            patterns = {}
            for a, z in tile:
                positions = tuple(np.flatnonzero(~store.seen[a, z, group]).tolist())
                patterns.setdefault(positions, []).append((a, z))
            for positions, selected_pairs in patterns.items():
                if not positions:
                    continue
                pairs = tuple(selected_pairs)
                external = group[list(positions)]
                tick = time.monotonic()
                x, y, legs = _operands(factors, q, counted_electron, pairs)
                manifest['operands_seconds'] += time.monotonic() - tick
                manifest['operand_pair_evaluations'] += len(pairs)
                manifest['operands_payload_peak'] = max(manifest['operands_payload_peak'], x.nbytes + y.nbytes)
                tick = time.monotonic()
                direct = np.empty((2, len(external), len(pairs), c.m, c.m), dtype=np.complex128)
                for k, (forward, reverse) in enumerate(legs):
                    direct[0, :, k] = direct_selected(factors, q, forward[0], reverse[1], external)
                    direct[1, :, k] = direct_selected(factors, q, forward[1], reverse[0], external)
                manifest['direct_seconds'] += time.monotonic() - tick
                manifest['direct_output_cells'] += len(external) * len(pairs)
                for start in range(0, len(external), external_batch):
                    end = min(start + external_batch, len(external))
                    tick = time.monotonic()
                    # Static rows are external-major, with every original h.
                    selected = positions[start:end]
                    if selected[-1] - selected[0] + 1 == len(selected):
                        rows = slice(selected[0] * c.nq, (selected[-1] + 1) * c.nq)
                    else:
                        rows = (np.asarray(selected)[:, None] * c.nq
                                + np.arange(c.nq)[None, :]).reshape(-1)
                    xb = x[:, b[rows]].transpose(1, 0, 2, 3, 4)
                    yb = y[:, cc[rows]].transpose(1, 0, 2, 3, 4)
                    value = exchange_superoperator(matrix[rows], xb, yb)
                    value = value.reshape(end-start, c.nq, 2*len(pairs), c.m, c.m).sum(axis=1)
                    greater = direct[0, start:end] + value[:, 0::2]
                    lesser = direct[1, start:end] + value[:, 1::2]
                    mask = np.asarray([times[a] >= times[z] for a, z in pairs])[None, :, None, None]
                    values = dict(greater=greater, lesser=lesser, retarded=(greater-lesser)*mask)
                    manifest['exchange_seconds'] += time.monotonic() - tick
                    manifest['exchange_output_cells'] += (end-start) * len(pairs)
                    manifest['output_payload_peak'] = max(manifest['output_payload_peak'],
                        direct.nbytes + sum(v.nbytes for v in values.values()))
                    store.commit(pairs, external[start:end], values)
                    manifest['blocks'] += 1
                    del xb, yb, value, values, greater, lesser
                del x, y, legs, direct
        del matrix, b, cc
    return manifest


def oracle_direct_selected(factors, q, forward, reverse, *, rtol=1.e-12, atol=0.):
    """Unrun narrow row-selection oracle; caller supplies actual full128 legs.

    Compare noncontiguous/reordered and final-tail rows to original direct_q.
    A deliberately wrong-row selection must fail the SAME tolerance. An input
    unable to distinguish that mutation is rejected, not counted as evidence.
    Call separately for greater and lesser physical operand orientations.
    """
    if factors.c.ns != 128:
        raise ValueError('original all128 inputs required for the row oracle')
    if (not np.isfinite(rtol) or not np.isfinite(atol) or rtol < 0 or atol < 0):
        raise ValueError('finite nonnegative oracle tolerances required')
    reference = direct_q(factors, q, forward, reverse)
    if not np.isfinite(reference).all():
        raise ValueError('nonfinite original direct output')
    cases = {}
    for name, external in (
            ('noncontiguous_reordered', np.array([127, 3, 64, 1])),
            ('tail', np.array([125, 126, 127]))):
        actual = direct_selected(factors, q, forward, reverse, external)
        expected = reference[external]
        if not np.isfinite(actual).all() or not np.allclose(actual, expected, rtol=rtol, atol=atol):
            raise AssertionError('selected direct output differs: ' + name)
        wrong_rows = (external + 1) % 128
        mutated = reference[wrong_rows]
        if np.allclose(mutated, expected, rtol=rtol, atol=atol):
            raise AssertionError('input lacks wrong-row mutation sensitivity: ' + name)
        cases[name] = dict(external=external.tolist(), wrong_rows=wrong_rows.tolist(),
            max_abs_error=float(np.max(np.abs(actual-expected))),
            reference_max_abs=float(np.max(np.abs(expected))),
            wrong_row_max_abs_defect=float(np.max(np.abs(mutated-expected))),
            wrong_row_mutation_rejected=True)
    return dict(schema='p254-a87-selected-direct-row-oracle-v1',
        original_reference='A46 all_slot_coulomb_blas.direct_q',
        original_q=int(q), input_slots=128, rtol=rtol, atol=atol, cases=cases)
