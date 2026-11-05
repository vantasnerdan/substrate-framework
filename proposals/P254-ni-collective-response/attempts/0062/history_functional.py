"""Arbitrary-history exact A56 convention, reduced directly into owned NPY banks.

API: produce(output, bundle=...|dataset=...|provider=..., times=None,
             owned_q=None, owned_pairs=None, partial=False, pair_batch=2,
             external_batch=2, static_slots=16).
No production allocation/launch is hidden in this callable.
"""
from pathlib import Path
import sys
import time
import resource
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE)] + [str(HERE.parent / n) for n in
    ('0063','0060','0056','0055','0054','0053','0052','0051','0050','0048','0046','0045','0044','0043','0042','0041','0030','0029','0028','0027','0025','0024')]
import numpy as np
from history_banks import FinalBanks, ORDERS, binding
from history_sources import open_reference
from streamed_ee import stream_ee, _pair_tiles
from dressed_functionals import common_pair, SOURCE
from electron_ion_functional import evaluate_components
from material_electron_ion_inputs import edge_pair
from joint_material_inputs import load_inputs
from coulomb_material_inputs import load_coulomb_inputs
from bare_coulomb_factors import BareCoulombFactors
from cubic_contour import CubicContour, ordered_matrix


def produce(output, *, bundle=None, dataset=None, provider=None, times=None,
            owned_q=None, owned_pairs=None, partial=False, pair_batch=2,
            external_batch=2, static_slots=16):
    """Compute all components on declared q/pair ownership, with exact inner64.

    Partial sums retain full physical slot/band/coordinate dimensions and an
    exact global time grid. Unowned zeros are storage only, never full kernels.
    Existing owning evaluators implement every contraction, including A56's
    joint-minus-ionic Sigma and unstacked common Pi. Cubic is freshly summed
    over owned q, equivalent to A56's paid-bare-minus-old0-plus-new0 convention.
    """
    started = time.monotonic()
    qs = tuple(range(64)) if owned_q is None else tuple(owned_q)
    if not qs or len(set(qs)) != len(qs) or any(type(q) is not int or not 0 <= q < 64 for q in qs):
        raise ValueError('nonempty unique original q ownership required')
    if any(type(n) is not int or n <= 0 for n in (pair_batch, external_batch, static_slots)):
        raise ValueError('positive integer tile sizes required')
    if not partial and set(qs) != set(range(64)):
        raise ValueError('q subset requires explicit partial=True')
    c, baseline, old, material_bindings = load_inputs()
    old.close()
    ref = open_reference(c, material_bindings['manifest'], bundle=bundle,
                         dataset=dataset, provider=provider, times=times)
    try:
        nt = len(ref.times)
        pairs = None if owned_pairs is None else tuple(tuple(p) for p in owned_pairs)
        if pairs is not None and (not pairs or len(set(pairs)) != len(pairs)):
            raise ValueError('nonempty unique pair ownership required')
        # Validate before output allocation/numeric side effects.
        for tile in _pair_tiles(nt, pair_batch, pairs):
            for a, b in tile:
                if type(a) is not int or type(b) is not int:
                    raise ValueError('integer time indices required')
        full_pairs = pairs is None or set(pairs) == {(a,b) for a in range(nt) for b in range(nt)}
        if not partial and not full_pairs:
            raise ValueError('pair subset requires explicit partial=True')
        pair_scope = 'all' if full_pairs else [list(p) for p in pairs]
        coverage = {name: dict(owned_q=list(qs), completed_q=[], owned_pairs=pair_scope,
                              completed_pairs=[]) for name in ('ee','mixed_ei','cubic')}
        coverage['joint'] = dict(owned_pairs=pair_scope, completed_pairs=[])
        source = dataset if dataset is not None else provider if provider is not None else bundle
        record = dict(times_inverse_ry=ref.times.tolist(), input_dataset=binding(source),
            source_bindings=ref.bindings, coefficient_bindings=material_bindings,
            implementation_bindings=[binding(p) for p in sorted(HERE.glob('*.py'))],
            component_implementation_bindings=[binding(Path(__import__(name).__file__)) for name in
                ('streamed_ee','dressed_functionals','electron_ion_functional','material_electron_ion_inputs',
                 'joint_gamma_adapter','cubic_contour','general_moving_reference','bare_coulomb_factors')],
            coverage=coverage, slots=128, bands=12, gamma_coordinates=114,
            global_q=64, ee_exchange_inner_h=64, partial_requested=partial,
            combination='EE total + all64 raw EI Sigma + joint Gamma Sigma - ionic Gamma Sigma; Pi = full common joint + cubic ionic block; represented control separate',
            pair_batch=pair_batch, external_batch=external_batch, static_slots=static_slots,
            axes=dict(sigma=['time','time','slot','band','band'], pi=['time','time','coordinate','coordinate']),
            timings=[])
        with FinalBanks(output, ref.times) as owned:
            sigma, pi = owned.banks['sigma'], owned.banks['pi']
            # Coulomb's all-q density bank is the existing bounded static input,
            # not a history evaluator operand, and is released before EI decode.
            began = time.monotonic()
            factors = BareCoulombFactors(load_coulomb_inputs())
            record['ee_coefficient_bindings'] = factors.c.bindings
            def consume(block):
                for j, (a,b) in enumerate(block.pairs):
                    for order in ORDERS:
                        sigma[order][a,b,list(block.external)] += block.values[order + '_total'][:,j]
            record['ee_stream'] = stream_ee(factors, ref.electron, ref.times, consume,
                owned_q=qs, owned_pairs=pairs, pair_batch=pair_batch,
                external_batch=external_batch, static_slots=static_slots,
                source_bindings=ref.bindings)
            del factors
            coverage['ee']['completed_q'] = list(qs)
            coverage['ee']['completed_pairs'] = pair_scope
            record['timings'].append(dict(component='ee', seconds=time.monotonic()-began))
            for q in qs:
                began = time.monotonic()
                data = c.q(q)
                for tile in _pair_tiles(nt, pair_batch, pairs):
                    for a,b in tile:
                        forward, reverse = ref.electron(a,b), ref.electron(b,a)
                        boson = ref.boson(q,a,b)
                        for source_slot in range(c.ns):
                            target = int(c.target[q,source_slot])
                            pair = edge_pair(data, source_slot, target, forward, reverse, boson)
                            values = evaluate_components(pair, n=c.n, t=ref.times[a], s=ref.times[b], components=('sigma',))
                            for order in ORDERS:
                                sigma[order][a,b,target] += values['sigma_' + order]
                coverage['mixed_ei']['completed_q'].append(q)
                record['timings'].append(dict(component='mixed_ei', q=q, seconds=time.monotonic()-began))
            coverage['mixed_ei']['completed_pairs'] = pair_scope
            began = time.monotonic()
            for tile in _pair_tiles(nt, pair_batch, pairs):
                for a,b in tile:
                    full = common_pair(c,ref,a,b)
                    ionic = common_pair(c,ref,a,b,True)
                    for order in ORDERS:
                        sigma[order][a,b] += full['sigma_' + order] - ionic['sigma_' + order]
                        pi[order][a,b] += full['pi_' + order]
            coverage['joint']['completed_pairs'] = pair_scope
            record['timings'].append(dict(component='joint', seconds=time.monotonic()-began))
            cubic = CubicContour.from_source(SOURCE)
            record['cubic_source'] = binding(SOURCE)
            record['cubic_source_arrays'] = cubic.source.source_record['arrays']
            for q in qs:
                began = time.monotonic()
                partner,left,right = cubic.pair(0,q)
                for tile in _pair_tiles(nt, pair_batch, pairs):
                    for a,b in tile:
                        first, second = ref.boson(q,a,b), ref.boson(partner,a,b)
                        greater = ordered_matrix(left,right,first[0],second[0])
                        lesser = ordered_matrix(left,right,first[1],second[1])
                        pi['greater'][a,b,:111,:111] += greater
                        pi['lesser'][a,b,:111,:111] += lesser
                        if ref.times[a] >= ref.times[b]:
                            pi['retarded'][a,b,:111,:111] += greater - lesser
                coverage['cubic']['completed_q'].append(q)
                record['timings'].append(dict(component='cubic', q=q, partner=int(partner), seconds=time.monotonic()-began))
                del left,right
            coverage['cubic']['completed_pairs'] = pair_scope
            record['elapsed_seconds'] = time.monotonic() - started
            record['peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
            return owned.publish(record)
    finally:
        ref.close()


def main():
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument('--bundle', type=Path)
    source.add_argument('--dataset', type=Path)
    source.add_argument('--provider', type=Path)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--times', type=float, nargs='+')
    p.add_argument('--q', type=int, nargs='+')
    p.add_argument('--pairs', help='JSON ordered [a,b] index list; no inferred reverse legs')
    p.add_argument('--partial', action='store_true')
    p.add_argument('--pair-batch', type=int, default=2)
    p.add_argument('--external-batch', type=int, default=2)
    p.add_argument('--static-slots', type=int, default=16)
    args = p.parse_args()
    import json
    record = produce(args.output, bundle=args.bundle, dataset=args.dataset, provider=args.provider,
        times=args.times, owned_q=args.q, owned_pairs=None if args.pairs is None else json.loads(args.pairs),
        partial=args.partial, pair_batch=args.pair_batch, external_batch=args.external_batch, static_slots=args.static_slots)
    print(json.dumps(dict(status=record['status'], output=str(args.output), elapsed_seconds=record['elapsed_seconds'])), flush=True)


if __name__ == '__main__':
    main()
