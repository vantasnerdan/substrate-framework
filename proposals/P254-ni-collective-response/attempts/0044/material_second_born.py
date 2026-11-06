"""Complete momentum/internal-spin sums, using separable bare material factors."""
import numpy as np
from second_born_functional import terms_from_vertices, evaluate, jvp, transpose


def local_legs(forward, backward, a, b, c):
    gg,gl = forward
    bg,bl = backward
    return dict(greater1=gg[b],greater2=gg[c],lesser1=gl[b],lesser2=gl[c],
                reverse_greater=bg[a],reverse_lesser=bl[a])


def empty_result(m):
    return {order:{channel:np.zeros((m,m),complex) for channel in ('direct','exchange','total')}
            for order in ('greater','lesser','retarded')}


def _add(result, value):
    for order in result:
        for channel in result[order]:
            result[order][channel] += value[order][channel]


def evaluate_material(factors, external, forward, backward, *, retarded=True):
    """All original internal slots/bands/features for ONE declared external slot.

    forward=(Ggreater(t,s),Glesser(t,s)), backward independent reversed pair;
    each array [128,12,12]. No momentum, spin or time-reversal reduction.
    """
    result = empty_result(factors.c.m)
    quartets = terms = 0
    for a,b,c in factors.quartets(external):
        block = terms_from_vertices(*factors.term_factors(external,a,b,c))
        _add(result,evaluate(block,local_legs(forward,backward,a,b,c),retarded=retarded))
        quartets += 1; terms += len(block)
    return result,dict(momentum_spin_quartets=quartets,expanded_vertex_terms=terms,
                       internal_band_sums='all12 at each of six indices',ncells=factors.c.n)


def jvp_material(factors, external, forward, backward, forward_direction,
                 backward_direction, *, factor_direction=None, retarded=True):
    """Common functional JVP. Optional callback receives (i,a,b,c,terms)."""
    result = empty_result(factors.c.m)
    for a,b,c in factors.quartets(external):
        block = terms_from_vertices(*factors.term_factors(external,a,b,c))
        df = None if factor_direction is None else factor_direction(external,a,b,c,block)
        value = jvp(block,local_legs(forward,backward,a,b,c),
            local_legs(forward_direction,backward_direction,a,b,c),
            factor_direction=df,retarded=retarded)
        _add(result,value)
    return result


def transpose_material(factors, external, forward, backward, dual, *,
                       factor_scatter=None, retarded=True):
    """Scatter every repeated source/target and reversed-time G occurrence.

    Factor cotangents can be consumed immediately by factor_scatter(i,a,b,c,
    terms,cotangents). They are independent factor-coordinate bilinear duals,
    not a claimed complex-linear derivative through D.conj(). Physical D is a
    real-linear coordinate: pullback through conjugation conjugates its dual
    only when the consumer pairs Re(sum(z*dD)). No Hilbert transpose is hidden.
    """
    gf = [np.zeros_like(x) for x in forward]
    gb = [np.zeros_like(x) for x in backward]
    for a,b,c in factors.quartets(external):
        block = terms_from_vertices(*factors.term_factors(external,a,b,c))
        value = transpose(block,local_legs(forward,backward,a,b,c),dual,retarded=retarded)
        for name,target in (('greater1',gf[0][b]),('greater2',gf[0][c]),
                            ('lesser1',gf[1][b]),('lesser2',gf[1][c]),
                            ('reverse_greater',gb[0][a]),('reverse_lesser',gb[1][a])):
            target += value['legs'][name]
        if factor_scatter is not None:
            factor_scatter(external,a,b,c,block,value['factors'])
    return tuple(gf),tuple(gb)


def evaluate_material_time_grid(factors, external, electrons, times, cache):
    """Outer-transfer-q chunks, quartet-outer all-time contraction within each q.

    Static factors are constructed once per quartet. Original quartet order is
    retained inside each transfer q; final sums follow ascending q. Only the
    current uncommitted q can be lost on interruption. ProjectionCache owns
    immutable archive authentication/publication; no original-bank copy.
    """
    from time import perf_counter
    nt,m = len(times),factors.c.m
    keys = [o+'_'+k for o in ('greater','lesser','retarded')
            for k in ('direct','exchange','total')]
    total = {k:np.zeros((nt,nt,m,m),complex) for k in keys}
    groups = [[] for _ in range(factors.c.nq)]
    transfers = [factors.transfer(external,b) for b in range(factors.c.ns)]
    for a,b,c in factors.quartets(external):
        groups[transfers[b]].append((a,b,c))
    counts = dict(quartets=0,new_quartets=0,reused_quartets=0,
                  new_q_chunks=0,reused_q_chunks=0,expanded_vertex_terms=0,q_chunks=[])
    for q,quartets in enumerate(groups):
        key = f'slot/{external:03d}/q/{q:02d}'
        reused = cache.has(key)
        if not reused:
            started = perf_counter()
            arrays = {k:np.zeros((nt,nt,m,m),complex) for k in keys}
            term_count = 0
            for a,b,c in quartets:
                block = terms_from_vertices(*factors.term_factors(external,a,b,c))
                term_count += len(block)
                for ti,t in enumerate(times):
                    for si,s in enumerate(times):
                        value = evaluate(block,local_legs(electrons[ti,si],electrons[si,ti],a,b,c),
                                         retarded=t>=s)
                        for order,channels in value.items():
                            for channel,entry in channels.items():
                                arrays[order+'_'+channel][ti,si] += entry
            if not all(np.isfinite(x).all() for x in arrays.values()):
                raise FloatingPointError('nonfinite raw transfer-q chunk, no repair')
            cache.put(key,arrays,metadata=dict(external=external,outer_transfer_q=q,
                quartets=len(quartets),quartet_order=[list(x) for x in quartets],
                expanded_vertex_terms=term_count),timings=dict(seconds=perf_counter()-started))
        # The final consumer reads authenticated committed q chunks, never RAM echoes.
        arrays = cache.get(key)
        record = cache.manifest['items'][key]
        metadata = record['metadata']
        if (metadata['external'] != external or metadata['outer_transfer_q'] != q
                or metadata['quartet_order'] != [list(x) for x in quartets]
                or metadata['quartets'] != len(quartets)
                or set(arrays) != set(keys) or any(x.shape != (nt,nt,m,m) for x in arrays.values())):
            raise ValueError('cached all-time transfer-q identity/layout differs')
        for name in keys:
            total[name] += arrays[name]
        counts['quartets'] += metadata['quartets']
        counts['reused_quartets' if reused else 'new_quartets'] += metadata['quartets']
        counts['reused_q_chunks' if reused else 'new_q_chunks'] += 1
        counts['expanded_vertex_terms'] += metadata['expanded_vertex_terms']
        counts['q_chunks'].append(dict(q=q,archive=record['archive'],metadata=metadata,
                                      timings=cache.manifest['timings'][key],reused=reused))
        print(dict(external=external,q=q,quartets=metadata['quartets'],
                   expanded_vertex_terms=metadata['expanded_vertex_terms'],reused=reused),flush=True)
    return total,counts


def source_bounds(external_count, pair_count):
    """Conservative combinatorial/storage bounds, not measured timings."""
    return dict(external_slots=external_count,time_pairs=pair_count,
        quartets_per_external_pair=12288,expanded_terms_per_external_pair_upper=32768,
        ordered_component_term_calls_upper=65536*external_count*pair_count,
        band_width=12,feature_width_upper=46,
        density_bytes_upper=64*128*46*12*12*16,
        static_factor_constructions_upper=12288*external_count,
        q_output_archives=64*external_count,
        q_payload_bytes_upper=64*external_count*pair_count*9*12*12*16,
        final_payload_bytes=external_count*pair_count*9*12*12*16,
        uncommitted_q_quartet_loss_upper=192,
        uncommitted_q_ordered_term_calls_upper=1024*pair_count,
        staging_payload_bytes_upper=pair_count*9*12*12*16,
        q_output_storage_policy='one immutable all-time9-component tensor per outer transfer-q/external; add NPZ/JSON/filesystem overhead; no old-bank copy',
        time_batch_resident_bytes=pair_count*9*12*12*16,
        einsum_intermediate_elements_cap=2_000_000,
        einsum_intermediate_bytes_upper=32_000_000,
        global_four_index_allocated=False,persistent_original_bank_copy_bytes=0,
        walltime_bound=None,walltime_policy='Root must grant actual source-bound run; no measured cost claim',
        arithmetic_bound='per expanded term, two ordered contractions, each at most46^2*12^8 naive scalar products; cached greedy paths reduce this but no timing inferred')
