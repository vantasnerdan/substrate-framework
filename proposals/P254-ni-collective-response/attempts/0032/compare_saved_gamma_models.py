"""Compare original and moving-hybrid Gamma readouts without another trajectory.

The old preparation-time dual is the executed A30 reconstruction of A29's
model, not recovered original adaptive history. Difference is a declared-model
enrichment contrast, NOT an all-source passive-noise error estimate.
"""
import argparse
import json
from pathlib import Path
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'),str(HERE.parent/'0029'),str(HERE.parent/'0030')]

import numpy as np
from clamped_scientific_readout import encode
from export_output_directed_actions import restore_filter
from reduced_correlation_transpose import wick_edge_pair
from reduced_coherence_readout import normalized_ordered
from substrate_framework.material_clamped_array import ClampedArray
from substrate_framework.streamed_material_state import authenticate, file_binding


def arrays(binding):
    with np.load(authenticate(binding),allow_pickle=False) as archive:
        return {key:archive[key] for key in archive.files}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hybrid-result',required=True)
    parser.add_argument('--filter-packet',required=True)
    parser.add_argument('--output',required=True)
    args = parser.parse_args()
    started = perf_counter()
    hybrid_binding = file_binding(args.hybrid_result)
    hybrid = json.loads(Path(args.hybrid_result).read_bytes())
    filter_,filter_binding,filter_archive = restore_filter(args.filter_packet)
    packet = json.loads(Path(args.filter_packet).read_bytes())
    if (hybrid['status'] != 'complete' or hybrid['q_index'] != 0
            or hybrid['source']['original_result']['sha256'] != packet['source']['original_result']['sha256']
            or hybrid['source']['coordinate'] != packet['filter']['coordinates'][0]
            or hybrid['source']['q_index'] != packet['filter']['q_index']):
        raise ValueError('complete same-original-input Gamma/Ni-coordinate models required')
    source = hybrid['source']
    bindings = source['coefficient_source_archives']
    ns = hybrid['model_metadata']['source_count']
    m = hybrid['model_metadata']['active_dimension']
    mechanical = ClampedArray(**arrays(bindings['mechanical']))
    rho = np.stack([arrays(bindings['slot/%03d' % slot])['rho0'] for slot in range(ns)])
    nx = hybrid['model_metadata']['mechanical_size']
    old_terminal = filter_.initial_dual
    if old_terminal.shape != (nx+ns*m*m,1):
        raise ValueError('original coupled dual layout differs')
    old_initial = np.zeros_like(old_terminal)
    old_initial[source['coordinate'],0] = 1
    def old_ordered(left,right):
        ionic = left[:nx].T@mechanical.w0[0]@right[:nx]
        electronic = wick_edge_pair(left[nx:].reshape(ns,m,m,-1),
            right[nx:].reshape(ns,m,m,-1),rho,rho)
        return (ionic+electronic).item()
    old_00 = old_ordered(old_initial,old_initial)
    old_tt = old_ordered(old_terminal,old_terminal)
    old_t0 = old_ordered(old_terminal,old_initial)
    saved = arrays(hybrid['artifacts']['readout'])
    times = saved['times_inverse_ry']
    index = np.flatnonzero(times == filter_.final_time)
    zero = np.flatnonzero(times == 0)
    if len(index) != 1 or len(zero) != 1:
        raise ValueError('actual common requested Gamma endpoints required')
    t,z = int(index[0]),int(zero[0])
    covariance = saved['ordered_covariance']
    old_normalized = normalized_ordered(old_t0,old_tt,old_00)
    new_normalized = normalized_ordered(covariance[t,z],covariance[t,t],covariance[z,z])
    rows = [dict(model='original_active_reconstructed_Gamma',equal_time_initial=old_00,
        equal_time_terminal=old_tt,terminal_initial_ordered=old_t0,normalized=old_normalized),
        dict(model='moving_source0_enriched_Gamma',equal_time_initial=covariance[z,z],
        equal_time_terminal=covariance[t,t],terminal_initial_ordered=covariance[t,z],normalized=new_normalized)]
    result = dict(schema='P254-A32-saved-Gamma-model-contrast-v1',status='complete',
        hybrid_result=hybrid_binding,original_filter_packet=filter_binding,
        original_filter_archive=filter_archive,hybrid_readout=hybrid['artifacts']['readout'],
        original_result=source['original_result'],time_inverse_ry=float(times[t]),
        coordinate=source['coordinate'],atom=1,axis='x',q_index=0,
        rows=rows,raw_initial_variance_difference=covariance[z,z]-old_00,
        terminal_variance_difference=covariance[t,t]-old_tt,
        terminal_initial_ordered_difference=covariance[t,z]-old_t0,
        normalized_complex_difference=(new_normalized['complex_degree']-old_normalized['complex_degree']
            if new_normalized['defined'] and old_normalized['defined'] else None),
        coefficient_source_archives=bindings,
        implementation=[file_binding(__file__)],elapsed_seconds=perf_counter()-started,
        scope='Same original preparation and Ni1/x Gamma observable;old executed reconstruction vs movinghybrid source0finitebath enrichment. Difference includes changed source0literalH and coupled dynamics/noise,NOT isolated passive noise. No spatial Fourier sum,allsource/fullPWerrorbar,irreversibleloss,lifetime,materialrankingorcampaignclosure.')
    with Path(args.output).open('x') as stream:
        json.dump(encode(result),stream,indent=2,sort_keys=True,allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(result=file_binding(args.output),elapsed_seconds=result['elapsed_seconds']),allow_nan=False))


if __name__ == '__main__':
    main()
