"""Paid q1 Gamma cubic first column and genuine bare source factors0..30."""
from pathlib import Path
import argparse
import json
import resource
from time import perf_counter
import sys
from column_common import (HERE, BUNDLE, CHANNELS, begin, authenticate,
                           bound_items, atomic_json, file_binding)
import numpy as np
from material_electron_ion_inputs import load_material
from general_moving_reference import SavedMovingReference
from cubic_contour import CubicContour, ordered_matrix


def produce(output):
    started = perf_counter()
    out, spec = begin(output, 8)
    endpoint = authenticate(spec['saved_endpoints_cubic'])
    paid = json.loads(Path('/var/tmp/p254-a65-time-factor-001/result.json').read_bytes())
    c, manifest, coefficient_binding = load_material()
    reference = SavedMovingReference(c, BUNDLE, coefficient_binding)
    if (not np.array_equal(reference.times, np.arange(31, dtype=float)) or
            file_binding(BUNDLE) not in list(bound_items(paid['input_bindings'])) or
            coefficient_binding not in list(bound_items(paid['input_bindings']))):
        raise ValueError('same actual source/material and genuine31 nodes required')
    source = authenticate(spec['original_source'])
    provider = CubicContour.from_source(source)
    partner, left, right = provider.pair(0,1)
    arrays = {name:np.empty((31,111,111), dtype=complex) for name in CHANNELS}
    seen = np.zeros(31, dtype=bool)
    with np.load(endpoint, allow_pickle=False) as archive:
        if (int(archive['external_q']) != 0 or int(archive['internal_q']) != 1 or
                int(archive['partner_q']) != partner or
                not np.array_equal(archive['times'], [0.,30.])):
            raise ValueError('paid original cubic route and genuine endpoints differ')
        for name in CHANNELS:
            old = archive['old_'+name]
            if old.shape != (2,2,111,111):
                raise ValueError('paid full111 endpoint axes differ')
            arrays[name][0], arrays[name][30] = old[0,0], old[1,0]
    seen[[0,30]] = True
    atomic_json(out/'actual-inputs.json', dict(coefficient_manifest=coefficient_binding,
        reference=reference.bindings, compact_source=file_binding(source),
        paid_endpoint=file_binding(endpoint),external_q=0,internal_q=1,partner_q=partner))
    tick = perf_counter()
    for a in range(1,30):
        dg1, dl1 = reference.boson(1,a,0)
        dg2, dl2 = reference.boson(partner,a,0)
        arrays['greater'][a] = ordered_matrix(left,right,dg1,dg2)
        arrays['lesser'][a] = ordered_matrix(left,right,dl1,dl2)
        arrays['retarded'][a] = arrays['greater'][a]-arrays['lesser'][a]
        seen[a] = True
        # Every newly paid pair is recoverable, not only an elapsed-time marker.
        np.savez(out/('paid-pair-%02d.npz'%a),ordered_pair=np.array([a,0]),
            **{name:arrays[name][a] for name in CHANNELS})
        atomic_json(out/'progress.json',dict(completed_owned_pairs=int(seen.sum()),
            required_owned_pairs=31,elapsed_seconds=perf_counter()-started))
    seconds = perf_counter()-tick
    if not seen.all() or not all(np.isfinite(value).all() for value in arrays.values()):
        raise ValueError('complete finite owned first cubic column absent')
    column = out/'column.npz'
    np.savez(column,times=reference.times,ordered_pairs=np.array([(a,0) for a in range(31)]),
        external_q=np.array(0),internal_q=np.array(1),partner_q=np.array(partner),**arrays)
    factors = out/'bare-prescribed-source-factors.npz'
    np.savez(factors,times=reference.times,U=reference.u,rho0=c.rho0,
        S_gamma=np.asarray([reference.S(0,t) for t in reference.times]),
        W0_gamma=c.mechanical.w0[0],n=np.array(c.n))
    with np.load(column,allow_pickle=False) as saved:
        probe = np.eye(111,dtype=complex)[:,0]
        consumer = dict(final_retarded_ni1_action_norm=float(np.linalg.norm(saved['retarded'][-1]@probe)))
    with np.load(factors,allow_pickle=False) as saved:
        if (not np.array_equal(saved['U'][0],np.broadcast_to(np.eye(12),(128,12,12))) or
                not np.array_equal(saved['S_gamma'][0],np.eye(228))):
            raise ValueError('literal original source preparation differs')
    record = dict(schema='p254-paid-cubic-long-history-column-v1',status='complete-first-column-q1',
        output=file_binding(column),bare_source_factors=file_binding(factors),
        inputs=file_binding(out/'actual-inputs.json'),times_inverse_ry=reference.times.tolist(),
        external_q=0,internal_q=1,partner_q=partner,coordinates=111,
        reused_time_pairs=[[0,0],[30,0]],new_time_pairs=[[a,0] for a in range(1,30)],
        completed_owned_pairs=int(seen.sum()),new_evaluation_seconds=seconds,
        elapsed_seconds=perf_counter()-started,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        consumer=consumer,full_time_square=False,full_internal_mesh=False,
        reference_source_kind='bare-mechanical-with-contacts',
        source_factor_semantics='Actual saved U31 and original Gaussian S/W0/rho0; no common mean-feedback trajectory or represented-control history is claimed.',
        source_bindings=[file_binding(module.__file__) for module in list(sys.modules.values())
            if getattr(module,'__file__',None) and str(HERE.parent) in str(module.__file__)
            and str(module.__file__).endswith('.py')],
        limits=['First column/q1 only; not full intrinsic evolution, conservation or physical ranking.'])
    atomic_json(out/'result.json',record)
    reference.close()
    print('A66_CUBIC_NATIVE_COMPLETE '+json.dumps(file_binding(out/'result.json')),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    produce(parser.parse_args().output)
