"""Actual A51 ordered cubic Pi: paid33 cells plus missing928 genuine pairs."""
from pathlib import Path
import argparse
import json
import os
import resource
import sys
from time import perf_counter
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'0066'))
from column_common import CHANNELS, authenticate, bound_items, atomic_json, file_binding
from material_electron_ion_inputs import load_material
from general_moving_reference import SavedMovingReference
from cubic_contour import CubicContour, ordered_matrix
from consume_component import consume


def produce(output):
    started = perf_counter()
    out = Path(output)
    out.mkdir(exist_ok=False)
    spec = json.loads((HERE/'preregistration.json').read_bytes())
    group = Path('/proc/self/cgroup').read_text().strip().split('::')[-1]
    runtime = dict(affinity=sorted(os.sched_getaffinity(0)),
        threads={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')},
        limits={k:(Path('/sys/fs/cgroup')/group.lstrip('/')/k).read_text().strip() for k in ('memory.max','memory.swap.max')},
        executable=sys.executable,pythonpath=os.environ.get('PYTHONPATH'),
        herdr={k:v for k,v in os.environ.items() if k.startswith('HERDR_')})
    atomic_json(out/'actual-runtime.json',runtime)
    if (runtime['affinity'] != [9] or set(runtime['threads'].values()) != {'1'} or
            runtime['limits'] != {'memory.max':'3221225472','memory.swap.max':'0'} or
            runtime['herdr'].get('HERDR_ENV') != '1'):
        raise ValueError('actual Native CPU9/one-thread/3GiB/noSwap/owned view required')
    bundle = authenticate(spec['genuine_reference'])
    paid_result_path = authenticate(spec['paid_reuse']['first_column_result'])
    paid = json.loads(paid_result_path.read_bytes())
    column_path = authenticate(paid['output'])
    endpoint_path = authenticate(spec['paid_reuse']['endpoint_archive'])
    paid_inputs = json.loads(authenticate(paid['inputs']).read_bytes())
    source_path = authenticate(paid_inputs['compact_source'])
    c, manifest, coefficients = load_material()
    ref = SavedMovingReference(c,bundle,coefficients)
    if (not np.array_equal(ref.times,np.arange(31,dtype=float)) or
            paid_inputs['coefficient_manifest'] != coefficients or
            spec['genuine_reference'] not in list(bound_items(paid_inputs['reference']))):
        raise ValueError('unchanged original A52/material31 history required')
    provider = CubicContour.from_source(source_path)
    partner,left,right = provider.pair(0,1)
    if partner != 3 or left.shape != (111,111,111) or right.shape != left.shape:
        raise ValueError('independent complete original cubic vertices required')
    raw = out/'actual-inputs.npz'
    qs = sorted({1,partner,int(c.minus[1]),int(c.minus[partner])})
    factors = dict(left=left,right=right,times=ref.times,minus=c.minus,U=ref.u,rho0=c.rho0)
    for q in qs:
        factors['generator_'+str(q)] = c.mechanical.generators[q]
        factors['W0_'+str(q)] = c.mechanical.w0[q]
        factors['S_'+str(q)] = np.asarray([ref.S(q,t) for t in ref.times])
    np.savez(raw,**factors)
    input_shapes = {k:dict(shape=list(v.shape),dtype=str(v.dtype),payload_bytes=int(v.nbytes)) for k,v in factors.items()}
    del factors
    raw_archives = {entry['archive']['path']:entry['archive'] for entry in manifest['items'].values()}
    atomic_json(out/'inputs.json',dict(reference=ref.bindings,coefficient_manifest=coefficients,
        material_archives=list(raw_archives.values()),cubic_source=file_binding(source_path),
        cubic_source_arrays=provider.source.source_record['arrays'],
        cubic_source_bindings=provider.source.source_record['source_bindings'],
        actual_independent_vertices_and_raw_factors=dict(**file_binding(raw),arrays=input_shapes),
        paid_first_column_result=spec['paid_reuse']['first_column_result'],paid_first_column=file_binding(column_path),
        paid_endpoints=file_binding(endpoint_path)))
    paths = {k:out/(k+'.npy') for k in CHANNELS}
    arrays = {k:np.lib.format.open_memmap(paths[k],mode='w+',dtype=np.complex128,shape=(31,31,111,111)) for k in CHANNELS}
    coverage_path = out/'coverage.npy'
    seen = np.lib.format.open_memmap(coverage_path,mode='w+',dtype=bool,shape=(31,31))
    seen[:] = False
    fds = {k:os.open(paths[k],os.O_RDWR) for k in CHANNELS}
    coverage_fd = os.open(coverage_path,os.O_RDWR)
    def commit(a,b):
        # Dirty-page flush, never a full-bank copy; coverage follows numeric fsync.
        for key in CHANNELS:
            arrays[key].flush()
            os.fsync(fds[key])
        seen[a,b] = True
        seen.flush()
        os.fsync(coverage_fd)
    try:
        with np.load(column_path,allow_pickle=False) as saved:
            if (saved['greater'].shape != (31,111,111) or not np.array_equal(saved['times'],ref.times) or
                    int(saved['external_q']) != 0 or int(saved['internal_q']) != 1 or int(saved['partner_q']) != partner):
                raise ValueError('paid first column source axes differ')
            for a in range(31):
                for name in CHANNELS:
                    arrays[name][a,0] = saved[name][a]
                commit(a,0)
        with np.load(endpoint_path,allow_pickle=False) as saved:
            if (not np.array_equal(saved['times'],[0.,30.]) or int(saved['external_q']) != 0 or
                    int(saved['internal_q']) != 1 or int(saved['partner_q']) != partner):
                raise ValueError('original ordered endpoint route differs')
            for a,i in ((0,0),(30,1)):
                for name in CHANNELS:
                    arrays[name][a,30] = saved['old_'+name][i,1]
                commit(a,30)
        reused = np.argwhere(seen).tolist()
        tick = perf_counter()
        computed = 0
        for a in range(31):
            for b in range(1,31):
                if seen[a,b]:
                    continue
                dg1,dl1 = ref.boson(1,a,b)
                dg2,dl2 = ref.boson(partner,a,b)
                arrays['greater'][a,b] = ordered_matrix(left,right,dg1,dg2)
                arrays['lesser'][a,b] = ordered_matrix(left,right,dl1,dl2)
                arrays['retarded'][a,b] = arrays['greater'][a,b]-arrays['lesser'][a,b] if a>=b else 0.
                commit(a,b)
                computed += 1
            progress = dict(completed_owned_pairs=int(seen.sum()),required_owned_pairs=961,
                reused_owned_pairs=len(reused),new_owned_pairs=computed,elapsed_seconds=perf_counter()-started,
                new_evaluation_seconds=perf_counter()-tick,last_complete_row=a)
            atomic_json(out/'progress.json',progress)
            print('A72_COMMITTED '+json.dumps(progress),flush=True)
        numeric_seconds = perf_counter()-tick
        if not seen.all() or computed != 928 or len(reused) != 33:
            raise ValueError('actual complete961 with only missing928 required')
        record = dict(schema='p254-additive-history-component-v1',status='complete-owned-component-q1',
            component='cubic',input_dataset=spec['genuine_reference'],times_inverse_ry=ref.times.tolist(),
            owned_q=[1],completed_q=[1],external_q=0,internal_q=1,partner_q=partner,
            full_mesh=False,full_time_square=True,axes=spec['axes'],shape=spec['shape'],
            arrays={k:dict(**file_binding(paths[k]),shape=list(arrays[k].shape),dtype=str(arrays[k].dtype)) for k in CHANNELS},
            coverage=dict(array=file_binding(coverage_path),owned_pairs='all',completed_pairs=np.argwhere(seen).tolist(),count=int(seen.sum())),
            paid_reuse=dict(column=file_binding(column_path),endpoints=file_binding(endpoint_path),reused_pairs=reused),
            newly_evaluated_pairs=computed,input_bindings=file_binding(out/'inputs.json'),
            implementation_bindings=[file_binding(module.__file__) for module in list(sys.modules.values())
                if getattr(module,'__file__',None) and str(HERE.parent) in str(module.__file__) and str(module.__file__).endswith('.py')],
            preregistration=file_binding(HERE/'preregistration.json'),raw_nonstationary=True,theta_diagonal=1,
            equation='A62lines126-144/A51 ordered_matrix; independently assembled left/right and greater/lesser opposite-q transpose legs',
            no_claims=['all64q','complete six-array functional','field-cubic','selfconsistency','physical ranking','rates','heat'],
            new_evaluation_seconds=numeric_seconds,elapsed_seconds_before_consumer=perf_counter()-started,
            peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        catalog = out/'component.json'
        atomic_json(catalog,record)
        consumer = consume(catalog)
        atomic_json(out/'result.json',dict(status=record['status'],catalog=file_binding(catalog),consumer=file_binding(out/'consumer.json'),
            saved_consumer=consumer,elapsed_seconds=perf_counter()-started,
            new_evaluation_seconds=numeric_seconds,completed_pairs=int(seen.sum()),reused_pairs=33,new_pairs=computed,
            peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            bounded_incomplete_boundary_loss='None on this complete exit; interruption before coverage commit can lose at most one numeric pair, never marks it complete.'))
        print('A72_NATIVE_COMPLETE '+json.dumps(file_binding(out/'result.json')),flush=True)
    finally:
        for fd in fds.values():
            os.close(fd)
        os.close(coverage_fd)
        ref.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    produce(parser.parse_args().output)
