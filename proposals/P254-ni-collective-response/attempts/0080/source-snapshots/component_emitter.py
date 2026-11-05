"""Exact original-mesh additive histories; numerical API never launches jobs.

produce_component defaults to all genuine 31 saved nodes and all 961 pairs.
Explicit pairs index the selected time axis, not the original saved-node axis.
An output directory is single-use: no implicit retry, refill or overwrite.
"""
from pathlib import Path
import hashlib
import json
import os
import resource
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent/n) for n in (
    '0064','0062','0063','0060','0056','0055','0054','0053','0052','0051',
    '0050','0048','0046','0045','0044','0043','0042','0041','0030','0029',
    '0028','0027','0025','0024')]
import numpy as np
from history_banks import ORDERS, binding, authenticate, save
from general_moving_reference import SavedMovingReference
from material_electron_ion_inputs import load_material, edge_pair
from electron_ion_functional import evaluate_components
from coulomb_material_inputs import load_coulomb_inputs
from bare_coulomb_factors import BareCoulombFactors
from gaussian_ee import stream_gaussian_ee
from cubic_contour import CubicContour, ordered_matrix

ORIGINAL_BUNDLE = Path('/var/tmp/p254-a52-requested-reference-001/producer/bundle.json')
CUBIC_SOURCE = Path('/var/tmp/p254-a47-full37-source-001/source.json')


def numeric_identity(value):
    value = np.asarray(value)
    digest = hashlib.sha256()
    for row in value:
        digest.update(np.ascontiguousarray(row).tobytes())
    return dict(shape=list(value.shape), dtype=str(value.dtype), bytes=value.nbytes,
                sha256=digest.hexdigest(), hash_scope='C-order numeric payload')


def produce_component(component, q, output, *, bundle=ORIGINAL_BUNDLE,
                      times=None, owned_pairs=None):
    """Persist one actual outer/internal-q contribution using owning primitives.

    EE retains every 1728 exchange mode and all 64 inner h. Mixed EI contains
    sigma only, no Gamma joint term. Cubic is external Gamma, internal q with
    its independently routed partner. All components retain original axes.
    Completion is restricted to declared pair/slot ownership; never full mesh.
    """
    if component not in ('ee', 'mixed_ei', 'cubic'):
        raise ValueError('component must be ee, mixed_ei or cubic')
    if type(q) is not int or not 0 <= q < 64:
        raise ValueError('q must be an actual original integer 0..63')
    prereg = json.loads((HERE/'preregistration.json').read_bytes())
    if binding(bundle) != prereg['reference']:
        raise ValueError('unchanged original A52 input identity required')
    authenticate(prereg['owning_combination'])
    authenticate(prereg['EE_kernel'])
    selected = np.arange(31, dtype=float) if times is None else np.asarray(times, float)
    if (selected.ndim != 1 or not len(selected) or not np.isfinite(selected).all()
            or np.any(np.diff(selected) <= 0)):
        raise ValueError('strictly increasing nonempty exact saved times required')
    nt = len(selected)
    pairs = tuple((a,b) for a in range(nt) for b in range(nt)) if owned_pairs is None else tuple(tuple(p) for p in owned_pairs)
    if (not pairs or len(set(pairs)) != len(pairs) or any(len(p) != 2 or
            any(type(x) is not int or not 0 <= x < nt for x in p) for p in pairs)):
        raise ValueError('unique nonempty local-index owned pairs required')
    out = Path(output).resolve()
    out.mkdir(exist_ok=False)
    began = time.monotonic()
    ref = None
    arrays, fds = {}, {}
    seen = None
    record = None
    failure = None
    try:
        if component == 'ee':
            c = load_coulomb_inputs()
            coefficient = c.bindings['manifest']
            inputs = dict(material=c.bindings)
            manifest = None
        else:
            c, manifest, coefficient = load_material()
            inputs = dict(coefficient_manifest=coefficient)
        ref = SavedMovingReference(c, bundle, coefficient, times=selected)
        source_times = np.asarray(ref.spec['times_inverse_ry'])
        original_indices = [int(np.flatnonzero(source_times == t)[0]) for t in ref.times]
        targets = [int(c.target[q,s]) for s in range(c.ns)]
        if sorted(targets) != list(range(128)):
            raise ValueError('original source-target routing is not a permutation')
        routes = dict(source=list(range(128)), target=targets,
            slot_spin=np.asarray(c.slot_spin).tolist(), slot_k=np.asarray(c.slot_k).tolist(),
            opposite_q=int(c.minus[q]), outer_q=q)
        inputs.update(reference=ref.bindings, rho0=numeric_identity(c.rho0))
        data = provider = left = right = None
        if component == 'mixed_ei':
            data = c.q(q)
            keys = ['mechanical', 'kernel/%02d'%q] + ['slot/%03d'%s for s in range(128)] + ['edge/%02d/%03d'%(q,s) for s in range(128)]
            archives = {manifest['items'][key]['archive']['path']:manifest['items'][key]['archive'] for key in keys}
            inputs.update(actual_consumed_item_keys=keys, authenticated_archives=list(archives.values()),
                ionic=numeric_identity(data.ionic), ionic_reverse=numeric_identity(data.ionic_reverse))
        if component == 'cubic':
            provider = CubicContour.from_source(CUBIC_SOURCE)
            partner, left, right = provider.pair(0,q)
            if left.shape != (111,111,111) or right.shape != left.shape:
                raise ValueError('full independent original cubic vertices required')
            routes.update(external_q=0, internal_q=q, partner_q=int(partner))
            inputs.update(cubic_source=binding(CUBIC_SOURCE),
                cubic_source_arrays=provider.source.source_record['arrays'],
                cubic_source_bindings=provider.source.source_record['source_bindings'],
                left=numeric_identity(left), right=numeric_identity(right))
        if component != 'ee':
            mechanical_qs = {q, int(c.minus[q])}
            if component == 'cubic':
                mechanical_qs.update((int(partner), int(c.minus[partner])))
            inputs['mechanical'] = {str(k):dict(generator=numeric_identity(c.mechanical.generators[k]),
                W0=numeric_identity(c.mechanical.w0[k])) for k in sorted(mechanical_qs)}
        del manifest
        save(out/'inputs.json',inputs)
        nslots = 1 if component == 'cubic' else 128
        shape = (nt,nt,111,111) if component == 'cubic' else (nt,nt,128,12,12)
        prefix = '' if component == 'cubic' else 'sigma-'
        for order in ORDERS:
            path = out/(prefix+order+'.npy')
            arrays[order] = np.lib.format.open_memmap(path,mode='w+',dtype=np.complex128,shape=shape)
            arrays[order][:] = 0
            fds[order] = os.open(path,os.O_RDWR)
        seen_path = out/'slot-coverage.npy'
        seen = np.lib.format.open_memmap(seen_path,mode='w+',dtype=bool,shape=(nt,nt,nslots))
        seen[:] = False
        fds['coverage'] = os.open(seen_path,os.O_RDWR)
        record = dict(schema='p254-additive-history-component-v1',component=component,
            status='partial-owned-component-q%d'%q,full_mesh=False,full_solve_eligible=False,
            input_dataset=binding(bundle), times_inverse_ry=ref.times.tolist(),
            selected_original_time_indices=original_indices,owned_q=[q],completed_q=[],
            axes=['time','time','ionic_coordinate','ionic_coordinate'] if component == 'cubic' else ['time','time','slot','band','band'],
            shape=list(shape),slots=nslots,bands=12,spins=2,ionic_directions=111,
            canonical_routes=routes, input_bindings=binding(out/'inputs.json'),
            source_bindings=ref.bindings,preregistration=binding(HERE/'preregistration.json'),
            coverage=dict(owned_pairs=[list(p) for p in pairs],completed_pairs=[],
                external_targets=list(range(nslots)), all111_ionic_directions=True,
                commit='three numeric slices flush/fsync, then slot bitmap flush/fsync, then atomic checkpoint'),
            full_time_square=False,raw_nonstationary=True,theta_diagonal=1,
            scope='Bare original Gaussian additive component only; no joint, all64q, closure, physical ranking or full solve',
            timings=[])
        if component == 'cubic':
            record.update(external_q=0,internal_q=q,partner_q=int(partner))
        save(out/'checkpoint.json',record)

        def commit(completed, slots):
            for order in ORDERS:
                arrays[order].flush()
                os.fsync(fds[order])
            for a,b in completed:
                if np.any(seen[a,b,list(slots)]):
                    raise ValueError('duplicate numeric ownership commit')
                seen[a,b,list(slots)] = True
            seen.flush()
            os.fsync(fds['coverage'])
            record['coverage']['completed_pairs'] = [list(p) for p in pairs if seen[p].all()]
            record['elapsed_seconds'] = time.monotonic()-began
            save(out/'checkpoint.json',record)

        if component == 'ee':
            def consume(block):
                if block.q != q:
                    raise ValueError('unexpected q from owning EE stream')
                for j,(a,b) in enumerate(block.pairs):
                    for order in ORDERS:
                        arrays[order][a,b,list(block.external)] = block.values[order+'_total'][:,j]
                commit(block.pairs,block.external)
            record['ee_stream'] = stream_gaussian_ee(BareCoulombFactors(c),ref,ref.times,consume,
                owned_q=(q,),owned_pairs=pairs,inner_batch=8,internal_batch=8,pair_batch=4,
                source_bindings=ref.bindings)
        else:
            for a,b in pairs:
                tick = time.monotonic()
                if component == 'mixed_ei':
                    forward,reverse = ref.electron(a,b),ref.electron(b,a)
                    boson = ref.boson(q,a,b)
                    for source_slot in range(c.ns):
                        target = int(c.target[q,source_slot])
                        values = evaluate_components(edge_pair(data,source_slot,target,forward,reverse,boson),
                            n=c.n,t=ref.times[a],s=ref.times[b],components=('sigma',))
                        for order in ORDERS:
                            arrays[order][a,b,target] = values['sigma_'+order]
                else:
                    first,second = ref.boson(q,a,b),ref.boson(partner,a,b)
                    arrays['greater'][a,b] = ordered_matrix(left,right,first[0],second[0])
                    arrays['lesser'][a,b] = ordered_matrix(left,right,first[1],second[1])
                    arrays['retarded'][a,b] = arrays['greater'][a,b]-arrays['lesser'][a,b] if ref.times[a]>=ref.times[b] else 0.
                commit(((a,b),),tuple(range(nslots)))
                record['timings'].append(dict(pair=[a,b],seconds=time.monotonic()-tick))
        if len(record['coverage']['completed_pairs']) != len(pairs):
            raise RuntimeError('owning primitive did not complete declared coverage')
        record['status'] = 'complete-owned-component-q%d'%q
        record['completed_q'] = [q]
        record['full_time_square'] = np.array_equal(ref.times,np.arange(31)) and len(pairs)==961
    except BaseException as exc:
        import traceback
        failure = dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
        raise
    finally:
        if record is not None:
            record['failure'] = failure
            record['elapsed_seconds'] = time.monotonic()-began
            record['native_cpu_seconds'] = time.process_time()
            record['peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
            record['u_integration_calls'] = ref.u_integration_calls
            record['gaussian_expm_calls'] = ref.expm_calls
            record['implementation_bindings'] = [binding(m.__file__) for m in list(sys.modules.values())
                if getattr(m,'__file__',None) and str(HERE.parent) in str(m.__file__) and str(m.__file__).endswith('.py')]
            descriptors = {}
            for order,arr in arrays.items():
                arr.flush()
                os.fsync(fds[order])
                descriptors[order] = dict(binding(out/(prefix+order+'.npy')),shape=list(arr.shape),dtype=str(arr.dtype))
            record['arrays' if component == 'cubic' else 'sigma'] = descriptors
            seen.flush()
            os.fsync(fds['coverage'])
            record['coverage']['array'] = dict(binding(seen_path),shape=list(seen.shape),dtype=str(seen.dtype))
            record['coverage']['completed_slot_cells'] = int(seen.sum())
            save(out/'catalog.json',record)
        for fd in fds.values():
            os.close(fd)
        for arr in arrays.values():
            arr._mmap.close()
        if seen is not None:
            seen._mmap.close()
        if ref is not None:
            ref.close()
    return out/'catalog.json'
