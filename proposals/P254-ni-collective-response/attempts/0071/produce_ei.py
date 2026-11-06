"""A62 lines100–115 only: actual q1 raw EI Sigma, complete-pair durable checkpoints."""
from pathlib import Path
import hashlib
import json
import os
import resource
import signal
import sys
import time
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent/n) for n in ('0062','0063','0060','0056','0055','0054','0053','0052','0051','0050','0048','0046','0045','0044','0043','0042','0041','0030','0029','0028','0027','0025','0024')]
import numpy as np
from history_banks import ORDERS, binding, authenticate, save
from history_sources import open_reference
from joint_material_inputs import load_inputs
from material_electron_ion_inputs import edge_pair
from electron_ion_functional import evaluate_components


def array_identity(value):
    value = np.asarray(value)
    digest = hashlib.sha256()
    if value.flags.c_contiguous:
        digest.update(memoryview(value).cast('B'))
    else:
        for row in value:
            digest.update(np.ascontiguousarray(row).tobytes())
    return dict(shape=list(value.shape), dtype=str(value.dtype), bytes=value.nbytes,
                sha256=digest.hexdigest(), hash_scope='C-order raw numeric payload')


def produce(output, deadline):
    started = time.monotonic()
    prereg = json.loads((HERE/'preregistration.json').read_bytes())
    source = authenticate(prereg['genuine_reference'])
    authenticate(prereg['owning_implementation'])
    output = Path(output)
    output.mkdir(exist_ok=False)
    stop = []
    signal.signal(signal.SIGTERM, lambda number,frame: stop.append('SIGTERM'))
    c, baseline, old, material_bindings = load_inputs()
    old.close()
    del old, baseline
    ref = open_reference(c, material_bindings['manifest'], bundle=source)
    banks = {}
    handles = {}
    record = dict(schema='p254-additive-history-component-v1', component='mixed_ei',
        status='partial-owned-component-q1', full_mesh=False, full_solve_eligible=False,
        input_dataset=binding(source), times_inverse_ry=ref.times.tolist(), owned_q=[1], completed_q=[],
        axes=['time','time','slot','band','band'], slots=128, bands=12, spins=2,
        ionic_directions=111, normalization=dict(original_N=c.n, factor='1/N once in evaluate_components'),
        source_bindings=ref.bindings, coefficient_bindings=material_bindings,
        implementation_bindings=[binding(HERE/name) for name in ('produce_ei.py','component_reader.py')],
        owning_implementation=prereg['owning_implementation'],
        component_implementation_bindings=[binding(Path(__import__(name).__file__)) for name in
            ('electron_ion_functional','material_electron_ion_inputs','joint_material_inputs','history_sources','general_moving_reference')],
        coverage=dict(owned_pairs=[[a,b] for a in range(31) for b in range(31)], completed_pairs=[],
                      external_targets=list(range(128)), all111_ionic_directions=True,
                      checkpoint='all three numeric pair slices flushed and fsynced before coverage commit'),
        paid_node_reuse=[], preparation=ref.spec['initial_conditions'],
        scope='Raw q1 EI Sigma only. No EE/cubic/joint Gamma correction; no full all64q claim.',
        source_identity=dict(U=dict(ref.spec['U'],shape=list(ref.u.shape),dtype=str(ref.u.dtype))),
        timings=[], boundary_loss=dict(complete_pairs_lost=0, incomplete_pair=None))
    failure = None
    active_pair = None
    try:
        if not np.array_equal(ref.times,np.arange(31)) or (c.ns,c.m,c.d,c.n)!=(128,12,111,64):
            raise ValueError('frozen full31/original physical axes differ')
        data = c.q(1)
        targets = [int(c.target[1,s]) for s in range(c.ns)]
        if sorted(targets) != list(range(128)):
            raise ValueError('q1 source-target routes do not cover each canonical target once')
        manifest = json.loads(Path(material_bindings['manifest']['path']).read_bytes())
        keys = ['mechanical','kernel/01']+['slot/%03d'%s for s in range(128)]+['edge/01/%03d'%s for s in range(128)]
        record['actual_consumed_items'] = {key:manifest['items'][key] for key in keys}
        del manifest
        record['source_identity'].update(ionic_q1=array_identity(data.ionic),
            ionic_reverse_q1=array_identity(data.ionic_reverse), rho0=array_identity(c.rho0),
            W0_q1=array_identity(c.mechanical.w0[1]), W0_opposite_q=array_identity(c.mechanical.w0[int(c.minus[1])]),
            generator_q1=array_identity(c.mechanical.generators[1]),
            generator_opposite_q=array_identity(c.mechanical.generators[int(c.minus[1])]))
        record['opposite_q'] = int(c.minus[1])
        record['canonical_routes'] = dict(source=list(range(128)),target=targets,
            slot_spin=np.asarray(c.slot_spin).tolist(), slot_k=np.asarray(c.slot_k).tolist())
        for order in ORDERS:
            path = output/('sigma-'+order+'.npy')
            banks[order] = np.lib.format.open_memmap(path, mode='w+', dtype=np.complex128, shape=(31,31,128,12,12))
            handles[order] = path.open('r+b')
        save(output/'checkpoint.json',record)
        for a in range(31):
            for b in range(31):
                if stop or time.monotonic() >= deadline:
                    failure = dict(type='bounded-native-boundary', reason=stop or ['predeclared production deadline'],
                        next_pair=[a,b], elapsed_seconds=time.monotonic()-started)
                    break
                active_pair = [a,b]
                began = time.monotonic()
                forward, reverse = ref.electron(a,b), ref.electron(b,a)
                boson = ref.boson(1,a,b)
                # The owning loop and evaluators are unchanged, including all111 directions.
                for source_slot in range(c.ns):
                    target = int(c.target[1,source_slot])
                    pair = edge_pair(data,source_slot,target,forward,reverse,boson)
                    values = evaluate_components(pair,n=c.n,t=ref.times[a],s=ref.times[b],components=('sigma',))
                    for order in ORDERS:
                        banks[order][a,b,target] = values['sigma_'+order]
                for order in ORDERS:
                    banks[order].flush()
                    os.fsync(handles[order].fileno())
                record['coverage']['completed_pairs'].append([a,b])
                record['timings'].append(dict(pair=[a,b],seconds=time.monotonic()-began))
                record['elapsed_seconds'] = time.monotonic()-started
                save(output/'checkpoint.json',record)
                active_pair = None
                print(json.dumps(dict(event='actual-pair-committed',pair=[a,b],completed_pairs=len(record['coverage']['completed_pairs']),
                    elapsed_seconds=record['elapsed_seconds'],native_cpu_seconds=time.process_time())),flush=True)
            if failure:
                break
    except BaseException as exc:
        import traceback
        failure = dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
        record['boundary_loss']['incomplete_pair'] = active_pair
    finally:
        record['failure'] = failure
        record['elapsed_seconds'] = time.monotonic()-started
        record['native_cpu_seconds'] = time.process_time()
        record['peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
        record['u_integration_calls'] = ref.u_integration_calls
        record['gaussian_expm_calls'] = ref.expm_calls
        if len(record['coverage']['completed_pairs']) == 961:
            record['status'] = 'complete-owned-component-q1'
            record['completed_q'] = [1]
        record['sigma'] = {}
        for order,arr in banks.items():
            arr.flush()
            os.fsync(handles[order].fileno())
            record['sigma'][order] = dict(binding(output/('sigma-'+order+'.npy')),
                shape=list(arr.shape),dtype=str(arr.dtype))
            handles[order].close()
            arr._mmap.close()
        ref.close()
        save(output/'catalog.json',record)
        print(json.dumps(dict(event='component-published',status=record['status'],pairs=len(record['coverage']['completed_pairs']),failure=failure)),flush=True)
    return output/'catalog.json'
