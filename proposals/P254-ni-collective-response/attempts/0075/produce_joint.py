"""A62 lines116–125 only: global full-Gamma minus ionic-Gamma Sigma, full Pi once."""
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
from dressed_functionals import common_pair
from joint_gamma_adapter import vertices


def identity(value):
    value = np.ascontiguousarray(value)
    return dict(shape=list(value.shape), dtype=str(value.dtype), bytes=value.nbytes,
        sha256=hashlib.sha256(memoryview(value).cast('B')).hexdigest(), hash_scope='C-order numeric payload')


def produce(output, deadline):
    started = time.monotonic()
    prereg = json.loads((HERE/'preregistration.json').read_bytes())
    source = authenticate(prereg['genuine_reference'])
    authenticate(prereg['existing_complete_joint_loop'])
    authenticate(prereg['source_descriptor'])
    output = Path(output)
    output.mkdir(exist_ok=False)
    stop = []
    signal.signal(signal.SIGTERM, lambda number, frame: stop.append('SIGTERM'))
    c, baseline, old, material_bindings = load_inputs()
    old.close()
    del old, baseline
    ref = open_reference(c, material_bindings['manifest'], bundle=source)
    banks, handles = {}, {}
    pairs = [[a,b] for a in range(31) for b in range(31)]
    record = dict(schema='p254-additive-history-component-v1', component='joint', q_additive=False,
        status='partial-global-joint-component', full_mesh=False, full_solve_eligible=False,
        full_time_square=False, input_dataset=binding(source), times_inverse_ry=ref.times.tolist(),
        slots=128, bands=12, spins=2, ionic_directions=111, gamma_coordinates=114,
        axes=dict(sigma=['time','time','slot','band','band'], pi=['time','time','coordinate','coordinate']),
        sigma_convention='joint_minus_ionic_Gamma', pi_convention='full_joint_Gamma_once',
        normalization='original canonical vertices: ionic/sqrt(N), field modes normalized; evaluate_components n=1',
        source_bindings=ref.bindings, coefficient_bindings=material_bindings,
        source_descriptor=prereg['source_descriptor'], owning_implementation=prereg['existing_complete_joint_loop'],
        implementation_bindings=[binding(HERE/name) for name in ('produce_joint.py','component_reader.py','run_native.py')],
        component_implementation_bindings=[binding(Path(__import__(name).__file__)) for name in
            ('dressed_functionals','joint_gamma_adapter','electron_ion_functional','joint_material_inputs','history_sources','general_moving_reference')],
        coverage=dict(owned_pairs=pairs, completed_pairs=[], external_targets=list(range(128)),
            all114_coordinates=True, checkpoint='all six numeric pair slices flushed and fsynced before coverage commit'),
        preparation=ref.spec['initial_conditions'], paid_node_reuse=[], timings=[],
        scope='Once-only global jointGamma correction Sigma and full jointGamma Pi. No q partition/sum, raw EI, EE, cubic, complete all64q functional or simulation.',
        boundary_loss=dict(complete_pairs_lost=0, incomplete_pair=None))
    failure = None
    active_pair = None
    try:
        if not np.array_equal(ref.times, np.arange(31)) or (c.ns,c.m,c.d,c.f,c.n)!=(128,12,111,3,64):
            raise ValueError('frozen full31 original physical axes differ')
        forward, reverse = vertices(c)
        record['physical_source_identity'] = dict(U=dict(ref.spec['U'],shape=list(ref.u.shape),dtype=str(ref.u.dtype)),
            rho0=identity(c.rho0), W0_Gamma=identity(c.mechanical.w0[0]),
            generator_Gamma=identity(c.mechanical.generators[0]),
            canonical_forward_vertex=identity(forward), independent_reverse_vertex=identity(reverse))
        del forward, reverse
        record['canonical_routes'] = dict(source=list(range(128)), target=np.asarray(c.target[0]).tolist(),
            slot_spin=np.asarray(c.slot_spin).tolist(), slot_k=np.asarray(c.slot_k).tolist())
        for kind, tail in (('sigma',(128,12,12)),('pi',(114,114))):
            banks[kind], handles[kind] = {}, {}
            for order in ORDERS:
                path = output/(kind+'-'+order+'.npy')
                banks[kind][order] = np.lib.format.open_memmap(path, mode='w+', dtype=np.complex128, shape=(31,31,*tail))
                handles[kind][order] = path.open('r+b')
        save(output/'checkpoint.json',record)
        for a,b in pairs:
            if stop or time.monotonic() >= deadline:
                failure = dict(type='bounded-native-boundary',reason=stop or ['production deadline'],next_pair=[a,b])
                break
            active_pair = [a,b]
            began = time.monotonic()
            # Exact A62 owning calls, preserving both independent electron legs and normalized vertices.
            full = common_pair(c,ref,a,b)
            ionic = common_pair(c,ref,a,b,True)
            for order in ORDERS:
                banks['sigma'][order][a,b] = full['sigma_'+order] - ionic['sigma_'+order]
                banks['pi'][order][a,b] = full['pi_'+order]
            del full, ionic
            for kind in banks:
                for order in ORDERS:
                    banks[kind][order].flush()
                    os.fsync(handles[kind][order].fileno())
            record['coverage']['completed_pairs'].append([a,b])
            record['timings'].append(dict(pair=[a,b],seconds=time.monotonic()-began))
            record['elapsed_seconds'] = time.monotonic()-started
            save(output/'checkpoint.json',record)
            active_pair = None
            print(json.dumps(dict(event='actual-global-joint-pair-committed',pair=[a,b],
                completed_pairs=len(record['coverage']['completed_pairs']),elapsed_seconds=record['elapsed_seconds'],
                native_cpu_seconds=time.process_time())),flush=True)
    except BaseException as exc:
        import traceback
        failure = dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
        record['boundary_loss']['incomplete_pair'] = active_pair
    finally:
        record.update(failure=failure, elapsed_seconds=time.monotonic()-started,
            native_cpu_seconds=time.process_time(),peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
            u_integration_calls=ref.u_integration_calls,gaussian_expm_calls=ref.expm_calls)
        if len(record['coverage']['completed_pairs']) == 961:
            record.update(status='complete-global-joint-component',full_time_square=True)
        for kind,bank in banks.items():
            record[kind] = {}
            for order,arr in bank.items():
                arr.flush()
                os.fsync(handles[kind][order].fileno())
                record[kind][order] = dict(binding(output/(kind+'-'+order+'.npy')),shape=list(arr.shape),dtype=str(arr.dtype))
                handles[kind][order].close()
                arr._mmap.close()
        ref.close()
        save(output/'catalog.json',record)
        print(json.dumps(dict(event='component-published',status=record['status'],pairs=len(record['coverage']['completed_pairs']),failure=failure)),flush=True)
    return output/'catalog.json'
