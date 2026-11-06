"""One bounded actual saved-source and paid-provider consumer exercise."""
from pathlib import Path
import json
import os
import resource
import sys
import time
import traceback
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent/n) for n in
    ('0063','0062','0060','0052','0051','0050','0046','0044','0045','0041','0030','0029')]
import numpy as np
from history_storage import binding, authenticate, atomic_json, StoredNPZ
from saved_mean import SavedMean
from history_provider import PreviousHistory, SavedHistory
from history_readout import consume, json_value, raw_difference
from continue_history import validate_reference_identity
from history_sources import open_reference
from material_electron_ion_inputs import load_material
from produce_source import produce

OUT = Path('/var/tmp/p254-a68-genuine-bare-source-001')
PAID = Path('/var/tmp/p254-a63-full128-gamma-step-004/provider.json')


def save(path, value):
    path.write_text(json.dumps(value, default=json_value, indent=2)+'\n')


def main():
    started = time.monotonic()
    OUT.mkdir(exist_ok=False)
    group = Path('/proc/self/cgroup').read_text().strip().split('::')[-1]
    cg = Path('/sys/fs/cgroup')/group.lstrip('/')
    runtime = dict(affinity=sorted(os.sched_getaffinity(0)),
        threads={key:os.environ.get(key) for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS',
            'MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')},
        limits={key:(cg/key).read_text().strip() for key in ('memory.max','memory.swap.max')},
        HERDR={key:value for key,value in os.environ.items() if key.startswith('HERDR_')},
        argv=sys.argv, pid=os.getpid(), cgroup=group)
    try:
        runtime['cpuset.cpus.effective'] = (cg/'cpuset.cpus.effective').read_text().strip()
    except FileNotFoundError as error:
        runtime['cpuset.cpus.effective'] = dict(unavailable=str(error), evidence='sched_getaffinity')
    save(OUT/'runtime.json', runtime)
    if (runtime['affinity'] != [8] or set(runtime['threads'].values()) != {'1'} or
            runtime['limits'] != {'memory.max':'3221225472','memory.swap.max':'0'} or
            runtime['HERDR'].get('HERDR_PANE_ID') != 'w1:p79' or os.environ.get('HERDR_ENV') != '1'):
        raise ValueError('actual owned one-thread CPU8/3GiB/noSwap view differs')
    findings = []
    inputs = [binding(HERE/'preregistration.json'), binding(PAID)]
    # Exercise the original paid full128/full114 provider, not another solve.
    spec = json.loads(PAID.read_bytes())
    authenticate(spec['input_previous'])
    authenticate(spec['input_mean'])
    authenticate(spec['input_kernels'])
    inputs += [spec['input_previous'], spec['input_mean'], spec['input_kernels']]
    with SavedHistory(PAID) as paid, PreviousHistory(spec['input_previous']['path']) as previous, \
            SavedMean(spec['input_mean']['path']) as mean:
        catalog = json.loads(Path(spec['input_kernels']['path']).read_bytes())
        old_identity = validate_reference_identity(catalog, mean, previous)
        actual = json.loads(json.dumps(consume(paid, previous), default=json_value))
        expected_path = PAID.parent/'readout.json'
        expected = json.loads(expected_path.read_bytes())
        save(OUT/'state004-readout.json', actual)
        if actual != expected:
            raise AssertionError('paid state004 full readout changed')
        for name in ('control_ret','control_less'):
            with StoredNPZ(spec['input_previous']['path']) as legacy:
                if not np.array_equal(paid.control[name], legacy[name]):
                    raise AssertionError('paid original controls changed')
        findings.append(dict(claim='Paid state004 full128/full114 consumer unchanged',
            all_readout_values_equal=True, original_control_arrays_equal=True,
            source_readout=binding(expected_path), exercised_readout=binding(OUT/'state004-readout.json'),
            identity=old_identity, native_solves=0))
    descriptor = OUT/'source.json'
    source_binding = produce(descriptor)
    c, manifest, coefficient_binding = load_material()
    with SavedMean(descriptor) as mean, PreviousHistory(descriptor) as previous:
        if (not mean.is_bare or previous.provider is not None or previous.control is not None or
                mean.binding == mean.reference_identity or previous.reference_identity != mean.reference_identity):
            raise AssertionError('bare zero-iterate artifact/reference/control semantics differ')
        if (not np.array_equal(mean.times, np.arange(31, dtype=float)) or
                not np.array_equal(mean.rho0, c.rho0) or
                not np.array_equal(mean.w0, c.mechanical.w0[0]) or
                int(previous.aux['n']) != c.n):
            raise AssertionError('original all31 source preparation/normalization differs')
        bundle = Path(mean.reference_identity['path'])
        with_reference = open_reference(c, coefficient_binding, bundle=bundle)
        try:
            if not np.array_equal(mean.u, with_reference.u):
                raise AssertionError('actual original saved U differs')
            s_equal = all(np.array_equal(mean.s[a], with_reference.S(0,t))
                          for a,t in enumerate(mean.times))
            if not s_equal:
                raise AssertionError('actual original full-phase S differs')
            electron_delta = 0.
            electron_reference_delta = 0.
            # Every original sector and every stored ordered pair passes both A63 readers.
            first_column = [with_reference.electron(a,0) for a in range(31)]
            for slot in range(128):
                mr, ml = mean.electron_slot(slot)
                pr, pl = previous.electron_slot(slot)
                electron_delta = max(electron_delta, raw_difference(mr,pr), raw_difference(ml,pl))
                for a in range(31):
                    greater, lesser = first_column[a]
                    electron_reference_delta = max(electron_reference_delta,
                        float(np.max(np.abs(ml[a,0]-lesser[slot]))),
                        float(np.max(np.abs(mr[a,0]-(greater[slot]-lesser[slot])))))
            del mr, ml, pr, pl
            gr, gl = previous.gamma_slot()
            mr, ml = mean.gamma_slot()
            gamma_reader_delta = max(raw_difference(gr,mr), raw_difference(gl,ml))
            del mr, ml
            gamma_reference_delta, gamma_scale = 0., 0.
            for a in range(31):
                for b in range(31):
                    greater, lesser = with_reference.ordered(0,a,b,'positions')
                    gamma_reference_delta = max(gamma_reference_delta,
                        float(np.max(np.abs(gl[a,b]-lesser))))
                    gamma_scale = max(gamma_scale, float(np.max(np.abs(lesser))))
                    if a >= b:
                        gamma_reference_delta = max(gamma_reference_delta,
                            float(np.max(np.abs(gr[a,b]-(greater-lesser)))))
                        gamma_scale = max(gamma_scale,float(np.max(np.abs(greater-lesser))))
            tolerance = 64*np.finfo(float).eps*228*max(1.,gamma_scale)
            if (electron_delta != 0. or electron_reference_delta > 1e-12 or
                    gamma_reader_delta != 0. or gamma_reference_delta > tolerance):
                raise AssertionError('real A63/A62 source transport differs beyond binary64 contraction allowance')
            bare_readout = consume(previous, previous)
            save(OUT/'bare-readout.json', bare_readout)
            if bare_readout['represented_control']['status'] != 'unavailable':
                raise AssertionError('absent control became synthetic data')
            findings.append(dict(claim='Actual all31 original bare zero-iterate saved/read through A63 and A62',
                source=source_binding, artifact_identity=mean.binding,
                numerical_reference_identity=mean.reference_identity,
                all128_slots_all961_pairs_A63_maxabs=electron_delta,
                all128_slots_all31_first_column_pairs_A62_maxabs=electron_reference_delta,
                full114_all961_pairs_A63_maxabs=gamma_reader_delta,
                full114_all961_pairs_A62_maxabs=gamma_reference_delta,
                gamma_maxabs_scale=gamma_scale, gamma_rounding_allowance=tolerance,
                raw_U_equal=True, raw_S_all31_equal=s_equal, raw_rho0_W0_n_equal=True,
                original_literal_U0_S0=True, U_integration_calls=with_reference.u_integration_calls,
                material_trajectory_replays=0, native_history_solves=0,
                represented_control=bare_readout['represented_control'],
                readout=binding(OUT/'bare-readout.json')))
            inputs += mean.bindings + [coefficient_binding]
        finally:
            with_reference.close()
    result = dict(schema='p254-a68-genuine-bare-source-result-v1',status='complete',
        findings=findings, input_bindings=inputs,
        source_bindings=[binding(path) for directory in (HERE,HERE.parent/'0063')
                         for path in sorted(directory.glob('*.py'))],
        output_bindings=[binding(OUT/name) for name in
            ('source.json','runtime.json','state004-readout.json','bare-readout.json')],
        elapsed_seconds=time.monotonic()-started,
        peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        limits=['Prescribed bare zero-iterate, not common mean feedback or a finalized solved bank.',
            'No new kernel computation or full intrinsic evolution; no convergence, temporal accuracy or conservation closure.',
            'No new represented control consumed; no physical ranking, rates, loss or heat.'])
    save(OUT/'result.json', result)
    print('A68_SOURCE_NATIVE_COMPLETE '+json.dumps(binding(OUT/'result.json')),flush=True)


if __name__ == '__main__':
    code = 0
    started = time.monotonic()
    try:
        main()
    except BaseException:
        code = 1
        traceback.print_exc()
    finally:
        if OUT.exists():
            save(OUT/'native-exit.json',dict(returncode=code, elapsed_seconds=time.monotonic()-started,
                maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
        print('A68_SOURCE_NATIVE_TERMINAL_'+str(code),flush=True)
    raise SystemExit(code)
