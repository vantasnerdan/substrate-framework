"""A79 CPU7 saved A68/A69 readout only; no kernel or control production."""
from pathlib import Path
import json
import os
import resource
import subprocess
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
PANE = 'w1:p7M'
LABEL = 'A79 actual A68bare31 vs original A69 projected control'


def dump(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def native():
    started = time.monotonic()
    findings = dict(schema='p254-a79-native-findings-v1', status='running', findings=[], blockers=[])
    try:
        sys.path.insert(0, str(HERE.parent/'0063'))
        from history_storage import authenticate, binding
        from history_provider import PreviousHistory
        from history_readout import consume, json_value
        runtime = dict(affinity=sorted(os.sched_getaffinity(0)),
            environment={key: value for key, value in os.environ.items()
                         if key.startswith('HERDR_') or key in (
                'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
                'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS', 'PYTHONPATH')})
        cg = Path('/sys/fs/cgroup')/Path('/proc/self/cgroup').read_text().strip().split('::')[-1].lstrip('/')
        runtime['cgroup'] = {key: (cg/key).read_text().strip() for key in ('memory.max', 'memory.swap.max')}
        if (runtime['affinity'] != [7] or os.environ.get('HERDR_ENV') != '1'
                or os.environ.get('HERDR_PANE_ID') != PANE):
            raise RuntimeError('actual owned Herdr view and CPU7 required')
        if runtime['cgroup'] != {'memory.max': '2147483648', 'memory.swap.max': '0'}:
            raise RuntimeError('assigned 2GiB/noSwap cgroup required')
        if any(os.environ.get(key) != '1' for key in (
                'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
                'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS')):
            raise RuntimeError('one-thread environment required')
        dump(HERE/'native-startup.json', runtime)
        frozen = json.loads((HERE/'preregistration.json').read_bytes())
        source = authenticate(frozen['actual_bare_history'])
        handoff_path = HERE.parent/'0069/handoff.json'
        control_handoff = json.loads(handoff_path.read_bytes())
        declared = {key: control_handoff['control'][key] for key in ('path', 'bytes', 'sha256')}
        declared['path'] = str((handoff_path.parent/declared['path']).resolve())
        control = authenticate(declared)
        sources = [binding(HERE/'preregistration.json'), binding(HERE/'run_native.py'),
            binding(HERE.parent/'0063/history_readout.py'), binding(HERE.parent/'0063/continue_history.py'),
            binding(HERE.parent/'0063/history_provider.py'), binding(HERE.parent/'0063/saved_mean.py'),
            binding(HERE.parent/'0063/history_storage.py'), binding(handoff_path),
            frozen['actual_bare_history'], declared]
        dump(HERE/'actual-invocation.json', dict(schema='p254-a79-saved-readout-invocation-v1',
            argv=sys.argv, source_bindings=sources,
            public_call='history_readout.consume(actual_A68_history, same_actual_A68_history, projected_control=actual_A69_control)',
            no_kernel_or_control_production=True))
        with PreviousHistory(source) as history:
            readout = consume(history, history, projected_control=control)
            if (len(readout['times_inverse_ry']) != 31 or
                    readout['projected_comparison']['status'] != 'available' or
                    readout['represented_control']['status'] != 'unavailable'):
                raise RuntimeError('actual31-node separate projected-only readout required')
            with (HERE/'saved-path-readout.json').open('x') as stream:
                json.dump(readout, stream, default=json_value, indent=2, allow_nan=False)
                stream.write('\n')
            findings['actual_observations'] = dict(
                times_inverse_ry=readout['times_inverse_ry'].tolist(),
                projected_comparison=json.loads(json.dumps(readout['projected_comparison'], default=json_value)),
                legacy_full114_control=readout['represented_control'],
                original_bare=json.loads(json.dumps(readout['original_bare'], default=json_value)),
                preparation=dict(fermion_maxabs=readout['fermion_initial_source_maxabs'],
                    gamma_maxabs=readout['gamma_initial_source_maxabs']))
        findings['readout'] = binding(HERE/'saved-path-readout.json')
        findings['invocation'] = binding(HERE/'actual-invocation.json')
        findings['findings'] = [
            dict(status='exercised', claim='Actual A68 bare31 saved factors consumed by PreviousHistory and original A53 formulas; actual A69 saved NPZ consumed through its public normalizer.'),
            dict(status='validated', claim='Actual common original reference/manifest/baseline, all259 consumed archive identities, source layout, full128 rho0, full114 position W0, literal endpoint probe and N64.'),
            dict(status='separate', claim='Original represented coupled Gaussian projected baseline only; absent legacyfull114 control stays unavailable. No interpolated or synthetic full114 history.'),
            dict(status='carried', claim='Original source-specific control proof and all original suites reused without replay; no new Gaussian control solve or kernel production.')]
        findings['limits'] = ['Not matched evolved-common-mean control or proof of intrinsic-only loss.',
            'A68 bare zero-iterate is not an evolved intrinsic/common-mean result.',
            'No kernel recomputation, fullmesh/fixedpoint, physical calibration or physical closure claim.',
            'No new time-refinement or error enclosure; exact endpoint readout is not temporal accuracy.',
            'Full continuation Dyson solve was not run; the authorized affected exercise is saved-path readout.']
        findings['status'] = 'complete-at-stated-scope'
    except BaseException as exc:
        findings['status'] = 'failed'
        findings['blockers'].append(dict(type=type(exc).__name__, message=str(exc), traceback=traceback.format_exc()))
        raise
    finally:
        findings['native_elapsed_seconds'] = time.monotonic()-started
        findings['peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
        dump(HERE/'native-findings.json', findings)
        print('A79_NATIVE_TERMINAL '+findings['status'], flush=True)


def launch():
    receipts = HERE/'receipts-001'
    receipts.mkdir(exist_ok=False)
    contract = json.loads((HERE.parent/'0045/implementation-contract-001.json').read_bytes())['common_runtime']
    context = {key: value for key, value in os.environ.items() if key.startswith('HERDR_')}
    if (context.get('HERDR_ENV') != '1' or context.get('HERDR_PANE_ID') != PANE
            or not all(context.get(key) for key in ('HERDR_TAB_ID', 'HERDR_WORKSPACE_ID'))):
        raise RuntimeError('actual owned view context missing; never synthesize Herdr environment')
    unit = 'p254-a79-projected-readout-001.service'
    argv = contract['manager_env']+['systemd-run', '--user', '--wait', '--pipe', '--unit='+unit,
        '--property=RuntimeMaxSec=235', '--property=MemoryMax=2147483648',
        '--property=MemorySwapMax=0', '--property=Restart=no', '--property=CPUAffinity=7',
        'env', contract['python_environment_assignment'],
        *[key+'='+value for key, value in contract['threads'].items()],
        *[key+'='+value for key, value in context.items()], 'PYTHONDONTWRITEBYTECODE=1',
        '/usr/bin/time', '-v', '-o', str(receipts/'time.txt'), 'taskset', '-c', '7',
        contract['python'], '-u', str(Path(__file__).resolve()), '--native']
    started, wall = time.monotonic(), time.time()
    dump(receipts/'launch.json', dict(argv=argv, start_unix=wall, context=context, unit=unit,
        grant_seconds=240, systemd_cap_seconds=235, owned_view=PANE, label=LABEL,
        expected_useful_output='actual A68 saved31-history vs original A69 one-probe raw/normalized endpoint comparison',
        prior_costs='Original A69 paid control and proofs carried unchanged; no old suite or producer replay.'))
    rc = None
    try:
        with (receipts/'stdout.txt').open('xb') as stdout, (receipts/'stderr.txt').open('xb') as stderr:
            rc = subprocess.run(argv, stdout=stdout, stderr=stderr, check=False).returncode
    finally:
        query = contract['manager_env']+['systemctl', '--user', 'show', unit,
            '--property=Id,InvocationID,Result,ExecMainCode,ExecMainStatus,CPUUsageNSec,MemoryPeak,MemorySwapPeak,MemoryMax,MemorySwapMax,Restart,RuntimeMaxUSec,ActiveState,SubState']
        manager = subprocess.run(query, capture_output=True, check=False)
        dump(receipts/'manager-terminal.json', dict(argv=query, exitcode=manager.returncode,
            stdout=manager.stdout.decode(), stderr=manager.stderr.decode()))
        dump(receipts/'exit.json', dict(argv=argv, exitcode=rc, start_unix=wall,
            end_unix=time.time(), elapsed_seconds=time.monotonic()-started))
        print('A79_OWNED_LAUNCH_TERMINAL rc='+str(rc), flush=True)
    return rc


if __name__ == '__main__':
    if '--native' in sys.argv:
        native()
    else:
        sys.exit(launch())
