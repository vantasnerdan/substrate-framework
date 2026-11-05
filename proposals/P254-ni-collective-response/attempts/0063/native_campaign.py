"""Bounded CPU10/1GiB/300s actual slot, using the existing A61 binary."""
from pathlib import Path
import os
import resource
import sys
import time
from history_storage import HERE, binding, atomic_json
from production_transport import run_recorded

PYTHON = '/home/administrator/projects/substrate-framework/.venv/bin/python'
BINARY = Path('/var/tmp/p254-a61-lean-002/build/solve_result')
ROOT = Path('/var/tmp/p254-a63-history-slot-001')
UNIT = 'p254-a63-history-slot-001.service'


def native(root):
    started = time.monotonic()
    environment = {key: os.environ.get(key) for key in (
        'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
        'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS')}
    receipt = dict(affinity=sorted(os.sched_getaffinity(0)), herdr_env=os.environ.get('HERDR_ENV'),
        environment=environment, sources=[binding(p) for p in sorted(HERE.glob('*.py'))],
        view=binding(HERE/'native-view.json'), binary=binding(BINARY),
        upstream=[binding(HERE.parent/n/file) for n, file in (
            ('0061', 'production_transport.py'), ('0062', 'history_banks.py'),
            ('0054', 'supplied_keldysh.hpp'), ('0049', 'supplied_vie2.hpp'),
            ('0049', 'command_receipts.py'), ('0053', 'consume_material_step.py'),
            ('0029', 'reduced_coherence_readout.py'))])
    atomic_json(root/'native-bindings.json', receipt)
    code = 1
    try:
        if receipt['affinity'] != [10] or receipt['herdr_env'] != '1' or any(v != '1' for v in environment.values()):
            raise RuntimeError('actual assigned grant/environment mismatch')
        from smoke_saved_slot import exercise
        findings = exercise(BINARY, root)
        receipt['findings'] = binding(root/'actual-findings.json')
        receipt['identity'] = dict(ret= findings['retarded_bitwise_equal_to_paid_solution'],
                                  less=findings['lesser_bitwise_equal_to_paid_solution'])
        code = 0
    except BaseException as exc:
        receipt['failure'] = repr(exc)
        raise
    finally:
        own, children = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
        receipt.update(exit=code, wall_seconds=time.monotonic()-started,
            cpu_seconds=own.ru_utime+own.ru_stime+children.ru_utime+children.ru_stime,
            owner_maxrss_kib=own.ru_maxrss, child_maxrss_kib=children.ru_maxrss,
            output_bytes=sum(p.stat().st_size for p in root.rglob('*') if p.is_file()))
        atomic_json(root/'native-outcome.json', receipt)
        print('[dft-done]' if code == 0 else '[dft-fail]', 'A63 native exit='+str(code), flush=True)


def launch(root):
    if root != ROOT:
        raise ValueError('undeclared staging root')
    root.mkdir(exist_ok=False)
    argv = ['env', 'XDG_RUNTIME_DIR=/run/user/1000',
        'DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/1000/bus',
        'systemd-run', '--user', '--wait', '--pipe', '--collect', '--unit='+UNIT,
        '--property=RuntimeMaxSec=295', '--property=MemoryMax=1073741824',
        '--property=MemorySwapMax=0', '--property=Restart=no', '--property=StartLimitBurst=1',
        '--property=CPUAffinity=10', 'env', 'HERDR_ENV=1', 'OPENBLAS_NUM_THREADS=1',
        'OMP_NUM_THREADS=1', 'MKL_NUM_THREADS=1', 'NUMEXPR_NUM_THREADS=1',
        'BLIS_NUM_THREADS=1', 'PYTHONDONTWRITEBYTECODE=1', 'taskset', '-c', '10',
        PYTHON, '-u', str(Path(__file__).resolve()), '--native', str(root)]
    code = 1
    try:
        run_recorded(argv, root, 'systemd-native', 300)
        code = 0
    finally:
        atomic_json(root/'terminal.json', dict(exit=code, root=str(root),
            unit=UNIT, marker='A63_NATIVE_TERMINAL_EXIT_'+str(code)))
        print('A63_NATIVE_TERMINAL_EXIT_'+str(code), flush=True)


if __name__ == '__main__':
    if sys.argv[1] == '--native':
        native(Path(sys.argv[2]))
    else:
        launch(Path(sys.argv[1]))
