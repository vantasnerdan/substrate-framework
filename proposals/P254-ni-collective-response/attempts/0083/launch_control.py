"""One actual explicitly leased A83 control case; no restarts or retries."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent


def authenticate(item):
    raw = Path(item['path']).read_bytes()
    if len(raw) != item['bytes'] or hashlib.sha256(raw).hexdigest() != item['sha256']:
        raise ValueError('actual declared control/lease/source binding differs: '+item['path'])
    return raw


def main(path):
    case = json.loads(Path(path).read_bytes())
    for item in [case['authorization'], case['preregistration'], *case['source_bindings']]:
        authenticate(item)
    context = {k: v for k, v in os.environ.items() if k.startswith('HERDR_')}
    if (context.get('HERDR_ENV') != '1' or context.get('HERDR_PANE_ID') != case['pane_id'] or
            context.get('HERDR_TAB_ID') != case['tab_id'] or case['cpu'] not in range(6, 11) or
            case['cap_seconds'] != 360):
        raise ValueError('actual explicitly assigned bounded control view/CPU required')
    contract = json.loads((HERE.parent/'0045/implementation-contract-001.json').read_bytes())['common_runtime']
    receipts = Path(case['receipts'])
    receipts.mkdir(parents=True, exist_ok=False)
    unit = case['unit']
    argv = contract['manager_env']+['systemd-run', '--user', '--wait', '--pipe', '--unit='+unit,
        '--property=RuntimeMaxSec=360', '--property=MemoryMax=3221225472',
        '--property=MemorySwapMax=0', '--property=Restart=no',
        '--property=CPUAffinity='+str(case['cpu']), 'env', contract['python_environment_assignment'],
        *[k+'='+v for k, v in contract['threads'].items()],
        *[k+'='+v for k, v in context.items()], 'PYTHONDONTWRITEBYTECODE=1',
        '/usr/bin/time', '-v', '-o', str(receipts/'time.txt'), contract['python'], '-u',
        str(HERE/'projected_mean_control.py'), str(Path(path).resolve())]

    def save(name, record):
        with (receipts/name).open('x') as stream:
            json.dump(record, stream, indent=2)
            stream.write('\n')

    started = time.monotonic()
    save('launch.json', dict(argv=argv, unit=unit, context=context,
        actual_start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())))
    code = None
    try:
        with (receipts/'stdout.txt').open('xb') as stdout, (receipts/'stderr.txt').open('xb') as stderr:
            process = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=stderr)
            for line in process.stdout:
                stdout.write(line)
                stdout.flush()
                print(line.decode(errors='replace'), end='', flush=True)
            code = process.wait()
    finally:
        manager = contract['manager_env']+['systemctl', '--user', 'show', unit,
            '--property=Id,InvocationID,Result,ExecMainCode,ExecMainStatus,CPUUsageNSec,MemoryPeak,MemorySwapPeak,MemoryMax,MemorySwapMax,Restart,RuntimeMaxUSec,ActiveState,SubState']
        terminal = subprocess.run(manager, capture_output=True, check=False)
        save('manager-terminal.json', dict(argv=manager, returncode=terminal.returncode,
            stdout=terminal.stdout.decode(), stderr=terminal.stderr.decode()))
        save('exit.json', dict(argv=argv, returncode=code, elapsed_seconds=time.monotonic()-started,
            actual_end_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())))
        print('A83_OWNER_TERMINAL_'+str(code)+' '+unit, flush=True)
    return code if code is not None else 1


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1]))
