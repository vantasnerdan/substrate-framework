"""Own one new physical-operator probe, preserving inputs and all-exit outcomes."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import threading
import time


def save(path, value):
    temporary = path.with_name(path.name + '.part')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def sample_native(executable, stop, records):
    # Only inspect this owned service leaf, never unrelated/protected processes.
    relative = next(line.split(':', 2)[2] for line in Path('/proc/self/cgroup').read_text().splitlines()
                    if line.startswith('0::'))
    members = Path('/sys/fs/cgroup') / relative.lstrip('/') / 'cgroup.procs'
    while not stop.wait(.02):
        try:
            pids = members.read_text().split()
        except FileNotFoundError:
            return
        for text in pids:
            pid = int(text)
            try:
                if Path(os.readlink(f'/proc/{pid}/exe')).resolve() != executable:
                    continue
                records[text] = {p.name: sorted(os.sched_getaffinity(int(p.name)))
                                 for p in Path(f'/proc/{pid}/task').iterdir()}
            except (FileNotFoundError, ProcessLookupError):
                continue


def run_case(case_path):
    raw = case_path.read_bytes()
    case = json.loads(raw)
    for entry in case['inputs']:
        data = Path(entry['path']).read_bytes()
        if len(data) != int(entry['bytes']) or hashlib.sha256(data).hexdigest() != entry['sha256']:
            raise ValueError('Changed physical probe input: ' + entry['path'])
    if sorted(os.sched_getaffinity(0)) != [case['cpu']]:
        raise ValueError('Owner affinity differs from fresh grant')
    stop, affinity = threading.Event(), {}
    sampler = threading.Thread(target=sample_native,
        args=(Path(case['native_executable']).resolve(), stop, affinity))
    started = time.monotonic()
    process = subprocess.Popen(case['argv'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    sampler.start()
    try:
        stdout, stderr = process.communicate(timeout=150)
    except subprocess.TimeoutExpired as error:
        Path('owner.stdout').write_bytes(error.stdout or b'')
        Path('owner.stderr').write_bytes(error.stderr or b'')
        save(Path('owner-result.json'), {'status': 'case_timeout', 'native_affinity': affinity})
        raise SystemExit('Native probe budget reached; systemd owns whole-cgroup cleanup.')
    finally:
        stop.set()
        sampler.join()
    Path('owner.stdout').write_bytes(stdout)
    Path('owner.stderr').write_bytes(stderr)
    print(stdout.decode(), end='', flush=True)
    if process.returncode:
        raise RuntimeError(stderr.decode() + f' producer exit {process.returncode}')
    if stdout.count(b'JOB DONE.') != 1:
        raise ValueError('Native successful-exit marker not observed exactly once')
    if not affinity or any(mask != [case['cpu']] for tasks in affinity.values() for mask in tasks.values()):
        raise ValueError('Actual Native singleton-affinity evidence missing or wrong')
    closed = [json.loads(line) for line in stdout.decode().splitlines()
              if line.startswith('{"closed_q_bundles":')]
    if len(closed) != 1:
        raise ValueError('Missing unique producer closure result')
    save(Path('owner-result.json'), {'status': 'complete', 'native_returncode': 0,
        'elapsed_seconds': time.monotonic()-started, 'native_affinity': affinity,
        'case_input_sha256': hashlib.sha256(raw).hexdigest(), 'case': case,
        'closed_q_bundles': closed[0]['closed_q_bundles'],
        'scope': 'Constructed operator data; scientific kernel comparison is separate.'})


def notify(unit):
    outcome = {'unit': unit, 'actual_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'invocation': os.environ.get('INVOCATION_ID'), 'result': os.environ.get('SERVICE_RESULT'),
        'exit_code': os.environ.get('EXIT_CODE'), 'exit_status': os.environ.get('EXIT_STATUS'),
        'scope': 'Physical response-operator probe exit, not coherence-domain completion.'}
    path = Path('outcome-' + str(outcome['invocation']) + '.json')
    save(path, outcome)
    tag = 'dft-done' if outcome['result'] == 'success' else 'dft-fail'
    result = subprocess.run(['/usr/bin/python3',
        '/home/administrator/projects/catalyst-dft/runs/campaign3/notify.py', tag, unit, str(path.resolve())],
        capture_output=True, timeout=30)
    Path('notification.stdout').write_bytes(result.stdout)
    Path('notification.stderr').write_bytes(result.stderr)
    save(Path('notification-result.json'), {'exit': result.returncode})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', type=Path)
    parser.add_argument('--notify-exit')
    args = parser.parse_args()
    if args.notify_exit:
        notify(args.notify_exit)
    elif args.case:
        run_case(args.case)
    else:
        parser.error('Supply --case or --notify-exit')
