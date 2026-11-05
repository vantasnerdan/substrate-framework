"""Durable exact-argv execution and byte bindings for the owned Native boundary."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time


def binding(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return {'path': str(path), 'bytes': path.stat().st_size, 'sha256': digest.hexdigest()}


def run_recorded(argv, directory, name, timeout=590):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    record = {'argv': [str(a) for a in argv], 'cwd': os.getcwd(),
              'started_unix': time.time(), 'timeout_seconds': timeout}
    receipt = directory / (name + '.json')
    receipt.write_text(json.dumps(record, indent=2) + '\n')
    started = time.monotonic()
    with (directory / (name + '.stdout')).open('wb') as out, (directory / (name + '.stderr')).open('wb') as err:
        try:
            result = subprocess.run(record['argv'], stdout=out, stderr=err, timeout=timeout, check=False)
            record['exit'] = result.returncode
        except subprocess.TimeoutExpired:
            record['exit'] = None
            record['failure'] = 'timeout'
        except OSError as exc:
            record['exit'] = None
            record['failure'] = repr(exc)
    record['wall_seconds'] = time.monotonic() - started
    record['stdout'] = binding(directory / (name + '.stdout'))
    record['stderr'] = binding(directory / (name + '.stderr'))
    receipt.write_text(json.dumps(record, indent=2) + '\n')
    if record['exit'] != 0:
        raise RuntimeError(f'command failed; exact receipt {receipt}')
    return record
