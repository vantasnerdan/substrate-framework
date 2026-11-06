"""Bound metadata writes through publication; owner outcome supplies exit authority.

A terminal is a pre-publication snapshot, not a measurement of its own later
fsync/exit. Consumers also require the successful post-exit manager outcome.
No failure fallback may exceed the same output/time grant. If even a compact
receipt cannot fit, journal the refusal and let the independent owner record
failure. Existing scientific checkpoints are never removed or rewritten.
"""
import json
import os
from pathlib import Path
import signal
import subprocess
import time
import uuid


BUSCTL = Path('/usr/bin/busctl')

FAILURE_RESERVE_BYTES = 4096

def manager_runtime_limit(context):
    """Observe only this process's own user-service deadline; retain uint64 text."""
    unit = context['cgroup'].rsplit('/', 1)[-1]
    if not context['INVOCATION_ID'] or not unit.endswith('.service'):
        return {'availability': 'unavailable', 'reason': 'not an own user-service invocation'}
    encoded = ''.join(chr(n) if 48 <= n <= 57 or 65 <= n <= 90 or 97 <= n <= 122
                      else '_'+format(n, '02x') for n in unit.encode())
    argv = [str(BUSCTL), '--user', 'get-property', 'org.freedesktop.systemd1',
            '/org/freedesktop/systemd1/unit/'+encoded,
            'org.freedesktop.systemd1.Service', 'RuntimeMaxUSec']
    try:
        result = subprocess.run(argv, capture_output=True, timeout=5, check=True)
        fields = result.stdout.decode().split()
        if len(fields) != 2 or fields[0] != 't' or not fields[1].isdigit():
            raise ValueError('unexpected typed RuntimeMaxUSec response')
        return {'availability': 'available', 'value': fields[1], 'argv': argv}
    except (OSError, subprocess.SubprocessError, ValueError) as error:
        return {'availability': 'unavailable', 'reason': str(error)}


def refusal(summary, error):
    """Bounded journal-only fallback; never fabricate a successful terminal."""
    record = {'event': 'publication_refused', 'run_id': summary.get('run_id'),
              'status': 'failed', 'exit_code': 1, 'error_type': type(error).__name__,
              'error': str(error)[:240]}
    print(json.dumps(record, separators=(',', ':')), flush=True)


class Publication:
    def __init__(self, target, summary, boot, grant, *, context_only=False):
        self.target, self.summary, self.boot = Path(target), summary, boot
        self.context_only = context_only
        # Explicitly diagnostic-only bounds; never authority for material work.
        self.grant = ({'remaining_seconds': 30, 'output_max_bytes': 8 << 20,
                       'first_block_seconds': 30} if context_only else grant)
        if self.grant is None:
            raise ValueError('publication requires admitted output/time authority')

    def elapsed(self):
        return time.monotonic()-self.boot

    def limit(self):
        limit = self.grant['remaining_seconds']
        if not self.context_only and self.summary.get('first_useful_block_seconds') is None:
            limit = min(limit, self.grant['first_block_seconds'])
        return limit

    def check_time(self):
        if self.elapsed() >= self.limit():
            raise TimeoutError('publication would exceed admitted invocation deadline')

    def arm(self):
        self.check_time()
        signal.setitimer(signal.ITIMER_REAL, self.limit()-self.elapsed())

    def used(self):
        return sum(p.stat().st_size for p in self.target.iterdir() if p.is_file())

    def encode(self, document, *, reserve=FAILURE_RESERVE_BYTES):
        self.check_time()
        raw = (json.dumps(document, sort_keys=True, separators=(',', ':'),
                          allow_nan=False)+'\n').encode()
        self.check_time()
        if self.used()+len(raw)+reserve > self.grant['output_max_bytes']:
            raise RuntimeError('metadata publication exceeds output allowance')
        return raw

    def write(self, path, document, *, rotate=False, terminal=False):
        path = Path(path)
        if path.parent != self.target:
            raise ValueError('metadata target outside admitted directory')
        if path.exists() and not rotate:
            raise FileExistsError('immutable metadata already exists: '+str(path))
        raw = self.encode(document, reserve=0 if terminal else FAILURE_RESERVE_BYTES)
        temporary = path.with_name(path.name+'.part-'+uuid.uuid4().hex)
        published = False
        try:
            with temporary.open('xb') as stream:
                stream.write(raw)
                stream.flush()
                os.fsync(stream.fileno())
            self.check_time()
            if rotate and path.exists():
                path.replace(path.with_name(path.name+'.prev'))
            temporary.replace(path)
            published = True
            self.check_time()
        except BaseException:
            # Withdraw only this new immutable receipt, never an old checkpoint.
            if published and not rotate:
                path.unlink()
            raise
        finally:
            if temporary.exists():
                temporary.unlink()

    def event(self, event, **fields):
        record = {'event': event, 'generation': self.summary['generation'],
                  'run_id': self.summary['run_id'], 'elapsed_seconds': self.elapsed(),
                  **fields}
        raw = self.encode(record)
        with (self.target/'telemetry.jsonl').open('ab') as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        self.check_time()
        print(raw.decode().rstrip(), flush=True)
        self.check_time()

    def finish(self, code):
        terminal = self.target/('terminal-'+self.summary['run_id']+'.json')
        published = False
        try:
            self.arm()  # Do not disarm before serialization, fsync or terminal publication.
            self.summary.update(elapsed_seconds=self.elapsed(), exit_code=code,
                time_accounting='Pre-publication snapshot; successful post-exit manager outcome also required.')
            if code:
                self.summary['status'] = 'failed'
            self.write(self.target/'summary.json', self.summary, rotate=True)
            self.event('finalizing', status=self.summary['status'], exit_code=code)
            self.write(terminal, self.summary, terminal=True)
            published = True
            print(json.dumps({'event': 'terminal', 'generation': self.summary['generation'],
                  'run_id': self.summary['run_id'], 'status': self.summary['status'],
                  'exit_code': code, 'path': str(terminal)}), flush=True)
            self.check_time()
            return code
        except BaseException as error:
            if published:
                terminal.unlink()
            compact = {'status': 'failed', 'exit_code': 1,
                       'generation': self.summary['generation'], 'run_id': self.summary['run_id'],
                       'failure': {'type': type(error).__name__, 'message': str(error)[:240]},
                       'scope': 'Publication failed; committed scientific checkpoints retained.'}
            try:
                self.write(terminal, compact, terminal=True)
            except BaseException as fallback_error:
                refusal(self.summary, fallback_error)
            return 1
