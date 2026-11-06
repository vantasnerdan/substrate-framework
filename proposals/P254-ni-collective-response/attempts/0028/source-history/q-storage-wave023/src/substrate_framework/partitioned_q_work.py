"""Pinned q-edge/frame work, content-addressed receipts and complete barriers.

Inputs are content identities, never a mutable scratch filename. On restart the
caller regenerates the same RK4 inputs and authenticates completed edge outputs.
No worker computes global K or publishes a partially reduced response.
"""
import json
import multiprocessing as mp
import os
from pathlib import Path
import resource
from time import monotonic, process_time
from uuid import uuid4

import numpy as np

from .partitioned_joint_slots import digest, publish, SlotExecutor
from .streamed_material_state import file_binding


class QWork:
    verify = SlotExecutor.verify

    def __init__(self, directory, budget, cpus, identity):
        if (not cpus or len(set(cpus)) != len(cpus)
                or any(type(c) is not int or c not in range(6, 11) for c in cpus)):
            raise ValueError('explicit unique allocated CPUs6-10 required')
        self.directory, self.budget = Path(directory).resolve(), budget
        self.cpus, self.identity = tuple(cpus), identity
        self.directory.mkdir(parents=True, exist_ok=True)
        self._verified = {}
        self.last = None

    def request(self, operation, identity, count):
        request = dict(schema='partitioned-q-work-v1', source_identity=self.identity,
                       operation=operation, input=identity, count=count)
        key = digest(request)
        directory = self.directory/key
        directory.mkdir(exist_ok=True)
        path = directory/'request.json'
        if path.exists():
            if json.loads(path.read_text()) != json.loads(json.dumps(request)):
                raise ValueError('modified q scientific request')
        else:
            publish(path, request, self.budget)
        return key, directory

    def receipt(self, directory, key, index):
        path = directory/('item-%04d' % index)/'complete.json'
        if not path.exists():
            return None
        record = json.loads(path.read_text())
        if (record.get('complete') is not True or record.get('request') != key
                or record.get('index') != index):
            raise ValueError('modified completed q work receipt')
        for binding in record['arrays'].values():
            self.verify(binding)
        return record

    def collect(self, directory, key, count):
        records = [self.receipt(directory, key, i) for i in range(count)]
        if any(r is None for r in records):
            raise ValueError('incomplete q barrier: global feedback forbidden')
        return records

    def _worker(self, cpu, indices, function, directory, key):
        os.sched_setaffinity(0, {cpu})
        from threadpoolctl import threadpool_limits
        with threadpool_limits(limits=1):
            for index in indices:
                self.budget.check()
                target = directory/('item-%04d' % index)
                target.mkdir(exist_ok=True)
                attempt = target/('attempt-'+uuid4().hex)
                attempt.mkdir()
                started, cpu_started = monotonic(), process_time()
                actions = self.budget.local_actions
                values = function(index)
                if any(not np.isfinite(v).all() for v in values.values()):
                    raise ValueError('nonfinite q work cannot be completed')
                arrays = {name: self.budget.array(attempt/(name+'.npy'), value)
                          for name, value in values.items()}
                publish(target/'complete.json', dict(complete=True, request=key, index=index,
                    arrays=arrays, cpu=cpu, affinity=sorted(os.sched_getaffinity(0)), pid=os.getpid(),
                    blas_threads=1, cpu_seconds=process_time()-cpu_started,
                    elapsed_seconds=monotonic()-started,
                    hamiltonian_actions=self.budget.local_actions-actions,
                    process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024), self.budget)
                del values, arrays

    def execute(self, operation, identity, count, function):
        key, directory = self.request(operation, identity, count)
        self.last = key, directory, count
        missing = [i for i in range(count) if self.receipt(directory, key, i) is None]
        if missing:
            ctx = mp.get_context('fork')
            cpus = self.cpus[:len(missing)]
            # Deterministic disjoint assignment guarantees every allocated worker
            # receives real work, including tiny finite-system oracle workloads.
            workers = [ctx.Process(target=self._worker, args=(cpu, missing[j::len(cpus)],
                       function, directory, key)) for j, cpu in enumerate(cpus)]
            journal = directory/('execution-'+uuid4().hex+'.json')
            self.budget.reserve(16384)
            publish(journal, dict(status='allocated', cpus=list(cpus), items=missing,
                    budget=self.budget.snapshot()), self.budget)
            failure = None
            try:
                self.budget.check_processes(os.getpid(), [])
                for worker in workers:
                    worker.start()
                while any(w.is_alive() for w in workers):
                    self.budget.check_processes(os.getpid(), [w.pid for w in workers if w.is_alive()])
                    if any(w.exitcode not in (None, 0) for w in workers):
                        raise RuntimeError('q worker failed; immutable completed work retained')
                    for worker in workers:
                        worker.join(timeout=.1/len(workers))
                if any(w.exitcode != 0 for w in workers):
                    raise RuntimeError('q worker failed; immutable completed work retained')
            except BaseException as error:
                failure = type(error).__name__+': '+str(error)
                raise
            finally:
                for worker in workers:
                    if worker.is_alive():
                        worker.terminate()
                for worker in workers:
                    if worker.pid is not None:
                        worker.join()
                value = dict(status='failed' if failure else 'workers-finished', error=failure,
                    allocation=file_binding(journal), budget=self.budget.snapshot(),
                    exitcodes=[w.exitcode for w in workers])
                data = (json.dumps(value, indent=2)+'\n').encode()
                if len(data) > 16384:
                    raise RuntimeError('q terminal receipt reservation exceeded')
                with journal.with_suffix('.terminal.json').open('xb') as stream:
                    stream.write(data)
        return self.collect(directory, key, count)


def arrays(record):
    return {name: np.load(binding['path'], mmap_mode='r', allow_pickle=False)
            for name, binding in record['arrays'].items()}
