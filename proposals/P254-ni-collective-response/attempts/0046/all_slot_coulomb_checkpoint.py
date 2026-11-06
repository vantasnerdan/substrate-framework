"""Specific two-generation all-slot/q-prefix accumulator, not a generic cache.

An archive and its binding sidecar are committed to an inactive numeric slot.
A killed write leaves the other generation intact; consumers authenticate bytes
before decoding. At most two committed archives and one temporary archive.
"""
import hashlib
import json
import os
from pathlib import Path
import numpy as np
from all_slot_coulomb_contractions import KEYS


def binding(path):
    path = Path(path).resolve()
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024**2),b''):
            h.update(chunk)
    return dict(path=str(path),bytes=path.stat().st_size,sha256=h.hexdigest())


def authenticate(record):
    # Historical original receipts encode exact byte counts as decimal strings.
    # Keep their original schema in provenance; compare losslessly, never float.
    count = record['bytes']
    if not (type(count) is int or isinstance(count,str) and count.isascii() and count.isdecimal()):
        raise ValueError('byte count must be an exact integer or decimal integer string')
    actual = binding(record['path'])
    if (actual['path'] != record['path']
            or actual['sha256'] != record['sha256'] or actual['bytes'] != int(count)):
        raise ValueError('original-byte authentication failed: '+record['path'])
    return Path(record['path'])


def atomic_json(path,value):
    path = Path(path)
    stage = path.with_suffix(path.suffix+'.stage')
    with stage.open('w') as f:
        json.dump(value,f,sort_keys=True,indent=2,allow_nan=False)
        f.flush(); os.fsync(f.fileno())
    os.replace(stage,path)


class AllSlotPrefix:
    def __init__(self,directory,specification,shape,nq):
        self.directory = Path(directory).resolve()
        self.directory.mkdir(parents=True,exist_ok=True)
        self.specification = specification
        self.shape,self.nq = tuple(shape),int(nq)
        self.owned_q = list(specification.get('outer_q_partition',range(self.nq)))
        self.global_nq = int(specification.get('global_nq',self.nq))
        if (len(self.owned_q) != self.nq or self.owned_q != sorted(set(self.owned_q))
                or not self.owned_q or min(self.owned_q)<0 or max(self.owned_q)>=self.global_nq):
            raise ValueError('exact outer-q partition/count/order differs')
        self.recoveries = []
        spec = self.directory/'specification.json'
        if spec.exists():
            if json.loads(spec.read_bytes()) != specification:
                raise ValueError('prefix source/input/time/batch binding differs')
        else:
            atomic_json(spec,specification)

    def _decode(self,record):
        if (record['schema'] != 'p254-a46-all-slot-q-prefix-v1'
                or record['specification'] != binding(self.directory/'specification.json')
                or record['external_slots'] != list(range(self.shape[0]))
                or record['shape'] != list(self.shape)
                or not 1 <= record['completed_q'] <= self.nq
                or record['owned_q'] != self.owned_q[:record['completed_q']]
                or record['ordered_pairs'] != self.shape[1]*self.shape[2]
                or record['exchange_inner_h_count'] != self.global_nq):
            raise ValueError('numeric prefix ownership/counts differ')
        with np.load(authenticate(record['archive']),allow_pickle=False) as a:
            if set(a.files) != set(KEYS):
                raise ValueError('prefix channels differ')
            out = {k:a[k] for k in KEYS}
        if any(x.shape != self.shape or x.dtype != np.dtype('complex128')
               or not np.isfinite(x).all() for x in out.values()):
            raise ValueError('prefix numeric axes/type/finiteness differ')
        return out

    def restore(self):
        candidates = []
        for slot in range(2):
            path = self.directory/f'prefix-{slot}.json'
            if not path.exists():
                continue
            try:
                record = json.loads(path.read_bytes())
                # Byte-authenticate before retaining a potentially largest generation.
                authenticate(record['archive'])
                candidates.append(record)
            except (OSError,ValueError,KeyError) as e:
                self.recoveries.append(dict(sidecar=str(path),rejected=str(e)))
        for record in sorted(candidates,key=lambda r:r['completed_q'],reverse=True):
            try:
                return self._decode(record),record['completed_q']
            except (OSError,ValueError,KeyError) as e:
                self.recoveries.append(dict(archive=record.get('archive'),rejected=str(e)))
        if candidates or self.recoveries:
            raise ValueError('no authenticated numeric generation survives')
        return None,0

    def commit(self,arrays,completed_q):
        if set(arrays) != set(KEYS) or any(x.shape != self.shape or not np.isfinite(x).all()
                                          for x in arrays.values()):
            raise ValueError('uncommittable numeric prefix')
        slot = completed_q%2
        path = self.directory/f'prefix-{slot}.npz'
        stage = self.directory/'prefix-stage.npz'
        with stage.open('wb') as f:
            np.savez(f,**arrays)
            f.flush(); os.fsync(f.fileno())
        os.replace(stage,path)
        record = dict(schema='p254-a46-all-slot-q-prefix-v1',
            specification=binding(self.directory/'specification.json'),archive=binding(path),
            shape=list(self.shape),external_slots=list(range(self.shape[0])),
            completed_q=completed_q,owned_q=self.owned_q[:completed_q],
            ordered_pairs=self.shape[1]*self.shape[2],exchange_inner_h_count=self.global_nq,
            external_batches=self.specification['external_batches'])
        atomic_json(self.directory/f'prefix-{slot}.json',record)


def publish_result(directory,specification,arrays,counts,recoveries):
    directory = Path(directory)
    path = directory/'all-slots-all-pairs.npz'
    stage = directory/'all-slots-all-pairs.stage.npz'
    with stage.open('wb') as f:
        np.savez(f,**arrays)
        f.flush(); os.fsync(f.fileno())
    os.replace(stage,path)
    record = dict(schema='p254-a46-all-slot-result-v1',status='complete',
        contribution_scope=specification.get('contribution_scope','FULL_ALL_Q'),
        specification=specification,archive=binding(path),counts=counts,
        recoveries=recoveries,norms={k:float(np.linalg.norm(x)) for k,x in arrays.items()})
    atomic_json(directory/'result.json',record)
    return consume_result(directory,specification)


def consume_result(directory,specification):
    """Decode final saved values, never authenticate metadata without using numbers."""
    path = Path(directory)/'result.json'
    if not path.exists():
        return None
    record = json.loads(path.read_bytes())
    if (record['schema'] != 'p254-a46-all-slot-result-v1' or record['status'] != 'complete'
            or record['specification'] != specification
            or record['counts']['completed_q'] != specification['nq']
            or record['contribution_scope'] != specification.get('contribution_scope','FULL_ALL_Q')):
        raise ValueError('complete source/input/count identity differs')
    for item in specification['authenticated_sources_and_inputs']:
        authenticate(item)
    with np.load(authenticate(record['archive']),allow_pickle=False) as a:
        arrays = {k:a[k] for k in a.files}
    if (set(arrays) != set(KEYS) or any(x.shape != tuple(specification['shape'])
            or x.dtype != np.dtype('complex128') or not np.isfinite(x).all() for x in arrays.values())):
        raise ValueError('saved full-slot output axes/type/finiteness differ')
    norms = {k:float(np.linalg.norm(x)) for k,x in arrays.items()}
    if norms != record['norms']:
        raise ValueError('decoded saved norms differ')
    return dict(record=record,arrays=arrays,decoded_norms=norms)
