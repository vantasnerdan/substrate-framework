"""Supply only path-and-SHA-bound byte sequences to scientific consumers.

Validation and parsing/execution share the same in-memory bytes. Original paths
are identities, never permission to reopen an unverified live file.
"""
import hashlib
import json
from pathlib import Path
import sys
import types


class SourceBindingError(ValueError):
    """A requested source is absent, ambiguous, or differs from its binding."""


def read_bound_bytes(path,expected):
    raw=Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=expected:
        raise SourceBindingError('Source bytes differ from registered SHA: '+str(path))
    return raw


def close_source_records(references,snapshots):
    """Close original identities onto independently verified physical objects."""
    verified=[]
    for original,expected in sorted(references):
        actual=str(snapshots.get((original,expected),original))
        raw=read_bound_bytes(actual,expected)
        verified.append({'path':original,'sha256':expected,'consumed_path':actual,'bytes':len(raw)})
    return verified


class VerifiedSources:
    """Resolve a registered identity and recheck the actual bytes on every read."""

    def __init__(self,records):
        self.records={}
        self.consumed=[]
        for record in records:
            key=(record['path'],record['sha256'])
            if key in self.records and self.records[key]!=record:
                raise SourceBindingError('Conflicting source resolution: '+str(key))
            self.records[key]=dict(record)

    def record(self,path,expected=None):
        path=str(path)
        if expected is None:
            candidates=[row for (original,_),row in self.records.items() if original==path]
            if len(candidates)!=1:
                raise SourceBindingError('Source identity absent or SHA required: '+path)
            return dict(candidates[0])
        try:
            return dict(self.records[(path,expected)])
        except KeyError as error:
            raise SourceBindingError('Unregistered path/SHA pair: '+path+' '+expected) from error

    def read_bytes(self,path,expected=None):
        record=self.record(path,expected)
        raw=read_bound_bytes(record['consumed_path'],record['sha256'])
        if len(raw)!=record['bytes']:
            raise SourceBindingError('Registered byte count differs: '+str(path))
        self.consumed.append(record)
        return raw

    def text(self,path,expected=None):
        return self.read_bytes(path,expected).decode('utf-8')

    def json(self,path,expected=None):
        return json.loads(self.read_bytes(path,expected))

    def matches(self,path,expected):
        try:
            self.read_bytes(path,expected)
        except SourceBindingError:
            return False
        return True

    def module(self,path,expected,name):
        raw=self.read_bytes(path,expected)
        module=types.ModuleType(name)
        module.__file__=self.record(path,expected)['consumed_path']
        module.__package__=name.rpartition('.')[0]
        # dataclasses and ordinary subsequent imports must see this exact module.
        sys.modules[name]=module
        exec(compile(raw,str(path),'exec'),module.__dict__)
        return module


def worker_sources(manifest,directory):
    own=[{'path':path,'sha256':sha,'consumed_path':str(directory/'oracle-source-snapshot'/sha),
          'bytes':(directory/'oracle-source-snapshot'/sha).stat().st_size}
         for path,sha in manifest['oracle_sources'].items()]
    return VerifiedSources(manifest['upstream_verified']+own)
