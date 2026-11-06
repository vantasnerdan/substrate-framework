"""Mmap persistent active preparation arrays; native payloads are admitted once."""
from collections import OrderedDict
import json
import os
from pathlib import Path
import numpy as np
from substrate_framework.streamed_material_state import file_binding
from reduced_projection_storage import ImmutableItems, _publish


class MappedPreparations:
    def __init__(self, cache):
        self.cache = cache
        self.reader = ImmutableItems(cache.manifest['items'])
        self.root = cache.directory/'prepared-maps'
        self.root.mkdir(exist_ok=True)
        self.recent = OrderedDict()

    def __call__(self, index):
        if index in self.recent:
            self.recent.move_to_end(index)
            return self.recent[index]
        key = 'preparation/%03d' % index
        if key not in self.cache.manifest['items']:
            raise KeyError(key)
        item = self.cache.manifest['items'][key]
        directory = self.root/('%03d' % index); directory.mkdir(exist_ok=True)
        receipt_path = directory/'completed.json'
        if receipt_path.exists():
            receipt = json.loads(receipt_path.read_text())
            if receipt['source_archive'] != item['archive']:
                raise ValueError('mapped preparation source bytes differ')
        else:
            arrays = self.reader.get(key)
            bindings = {}
            for name, value in arrays.items():
                path = directory/(name+'.npy')
                temporary = directory/(name+'.partial.npy')
                with temporary.open('wb') as stream:
                    np.save(stream, value, allow_pickle=False)
                    stream.flush(); os.fsync(stream.fileno())
                os.replace(temporary, path)
                bindings[name] = file_binding(path)
            receipt = dict(schema='p254-mapped-active-preparation-v1',
                           source_archive=item['archive'], arrays=bindings)
            _publish(receipt_path, receipt)
        arrays = {name: np.load(self.reader.authenticate(binding), mmap_mode='r', allow_pickle=False)
                  for name, binding in receipt['arrays'].items()}
        result = arrays, item['metadata']
        self.recent[index] = result
        while len(self.recent) > 2:
            self.recent.popitem(last=False)
        return result
