"""Authenticated NPY descriptors and zero-extraction legacy stored-NPZ views."""
from pathlib import Path
import json
import struct
import sys
import zipfile
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'0061'))
from production_transport import binding, atomic_json

AUXILIARY = ('source_electron_less0', 'control_ret', 'control_less',
             'source_gamma_ret', 'source_gamma_less', 'n')


def authenticate(item):
    path = Path(item['path'])
    if not path.is_absolute() or binding(path) != {k: item[k] for k in ('path', 'bytes', 'sha256')}:
        raise ValueError('source binding differs: '+str(path))
    return path


def descriptor(path):
    item = binding(path)
    a = np.load(path, mmap_mode='r', allow_pickle=False)
    try:
        item.update(shape=list(a.shape), dtype=a.dtype.str)
    finally:
        a._mmap.close()
    return item


def mapped(item):
    a = np.load(authenticate(item), mmap_mode='r', allow_pickle=False)
    if list(a.shape) != item['shape'] or a.dtype.str != item['dtype']:
        a._mmap.close()
        raise ValueError('NPY descriptor layout differs')
    return a


class StoredNPZ:
    """Map original ZIP_STORED NPY members, never extract a history copy.

    Compressed archives are rejected rather than silently expanding all-history
    arrays. The original paid A58/A57 np.savez archives use ZIP_STORED.
    """
    def __init__(self, path):
        self.path = Path(path)
        self.archive = zipfile.ZipFile(self.path)
        self.arrays = {}

    def __getitem__(self, key):
        if key not in self.arrays:
            entry = self.archive.getinfo(key+'.npy')
            if entry.compress_type != zipfile.ZIP_STORED or entry.flag_bits & 1:
                raise ValueError('only original uncompressed unencrypted NPY members are admitted')
            with self.path.open('rb') as stream:
                stream.seek(entry.header_offset)
                header = stream.read(30)
                if header[:4] != b'PK\x03\x04':
                    raise ValueError('invalid ZIP local header')
                name_size, extra_size = struct.unpack_from('<HH', header, 26)
                start = entry.header_offset+30+name_size+extra_size
                stream.seek(start)
                version = np.lib.format.read_magic(stream)
                if version == (1, 0):
                    shape, fortran, dtype = np.lib.format.read_array_header_1_0(stream)
                elif version == (2, 0):
                    shape, fortran, dtype = np.lib.format.read_array_header_2_0(stream)
                else:
                    raise ValueError('unsupported NPY member version')
                offset = stream.tell()
                if dtype.hasobject or offset-start+int(np.prod(shape))*dtype.itemsize != entry.file_size:
                    raise ValueError('unsafe or inconsistent NPY member extent')
            self.arrays[key] = np.memmap(self.path, mode='r', dtype=dtype,
                shape=shape, offset=offset, order='F' if fortran else 'C')
        return self.arrays[key]

    def close(self):
        for a in self.arrays.values():
            a._mmap.close()
        self.arrays.clear()
        self.archive.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def finite_blocks(a):
    if a.ndim < 2:
        return bool(np.isfinite(a).all())
    return all(np.isfinite(a[index]).all() for index in np.ndindex(a.shape[:-2]))


def save_array(path, source):
    """One final NPY; matrixwise copy without an all-history temporary."""
    target = np.lib.format.open_memmap(path, mode='w+', dtype=source.dtype, shape=source.shape)
    try:
        if source.ndim < 2:
            target[...] = source
        else:
            for index in np.ndindex(source.shape[:-2]):
                target[index] = source[index]
        target.flush()
    finally:
        target._mmap.close()
    return descriptor(path)
