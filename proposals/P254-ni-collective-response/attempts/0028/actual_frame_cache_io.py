"""Owned frame-cache filesystem, authentication and durable publication helpers."""
from contextlib import contextmanager
import fcntl
from hashlib import sha256
import json
import os
from pathlib import Path
import stat

import numpy as np
from substrate_framework.partitioned_joint_slots import digest
from substrate_framework.partitioned_q_storage import durable_publish, sync_directory, sync_file
from substrate_framework.streamed_material_state import authenticate, file_binding

STORAGE = ('Hardlinks, not copies; each surviving cache/adoption alias pins physical '
           'storage after producer retirement. Logical retired bytes are not '
           'physical reclaimed bytes while any alias survives.')


class _JsonWriter:
    """Metadata-only writer for the existing durable_publish protocol."""
    def json(self, path, value):
        data = (json.dumps(value, indent=2, allow_nan=False)+'\n').encode()
        with Path(path).open('xb') as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        return file_binding(path)


def _publish(path, value):
    _safe(path, missing=True)
    return durable_publish(path, value, _JsonWriter())


def _safe(path, *, missing=False):
    """Reject symlinks before resolving, including every existing ancestor."""
    path = Path(path)
    if '..' in path.parts:
        raise ValueError('parent traversal is not an admitted path')
    path = Path(os.path.abspath(path))
    parts = list(reversed(path.parents))+[path]
    for index, part in enumerate(parts):
        try:
            mode = part.lstat().st_mode
        except FileNotFoundError:
            if missing:
                return path
            raise
        if stat.S_ISLNK(mode):
            raise ValueError('symlink path refused: '+str(part))
        if index < len(parts)-1 and not stat.S_ISDIR(mode):
            raise ValueError('non-directory ancestor: '+str(part))
        if index == len(parts)-1 and not (stat.S_ISDIR(mode) or stat.S_ISREG(mode)):
            raise ValueError('non-regular path refused: '+str(part))
    return path


def _snapshot(path):
    path = _safe(path)
    data = path.read_bytes()
    value = json.loads(data, parse_constant=lambda x: (_ for _ in ()).throw(
        ValueError('nonfinite JSON: '+x)))
    digest(value)  # Match QWork's canonical finite JSON convention.
    binding = dict(path=str(path), bytes=len(data), sha256=sha256(data).hexdigest())
    authenticate(binding)
    return value, binding


def _json(path):
    return _snapshot(path)[0]


def _same(left, right):
    return digest(left) == digest(right)


def _bound_json(binding):
    value, actual = _snapshot(binding['path'])
    if actual['sha256'] != binding['sha256'] or actual['bytes'] != int(binding['bytes']):
        raise ValueError('metadata differs from its original expected SHA/size')
    return value


def _directory(path, *, new=False):
    path = _safe(path, missing=True)
    if new:
        # Parents must already exist: never create an unowned tree implicitly.
        _safe(path.parent)
        path.mkdir(mode=0o700, exist_ok=False)
        sync_directory(path.parent)
    info = path.stat()
    if (not stat.S_ISDIR(info.st_mode) or info.st_uid != os.geteuid()
            or info.st_mode & 0o077):
        raise ValueError('owned private directory (mode 0700) required: '+str(path))
    return path


@contextmanager
def _locked(path):
    path = _directory(path)
    fd = os.open(path, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX)
        yield path
    finally:
        os.close(fd)


def _mkdir(path):
    path = _safe(path, missing=True)
    try:
        path.mkdir(mode=0o700)
        sync_directory(path.parent)
    except FileExistsError:
        _directory(path)
    return path


def _checked_binding(binding):
    if (set(binding) != {'path', 'bytes', 'sha256'} or type(binding['bytes']) not in (int, str)
            or int(binding['bytes']) <= 0 or not isinstance(binding['sha256'], str)
            or len(binding['sha256']) != 64
            or any(c not in '0123456789abcdef' for c in binding['sha256'])):
        raise ValueError('explicit regular file SHA256/size binding required')
    path = _safe(binding['path'])
    if str(path) != binding['path'] or not path.is_file():
        raise ValueError('absolute canonical regular file binding required')
    return binding


def _frame(binding, expected, shape):
    _checked_binding(binding)
    if binding['sha256'] != expected['sha256'] or binding['bytes'] != expected['bytes']:
        raise ValueError('alias differs from ORIGINAL expected frame SHA/size')
    authenticate(binding)
    value = np.load(binding['path'], mmap_mode='r', allow_pickle=False)
    try:
        if list(value.shape) != shape or value.dtype != np.complex128:
            raise ValueError('original full complex128 frame shape differs')
    finally:
        if isinstance(value, np.memmap):
            value._mmap.close()


def _link(source, target, expected, shape):
    source, target = _safe(source), _safe(target, missing=True)
    if target.exists():
        # Only interrupted writes in this explicitly owned item may resume.
        # Same content alone is insufficient: no copy substitution is allowed.
        if not os.path.samefile(source, target):
            raise ValueError('existing target is not the original hardlink')
    else:
        os.link(source, target, follow_symlinks=False)
    binding = dict(expected, path=str(target))
    _frame(binding, expected, shape)
    sync_file(target)
    sync_directory(target.parent)
    return binding


def _disjoint(left, right):
    if left == right or left.is_relative_to(right) or right.is_relative_to(left):
        raise ValueError('source/cache/target ownership trees must be distinct')


def _storage(paths):
    inodes, logical = {}, 0
    for path in paths:
        info = _safe(path).stat()
        logical += info.st_size
        inodes[(info.st_dev, info.st_ino)] = info.st_size
    return dict(logical_payload_link_bytes=logical,
                unique_pinned_inode_bytes=sum(inodes.values()),
                accounting='Payload bytes only, not metadata or newly allocated physical bytes. '+STORAGE)
