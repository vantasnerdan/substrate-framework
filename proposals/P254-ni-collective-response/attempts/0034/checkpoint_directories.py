"""Publish checkpoint directories and their containing entries durably."""
import os
from pathlib import Path


def sync_directory(path):
    """Synchronize directory entries without changing existing file contents."""
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def make_directory(path, *, exist_ok=False):
    """Create missing parents and sync each newly published directory entry."""
    path = Path(path)
    try:
        path.mkdir()
    except FileNotFoundError:
        make_directory(path.parent, exist_ok=True)
        path.mkdir(exist_ok=exist_ok)
    except FileExistsError:
        if not exist_ok or not path.is_dir():
            raise
    # Sync the directory itself as well as the parent which owns its name.
    # Recursive creation has already persisted every newly created ancestor.
    sync_directory(path)
    sync_directory(path.parent)
    return path
