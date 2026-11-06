"""Finite inotify observer of a new producer-owned bank lease; no polling."""
from pathlib import Path
import ctypes
import hashlib
import json
import os
import select
import sys
import time
HERE = Path(__file__).resolve().parent
PRODUCER = HERE.parent/'0081'


def wait(after_release, seconds):
    release_ns = Path(after_release).stat().st_mtime_ns
    release_bytes=Path(after_release).read_bytes()
    release_binding=dict(path=str(Path(after_release).resolve()),bytes=len(release_bytes),sha256=hashlib.sha256(release_bytes).hexdigest())
    libc = ctypes.CDLL(None,use_errno=True)
    fd = libc.inotify_init1(os.O_NONBLOCK|os.O_CLOEXEC)
    if fd < 0:
        raise OSError(ctypes.get_errno(),'inotify_init1')
    try:
        for path in (PRODUCER,PRODUCER/'handoff'):
            if libc.inotify_add_watch(fd,os.fsencode(path),0x00000008|0x00000080|0x00000100) < 0:
                raise OSError(ctypes.get_errno(),'inotify_add_watch')
        deadline = time.monotonic()+seconds
        while True:
            request_path,state_path=PRODUCER/'bank-core-request.json',PRODUCER/'bank-core-state.json'
            if request_path.exists() and state_path.exists():
                request_bytes=request_path.read_bytes()
                request=json.loads(request_bytes)
                request_binding=dict(path=str(request_path.resolve()),bytes=len(request_bytes),sha256=hashlib.sha256(request_bytes).hexdigest())
                state=json.loads(state_path.read_bytes())
                if request_path.stat().st_mtime_ns > release_ns and request.get('after_release')==release_binding and state.get('request_binding')==request_binding and request.get('pending') and state.get('native_busy') is False and state.get('cpu')==request.get('cpu') and state.get('pane_id')==request.get('pane_id'):
                    pending=[]
                    for handoff in sorted((PRODUCER/'handoff').iterdir()):
                        if handoff.suffix=='.json' and not (HERE/'intake'/handoff.name).exists():
                            pending.append(str(handoff))
                    if pending:
                        print(json.dumps(dict(event='actual-new-bank-lease-ready',request=request,state=state,pending_handoffs=pending)),flush=True)
                        return 0
            remaining=deadline-time.monotonic()
            if remaining <= 0:
                print(json.dumps(dict(event='finite-lease-observer-deadline',seconds=seconds)),flush=True)
                return 2
            ready,_,_=select.select([fd],[],[],remaining)
            if ready:
                os.read(fd,262144)
    finally:
        os.close(fd)


if __name__=='__main__':
    raise SystemExit(wait(sys.argv[1],float(sys.argv[2])))
