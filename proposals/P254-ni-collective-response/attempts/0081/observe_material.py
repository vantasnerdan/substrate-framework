"""Observe actual material terminals, errors, inactivity, or an absolute decision deadline."""
from pathlib import Path
import ctypes
from datetime import datetime, timezone
import json
import os
import select
import sys
import time

root=Path(__file__).resolve().parent
deadline=datetime.fromisoformat(sys.argv[1]).timestamp()
campaign=len(sys.argv)>2 and sys.argv[2]=='campaign'
last_material=time.time()
def receipt_generations():
    paths={p for pattern in ('settled-*.json','failed-*.json','blocked-*.json','launch-blocked-*.json','dispatch-blocked.json','backpressure.json') for p in root.glob(pattern)}
    return {p:p.stat().st_mtime_ns for p in paths}

initial=receipt_generations()
libc=ctypes.CDLL(None,use_errno=True)
fd=libc.inotify_init1(os.O_NONBLOCK|os.O_CLOEXEC)
if fd<0: raise OSError(ctypes.get_errno(),'inotify_init1')
try:
    if libc.inotify_add_watch(fd,os.fsencode(root),0x80|0x100|0x8)<0:
        raise OSError(ctypes.get_errno(),'inotify_add_watch')
    while True:
        current=receipt_generations()
        fresh=sorted(p for p,generation in current.items() if initial.get(p)!=generation)
        if fresh:
            print(json.dumps(dict(event='material-terminal',receipts=[str(p) for p in fresh])),flush=True)
            if not campaign or any(p.name.startswith(('failed-','blocked-','launch-blocked-','dispatch-blocked')) or p.name=='backpressure.json' for p in fresh):
                break
            initial=current
            last_material=time.time()
        ledger_path=root/'actual-cost-ledger.json'
        if campaign and ledger_path.exists():
            ledger=json.loads(ledger_path.read_text())
            if not ledger['active'] and not any(ledger['remaining_queues'].values()):
                print(json.dumps(dict(event='producer-all-dispatched-cases-terminal',ledger=str(ledger_path))),flush=True)
                break
        remaining=min(deadline,last_material+2700 if campaign else deadline)-time.time()
        if remaining<=0:
            print(json.dumps(dict(event='no-required-component-45min' if campaign and time.time()<deadline else 'decision-deadline',actual_utc=datetime.now(timezone.utc).isoformat())),flush=True)
            break
        select.select([fd],[],[],remaining)
        try: os.read(fd,65536)
        except BlockingIOError: pass
finally:
    os.close(fd)
