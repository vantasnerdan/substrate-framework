"""Finite event-driven view of original Native stage bytes; no scientific replay."""
import ctypes
import os
from pathlib import Path
import select
import sys
import time
roots=[Path('/var/tmp/p254-a45-moving-reference-001'),Path('/var/tmp/p254-a45-native-receipts-001')]
libc=ctypes.CDLL(None,use_errno=True)
fd=libc.inotify_init1(os.O_NONBLOCK|os.O_CLOEXEC)
if fd<0:raise OSError(ctypes.get_errno(),'inotify_init1')
for root in roots:
    if libc.inotify_add_watch(fd,os.fsencode(root),0x00000002|0x00000008|0x00000100|0x00000080)<0:
        raise OSError(ctypes.get_errno(),'inotify_add_watch '+str(root))
offsets={};last=None;deadline=time.monotonic()+1800
print('A45 GENUINE MOVING REFERENCE | original stage stdout/stderr | p254-a45-moving-wave-001.service',flush=True)
def display():
    global last
    for root in roots:
        for path in sorted(root.glob('*-std*.txt')):
            offset=offsets.get(path,0)
            with path.open('rb') as stream:
                stream.seek(offset);data=stream.read();offsets[path]=stream.tell()
            if data:
                if last!=path:
                    print('\n--- ORIGINAL '+str(path)+' ---',flush=True);last=path
                sys.stdout.buffer.write(data);sys.stdout.buffer.flush()
try:
    while True:
        display()
        terminal=roots[1]/'actual-native-exit.json'
        if terminal.exists():
            print('\n--- ORIGINAL NATIVE TERMINAL RECEIPT ---\n'+terminal.read_text(),flush=True)
            print('[A45-view-terminal] original unit terminal, interactive pane retained',flush=True)
            break
        remaining=deadline-time.monotonic()
        if remaining<=0:
            print('[A45-view-timeout] finite observer ended; Native is not declared failed',flush=True);break
        ready,_,_=select.select([fd],[],[],remaining)
        if ready:
            try:os.read(fd,65536)
            except BlockingIOError:pass
finally:os.close(fd)
