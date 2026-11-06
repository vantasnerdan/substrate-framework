"""Immutable JSON receipts and postexec Linux worker observations."""
import json
import os
from pathlib import Path
import time
from source_graph import digest


def save(path,value):
    with Path(path).open('x') as stream:
        json.dump(value,stream,indent=2,sort_keys=True,allow_nan=False)
        stream.write('\n');stream.flush();os.fsync(stream.fileno())


def observe(pid):
    proc=Path('/proc')/str(pid)
    stat=(proc/'stat').read_text(); fields=stat[stat.rfind(')')+2:].split()
    cgroup=(proc/'cgroup').read_text().strip().split('::')[-1]
    kernel=Path('/sys/fs/cgroup')/cgroup.lstrip('/')
    limits={key:(kernel/key).read_text().strip() for key in
            ('memory.max','memory.swap.max','memory.current','memory.peak','memory.events')}
    cpuset=kernel/'cpuset.cpus.effective'
    limits['cpuset.cpus.effective']=cpuset.read_text().strip() if cpuset.exists() else None
    limits['cpu_proof']='actual sched_getaffinity and /proc status; cpuset controller may not be delegated'
    executable=str((proc/'exe').resolve())
    return {'pid':pid,'start_ticks':fields[19],'boot_id':Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
            'observed_time_ns':str(time.time_ns()),'cmdline':(proc/'cmdline').read_bytes().decode().split('\0')[:-1],
            'executable':executable,'executable_sha256':digest(executable),
            'affinity':sorted(os.sched_getaffinity(pid)),'cgroup':cgroup,'live_kernel':limits,
            'status':(proc/'status').read_text(),'limits':(proc/'limits').read_text()}


def checkpoint(directory,block,content):
    path=directory/(block+'-checkpoint.json');save(path,content)
    state={'last_complete_block':block,'checkpoint':str(path),'sha256':digest(path)}
    temporary=directory/'checkpoint-next.json';save(temporary,state)
    os.replace(temporary,directory/'checkpoint.json')
