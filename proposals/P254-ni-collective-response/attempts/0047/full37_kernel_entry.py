"""Capture actual owned numerical kernel/manager state before A47 execution."""
import hashlib
import json
import os
from pathlib import Path
import runpy
import subprocess
import sys

record_path,unit,cpu,script,*arguments=sys.argv[1:]
argv=['env','XDG_RUNTIME_DIR=/run/user/1000','DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/1000/bus',
    'systemctl','--user','show',unit,
    '--property=InvocationID,ActiveState,SubState,MainPID,RuntimeMaxUSec,MemoryMax,MemorySwapMax,Restart,CPUAffinity,AllowedCPUs,EffectiveCPUs,WatchdogUSec,ExecMainStartTimestamp,ExecStart']
manager=subprocess.run(argv,capture_output=True)
properties=dict(line.split('=',1) for line in manager.stdout.decode().splitlines() if '=' in line)
record=dict(schema='p254-a47-actual-owned-kernel-start-v1',pid=os.getpid(),python=sys.executable,cwd=os.getcwd(),
    actual_cpu_affinity=sorted(os.sched_getaffinity(0)),actual_proc_status=Path('/proc/self/status').read_text(),
    actual_cgroup=Path('/proc/self/cgroup').read_text(),invocation_id=os.environ.get('INVOCATION_ID'),
    thread_environment={name:os.environ.get(name) for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')},
    executed_script=dict(path=script,bytes=Path(script).stat().st_size,sha256=hashlib.sha256(Path(script).read_bytes()).hexdigest()),
    executed_argv=[script,*arguments],manager=dict(argv=argv,stdout=manager.stdout.decode(),stderr=manager.stderr.decode(),exit_code=manager.returncode),
    optional_interfaces={name:({'status':'available','value':properties[name]} if properties.get(name) else {'status':'unavailable','reason':'actual manager emits no delegated/effective CPU value'}) for name in ('AllowedCPUs','EffectiveCPUs')})
with Path(record_path).open('x') as stream:
    json.dump(record,stream,indent=2,sort_keys=True);stream.write('\n');stream.flush();os.fsync(stream.fileno())
if record['actual_cpu_affinity']!=[int(cpu)] or any(value!='1' for value in record['thread_environment'].values()):
    raise RuntimeError('actual one-thread singleton CPU grant differs')
sys.argv=[script,*arguments];sys.path.insert(0,str(Path(script).resolve().parent))
runpy.run_path(script,run_name='__main__')
