"""Package observed findings/bytes/cost; no build or physical solver replay."""
from pathlib import Path
import json
import os
import resource
import sys
import time
import numpy as np
from production_transport import binding, array_digest, atomic_json, run_recorded
HERE=Path(__file__).resolve().parent
ROOT=Path('/var/tmp/p254-a61-lean-001')
PYTHON='/home/administrator/projects/substrate-framework/.venv/bin/python'


def finish():
    start=time.monotonic()
    if sorted(os.sched_getaffinity(0)) != [9]:raise RuntimeError('wrong finalizer affinity')
    actual=json.loads((ROOT/'actual-material-findings.json').read_text())
    comparisons=[]
    for result in actual['findings']:
        bank=ROOT/result['label']
        saved=json.loads((bank/'receipts'/'slot-000.json').read_text())
        with np.load(result['old_solution']['path'],allow_pickle=False) as prior:
            r=np.load(bank/'ret.npy',mmap_mode='r');l=np.load(bank/'less.npy',mmap_mode='r')
            oldr,oldl=array_digest(prior['ret']),array_digest(prior['less'])
            newr,newl=array_digest(r[0]),array_digest(l[0])
            comparisons.append(dict(label=result['label'],old_ret=oldr,new_ret=newr,old_less=oldl,new_less=newl,canonical_byte_hashes_equal=(oldr==newr and oldl==newl),original_packed_input_sha256=json.loads((Path(result['old_solution']['path']).parent/'manifest.json').read_text())['transport']['sha256'],new_packed_input_sha256=saved['provenance']['input']['sha256'],complete=binding(bank/'complete.json')))
            r._mmap.close();l._mmap.close()
    outcome=json.loads((ROOT/'native-outcome.json').read_text())
    systemd=json.loads((ROOT/'systemd-native.json').read_text())
    expected=json.loads((ROOT/'expected-native-extent-rejection.json').read_text())
    record=dict(schema='p254-a61-complete-handoff-v1',api=binding(HERE/'production-api.json'),contract=binding(HERE/'contract.json'),implementation=[binding(p) for p in sorted(HERE.iterdir()) if p.is_file()],binary=binding(ROOT/'build'/'solve_result'),original_bound_native_inputs=binding(ROOT/'native-bindings.json'),boundary=json.loads((ROOT/'boundary-findings.json').read_text()),actual_material=actual,Nt31=json.loads((ROOT/'Nt31-findings.json').read_text()),canonical_byte_comparisons=comparisons,all_finding_arrays=[binding(p) for p in sorted(ROOT.glob('*-complete-findings.npz'))],all_exit_receipts=[binding(p) for p in sorted(ROOT.rglob('native.json'))]+[binding(ROOT/p) for p in ('configure.json','build.json','exercise.json','expected-native-extent-rejection.json','systemd-native.json','native-outcome.json','terminal.json')],native_cost=dict(cumulative_service_wall_seconds=outcome['wall_seconds'],cumulative_manager_launch_wall_seconds=systemd['wall_seconds'],native_cpu_seconds=outcome['cpu_seconds'],child_peak_rss_kib=outcome['child_maxrss_kib'],owner_peak_rss_kib=outcome['owner_maxrss_kib'],build_count=1,unplanned_failures=0,expected_rejections=[dict(exit=expected['exit'],receipt=binding(ROOT/'expected-native-extent-rejection.json'))],finalizer_wall_before_package_seconds=time.monotonic()-start,limit_seconds=600),profiles=[binding(p) for p in sorted(ROOT.rglob('native.time'))]+[binding(ROOT/'build.time')],manager_memory_warning='systemd reports 256KiB peak despite GNU time measured 723108KiB build and 255276KiB Gamma; this manager peak is not credible application-memory evidence, use per-native GNU-time and getrusage measurements',performance_scope='Gamma new27.876272717s vs paid old26.317892029s: no measured solver speedup; reduced serialized result bytes and final numeric storage, no claim of internal A54 scratch reduction',storage_scope='One active slot transport; no full F/Fcc serialization/decoding, input.npz/output solution.npz duplication removed; existing A54 math/scratch unchanged; all staging released after durable canonical final hashes; complete-slot-only recovery, single writer; no original scientific files removed',resource_scope='CPU9/BLAS1/3GiB/noSwap; protected CPUs0-5/11, RootCPU10 and root view untouched',remaining_physics_scope='Only original material fermion slot000/full114Gamma3node supplied solve was executed here; no all128 replay or missing material31 collision history. Nt31 is explicitly free finite2level zero-selfenergy control.',observed_payload_bytes_before_package=sum(p.stat().st_size for p in ROOT.rglob('*') if p.is_file()))
    atomic_json(HERE/'complete-handoff.json',record)
    atomic_json(ROOT/'handoff-outcome.json',dict(exit=0,wall_seconds=time.monotonic()-start,cpu_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime,handoff=binding(HERE/'complete-handoff.json'),affinity=sorted(os.sched_getaffinity(0))))
    print('A61_HANDOFF_PACKAGED_EXIT_0',flush=True)


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='--native':finish()
    else:
        argv=['env','XDG_RUNTIME_DIR=/run/user/1000','DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/1000/bus','systemd-run','--user','--wait','--pipe','--collect','--unit=p254-a61-handoff-001.service','--property=RuntimeMaxSec=20','--property=MemoryMax=3221225472','--property=MemorySwapMax=0','--property=Restart=no','--property=CPUAffinity=9','env','HERDR_ENV=1','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1','MKL_NUM_THREADS=1','NUMEXPR_NUM_THREADS=1','BLIS_NUM_THREADS=1','PYTHONDONTWRITEBYTECODE=1','taskset','-c','9',PYTHON,'-u',str(Path(__file__).resolve()),'--native']
        code=1
        try:
            run_recorded(argv,ROOT,'systemd-handoff',25);code=0
        finally:
            print('A61_HANDOFF_TERMINAL_EXIT_'+str(code),flush=True)
