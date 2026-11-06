"""Durably capture exactly one requested systemd Native invocation."""
from pathlib import Path
import json
import os
import subprocess
import time
HERE=Path(__file__).resolve().parent
contract=json.loads((HERE/'implementation-contract-001.json').read_bytes())
python=contract['common_runtime']['python']
output=Path('/var/tmp/p254-a45-moving-reference-001')
receipts=Path('/var/tmp/p254-a45-native-receipts-001');receipts.mkdir(exist_ok=False)
unit='p254-a45-moving-wave-001.service'
argv=contract['common_runtime']['manager_env']+['systemd-run','--user','--wait','--pipe',
    '--unit='+unit,'--property=RuntimeMaxSec=1800','--property=MemoryMax=4294967296',
    '--property=MemorySwapMax=0','--property=Restart=no','--property=CPUAffinity=7',
    'env',contract['common_runtime']['python_environment_assignment']]+[
    key+'='+value for key,value in contract['common_runtime']['threads'].items()]+[
    'taskset','-c','7',python,'-u',str(HERE/'native_wave.py'),str(output)]
def save(name,value):
    with (receipts/name).open('x') as stream:
        json.dump(value,stream,indent=2,sort_keys=True);stream.write('\n');stream.flush();os.fsync(stream.fileno())
start=time.time();save('actual-native-argv.json',dict(unit=unit,argv=argv,start_unix=start))
with (receipts/'native-stdout.txt').open('xb') as stdout,(receipts/'native-stderr.txt').open('xb') as stderr:
    result=subprocess.run(argv,stdout=stdout,stderr=stderr,check=False)
manager_argv=contract['common_runtime']['manager_env']+['systemctl','--user','show',unit,
    '--property=Id,InvocationID,Result,ExecMainCode,ExecMainStatus,ExecMainStartTimestamp,ExecMainExitTimestamp,CPUUsageNSec,MemoryPeak,MemorySwapPeak,MemoryMax,MemorySwapMax,Restart,RuntimeMaxUSec,ActiveState,SubState']
manager=subprocess.run(manager_argv,capture_output=True,check=False)
save('actual-manager-terminal.json',dict(argv=manager_argv,exitcode=manager.returncode,
    stdout=manager.stdout.decode(),stderr=manager.stderr.decode()))
save('actual-native-exit.json',dict(argv=argv,exitcode=result.returncode,start_unix=start,
    end_unix=time.time(),elapsed_seconds=time.time()-start))
print(json.dumps(dict(unit=unit,exitcode=result.returncode,receipts=str(receipts),output=str(output))),flush=True)
