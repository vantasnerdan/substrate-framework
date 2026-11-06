"""One bounded paid component invocation; original argv/stdout/stderr/all exits."""
from pathlib import Path
import argparse
import json
import os
import subprocess
import time

HERE = Path(__file__).resolve().parent


def launch(component):
    contract = json.loads((HERE.parent/'0045/implementation-contract-001.json').read_bytes())['common_runtime']
    cpu, cap = {'ee':(6,420),'cubic':(8,90)}[component]
    output = Path('/var/tmp')/('p254-a66-'+component+'-first-column-001')
    receipts = HERE/(component+'-receipts-001')
    receipts.mkdir(exist_ok=False)
    unit = 'p254-a66-'+component+'-first-column-001.service'
    context = {key:value for key,value in os.environ.items() if key.startswith('HERDR_')}
    if os.environ.get('HERDR_ENV') != '1' or not all(context.get(key) for key in ('HERDR_PANE_ID','HERDR_TAB_ID','HERDR_WORKSPACE_ID')):
        raise ValueError('verified actual owner Herdr context required')
    argv = contract['manager_env']+['systemd-run','--user','--wait','--pipe','--unit='+unit,
        '--property=RuntimeMaxSec='+str(cap),'--property=MemoryMax=3221225472',
        '--property=MemorySwapMax=0','--property=Restart=no','--property=CPUAffinity='+str(cpu),
        'env',contract['python_environment_assignment'],
        *[key+'='+value for key,value in contract['threads'].items()],
        *[key+'='+value for key,value in context.items()], 'PYTHONDONTWRITEBYTECODE=1',
        '/usr/bin/time','-v','-o',str(receipts/'time.txt'),contract['python'],'-u',
        str(HERE/('produce_'+component+'_column.py')),'--output',str(output)]
    def save(name,value):
        (receipts/name).write_text(json.dumps(value,indent=2)+'\n')
    started = time.monotonic()
    save('launch.json',dict(argv=argv,unit=unit,cap_seconds=cap,
        actual_start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),context=context))
    print('A66_'+component.upper()+'_NATIVE_START '+unit,flush=True)
    code = None
    try:
        with (receipts/'stdout').open('xb') as stdout, (receipts/'stderr').open('xb') as stderr:
            process = subprocess.Popen(argv,stdout=subprocess.PIPE,stderr=stderr)
            for line in process.stdout:
                stdout.write(line)
                stdout.flush()
                print(line.decode(errors='replace'),end='',flush=True)
            code = process.wait()
    finally:
        manager = contract['manager_env']+['systemctl','--user','show',unit,
            '--property=Id,InvocationID,Result,ExecMainCode,ExecMainStatus,CPUUsageNSec,MemoryPeak,MemorySwapPeak,MemoryMax,MemorySwapMax,Restart,RuntimeMaxUSec,ActiveState,SubState']
        terminal = subprocess.run(manager,capture_output=True,check=False)
        save('manager-terminal.json',dict(argv=manager,returncode=terminal.returncode,
            stdout=terminal.stdout.decode(),stderr=terminal.stderr.decode()))
        save('exit.json',dict(argv=argv,returncode=code,elapsed_seconds=time.monotonic()-started))
        print('A66_'+component.upper()+'_TERMINAL_'+str(code),flush=True)
    return code if code is not None else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('component',choices=('ee','cubic'))
    raise SystemExit(launch(parser.parse_args().component))
