"""Run exactly three disjoint EE workers; no further work before their exits."""
from pathlib import Path
import subprocess,sys,time,json
HERE=Path(__file__).resolve().parent
PYTHON='/home/administrator/projects/substrate-framework/.venv/bin/python'
def main():
    jobs=[]
    for name,cpu,left,right in (('ee-a-001',6,0,21),('ee-b-001',7,22,42),('ee-c-001',8,43,63)):
        argv=[PYTHON,str(HERE/'launch_component.py'),'ee',name,str(cpu),'1800',str(left),str(right)]
        jobs.append((name,subprocess.Popen(argv)))
    codes={name:p.wait() for name,p in jobs}
    print(json.dumps(codes),flush=True);print('A55_EE_WAVE_001_TERMINAL',flush=True)
    return int(any(codes.values()))
if __name__=='__main__':sys.exit(main())
