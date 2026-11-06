"""Run an admitted dressed map with explicit campaign-level cost context."""
from pathlib import Path
import argparse
import concurrent.futures
import hashlib
import json
import sys
import subprocess

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'0056'))
from dressed_reference import binding, authenticate, save
from finalize_functional_iteration import finalize


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('dataset', 'baseline-catalog', 'output', 'ledger', 'carried-oracle'):
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    prior = json.loads(a.ledger.read_bytes())
    authenticate(prior['previous_catalog'])
    if binding(a.dataset) != prior['input_dataset']:
        raise ValueError('allocation and actual dressed input differ')
    for field, limit, additional in (
            ('native_cumulative_seconds', 'native_limit_seconds', 4710),
            ('ee_cumulative_seconds', 'ee_limit_seconds', 3900),
            ('cubic_cumulative_seconds', 'cubic_limit_seconds', 240)):
        from decimal import Decimal
        if Decimal(prior[field])+additional > Decimal(prior[limit]):
            raise ValueError('declared allocation cannot cover bounded workers: '+field)
    a.output.mkdir(parents=True, exist_ok=False)
    carried = json.loads(authenticate(binding(a.carried_oracle)).read_bytes())
    if carried['status'] != 'passed':
        raise ValueError('an actually passed unchanged callback boundary is required')
    save(a.output/'boundary/oracle.json', dict(
        schema='p254-a56-carried-boundary-v1', status='passed',
        original_oracle=binding(a.carried_oracle),
        scope='Unchanged supplied callback boundary; no old proof replay.'))

    def launch(mode, name, cpu, cap, *extra):
        suffix = hashlib.sha256(str(a.output.resolve()).encode()).hexdigest()[:8]
        argv = [sys.executable, str(HERE.parent/'0056/launch.py'),
                mode, name+'-'+suffix, str(cpu), str(cap),
                str(a.dataset), str(a.baseline_catalog), str(a.output), *map(str, extra)]
        if subprocess.call(argv):
            raise RuntimeError('native component failed; paid checkpoint retained: '+name)

    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        work = [pool.submit(launch, 'ee', 'ee-'+chr(97+i)+'-001', 6+i, 1300, lo, hi)
                for i, (lo, hi) in enumerate(((0, 22), (22, 43), (43, 64)))]
        for task in work:
            task.result()
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        work = [pool.submit(launch, 'mixed', 'mixed-001', 6, 450),
                pool.submit(launch, 'cubic', 'cubic-001', 8, 240)]
        for task in work:
            task.result()
    launch('assemble', 'assemble-001', 6, 120)
    finalize(a.output, a.ledger)


if __name__ == '__main__':
    main()
