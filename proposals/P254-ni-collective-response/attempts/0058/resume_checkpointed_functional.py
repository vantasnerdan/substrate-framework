"""Finish remaining q work after the registered per-q checkpoint boundary."""
from pathlib import Path
import argparse
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal
import hashlib
import json
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent/'0056'
sys.path.insert(0, str(SOURCE))
from dressed_reference import binding, authenticate, save
from finalize_functional_iteration import finalize


def resume(output, allocation, dataset, baseline, native_limit):
    output = Path(output).resolve()
    plan = json.loads(Path(allocation).read_bytes())
    if binding(dataset) != plan['input_dataset']:
        raise ValueError('checkpoint continuation input changed')
    suffix = '-'+hashlib.sha256(str(output).encode()).hexdigest()[:8]
    archive = output/'checkpoint-inputs-001'
    archive.mkdir(exist_ok=False)
    selected, work, preserved = [], [], []
    prior_ee = Decimal(plan['ee_cumulative_seconds'])
    for group, cpu, lo, hi in (('a', 6, 0, 22), ('b', 7, 22, 43), ('c', 8, 43, 64)):
        old_name = 'ee-'+group+'-001'+suffix
        receipt = SOURCE/'invocations'/old_name
        outcome = json.loads((receipt/'exit.json').read_bytes())
        scope_path = output/('ee-'+group)/'scope.json'
        scope = json.loads(scope_path.read_bytes())
        if scope['input_dataset'] != plan['input_dataset'] or scope['owned_outer_q'] != list(range(lo, hi)):
            raise ValueError('checkpoint ownership/input differs')
        checkpoint = authenticate(scope['archive'])
        target = archive/('ee-'+group+'.npz')
        shutil.copy2(checkpoint, target)
        copied = binding(target)
        if (copied['bytes'], copied['sha256']) != (scope['archive']['bytes'], scope['archive']['sha256']):
            raise ValueError('paid checkpoint snapshot changed')
        preserved.append(dict(original=scope['archive'], snapshot=copied,
                              completed_outer_q=scope['completed_outer_q'], old_exit=binding(receipt/'exit.json')))
        if outcome['exitcode'] == 0 and scope['status'] == 'complete':
            selected.append(old_name)
            continue
        error = (receipt/'stderr').read_text()
        if outcome['exitcode'] != 1 or 'RuntimeError: preserved q checkpoint before cap' not in error:
            raise RuntimeError('different failure requires diagnosis, not unchanged retry: '+old_name)
        if scope['status'] != 'partial' or not scope['completed_outer_q']:
            raise ValueError('no committed partial q work to continue')
        prior_ee += Decimal(str(outcome['elapsed_seconds']))
        name = 'ee-'+group+'-002'+suffix
        selected.append(name)
        work.append(('ee', name, cpu, 600, lo, hi))
    selected += [name+'-002'+suffix for name in ('mixed', 'cubic', 'assemble')]
    plan.update(native_limit_seconds=str(native_limit), ee_cumulative_seconds=str(prior_ee),
                current_map_invocation_names=selected, paid_checkpoint_inputs=preserved,
                checkpoint_continuation_source=binding(__file__))
    final_allocation = output/'checkpoint-continuation-allocation.json'
    save(final_allocation, plan)

    def launch(item):
        mode, name, cpu, cap, *extra = item
        argv = [sys.executable, str(SOURCE/'launch.py'), mode, name, str(cpu), str(cap),
                str(dataset), str(baseline), str(output), *map(str, extra)]
        code = subprocess.call(argv)
        if code:
            raise RuntimeError('bounded checkpoint continuation failed; retain all data: '+name)
    with ThreadPoolExecutor(max_workers=3) as pool:
        for task in [pool.submit(launch, item) for item in work]:
            task.result()
    with ThreadPoolExecutor(max_workers=2) as pool:
        tasks = [pool.submit(launch, (mode, mode+'-002'+suffix, cpu, cap))
                 for mode, cpu, cap in (('mixed', 6, 450), ('cubic', 8, 240))]
        for task in tasks:
            task.result()
    launch(('assemble', 'assemble-002'+suffix, 6, 120))
    return finalize(output, final_allocation)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('output', 'allocation', 'dataset', 'baseline'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--native-limit', type=Decimal, required=True)
    args = parser.parse_args()
    resume(args.output, args.allocation, args.dataset, args.baseline, args.native_limit)
