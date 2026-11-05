"""Finish paid component data using the declared cumulative allocation.

Historical A56 finalization has a first-map 5000s ceiling. Its immutable
receipts remain authoritative; this owner applies the actual amended ceiling
without repeating kernels or adding a prior total to already-global counters.
"""
from pathlib import Path
import argparse
from collections import Counter
from decimal import Decimal
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent/'0056'
sys.path.insert(0, str(SOURCE))
from dressed_reference import binding, authenticate, save


def finalize(output, allocation):
    output = Path(output).resolve()
    plan = json.loads(Path(allocation).read_bytes())
    previous = json.loads(authenticate(plan['previous_catalog']).read_bytes())
    cat = json.loads((output/'catalog.json').read_bytes())
    authenticate(cat['archive'])
    authenticate(cat['actual_saved_numeric_consumer'])
    if cat['input_dataset'] != plan['input_dataset']:
        raise ValueError('allocation and paid map input differ')
    records = []
    for path in sorted((SOURCE/'invocations').glob('*/exit.json')):
        item = json.loads(path.read_bytes())
        item['manager_exit'] = binding(path)
        item['launch'] = binding(path.parent/'launch.json')
        records.append(item)
    suffix = '-'+hashlib.sha256(str(output).encode()).hexdigest()[:8]
    current = [r for r in records if r['name'].endswith(suffix)]
    if (Counter(r['mode'] for r in current) != Counter(ee=3, mixed=1, cubic=1, assemble=1)
            or any(r['exitcode'] != 0 for r in current)):
        raise ValueError('all six actual component exits required; no replay/fallback')
    total = lambda rows: sum((Decimal(str(r['elapsed_seconds'])) for r in rows), Decimal(0))
    native = total(records)
    ee = Decimal(str(previous['costs']['ee_prior_seconds']))+total(r for r in records if r['mode']=='ee')
    cubic = Decimal(str(previous['costs']['cubic_prior_seconds']))+total(r for r in records if r['mode']=='cubic')
    for value, key in [(native, 'native_limit_seconds'), (ee, 'ee_limit_seconds'), (cubic, 'cubic_limit_seconds')]:
        if value > Decimal(plan[key]):
            raise RuntimeError('actual declared cumulative cap exceeded: '+key)
    ledger = dict(native_cumulative_seconds=str(native), ee_cumulative_seconds=str(ee),
        cubic_cumulative_seconds=str(cubic), current_map_manager_seconds=str(total(current)),
        limits={key: plan[key] for key in ('native_limit_seconds', 'ee_limit_seconds', 'cubic_limit_seconds')},
        all_invocations=records, current_map_invocations=current,
        prior_failures=[r for r in records if r['exitcode'] != 0],
        counter_epoch='All A56 invocation manager receipts; never add the previous global total again.')
    original = binding(output/'catalog.json')
    cat['original_assembly_catalog'] = original
    cat['original_assembly_cost_estimate'] = cat['costs']
    cat['costs'] = ledger
    cat['campaign_ledger'] = ledger
    cat['carried_first_map_preregistration'] = cat['preregistered_source_inventory']
    cat['preregistered_source_inventory'] = binding(allocation)
    cat['iteration_finalizer'] = binding(__file__)
    cat['findings'].append(dict(status='resolved-execution-boundary',
        claim='Physical EE/mixed/cubic/assembly all exited0. The controller failed only its obsolete first-map5000s finalizer ceiling; Root amended9000s before execution. Finalized from paid data without kernel/proof replay.',
        actual_current_component_exits=[r['exitcode'] for r in current]))
    target = output/'iteration-catalog.json'
    save(target, cat)
    handoff = dict(schema='p254-complete-further-dressed-map-v1', status='complete',
        catalog=binding(target), archive=cat['archive'], input_dataset=cat['input_dataset'],
        actual_saved_numeric_consumer=cat['actual_saved_numeric_consumer'],
        campaign_ledger=ledger, findings=cat['findings'],
        original_assembly_catalog=original, allocation=binding(allocation),
        source=binding(__file__),
        scope='One complete further common functional map, not converged evolution or temporal accuracy.')
    save(output/'iteration-handoff.json', handoff)
    print(json.dumps(dict(stage='paid-further-map-finalized', catalog=binding(target),
        cumulative_native=str(native), cumulative_ee=str(ee), cumulative_cubic=str(cubic),
        current_map_manager_seconds=str(total(current)), new_kernel_evaluations=0)), flush=True)
    return handoff


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--allocation', type=Path, required=True)
    a = p.parse_args()
    finalize(a.output, a.allocation)
