"""Assemble original evidence without editing production/cache artifacts."""
from pathlib import Path
import json
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent/'0063')]
from history_storage import binding, atomic_json


def main():
    producer = json.loads((HERE/'producer-result.json').read_bytes())
    consumer = json.loads((HERE/'consumer-result.json').read_bytes())
    published = Path(producer['result']['path']).with_name('published-mean.json')
    spec = json.loads(published.read_bytes())
    receipts = {}
    seconds = 0.
    for mode in ('producer','consumer'):
        folder = HERE/('receipts-'+mode)
        receipt = {}
        for name in ('launch.json','runtime.json','native-exit.json','exit.json','manager-terminal.json','stdout','stderr','time.txt'):
            path = folder/name
            receipt[name] = binding(path)
            if path.suffix == '.json':
                receipt[name+'_decoded'] = json.loads(path.read_bytes())
        if (receipt['exit.json_decoded']['returncode'] != 0 or
                receipt['native-exit.json_decoded']['returncode'] != 0):
            raise ValueError('terminal failure must not be reported complete')
        seconds += receipt['native-exit.json_decoded']['elapsed_seconds']
        receipts[mode] = receipt
    if seconds > 900:
        raise ValueError('aggregate actual native allocation exceeds900s')
    inventory = [binding(p) for p in sorted(HERE.rglob('*')) if p.is_file() and
        '__pycache__' not in p.parts and p.name not in ('complete-handoff.json','evidence-inventory.json')]
    atomic_json(HERE/'evidence-inventory.json',dict(files=inventory))
    handoff = dict(schema='p254-fullhistory-common-mean-handoff-v1',status='complete',
        preregistration=binding(HERE/'preregistration.json'),result=binding(published),archive=producer['archive'],
        original_executed_result=producer['result'],source_restoration=spec['source_restoration'],
        historical_source_snapshots=binding(HERE/'historical-source-snapshots.json'),
        input_identity=spec['input_identity'],coefficient_manifest=spec['coefficient_manifest'],
        physical_history_bindings=spec['physical_history_bindings'],source_bindings=spec['source_bindings'],
        producer_original_reader=spec['producer_reader'],consumer=consumer['consumer'],
        findings=producer['findings'],consumer_findings=consumer['findings'],receipts=receipts,
        actual_total_native_seconds=seconds,
        allocation=dict(cpu=6,memory_bytes=3221225472,threads=1,no_swap=True,no_restart=True,
            producer_cap_seconds=780,consumer_cap_seconds=120,aggregate_cap_seconds=900),
        timing_gap='Producer and consumer are separate bounded units. Between them only source-reader cutover/snapshot compatibility implementation occurs; no ODE or physics recomputation in gap. Wall-clock launch/exit UTC timestamps retained.',
        force_cache_source_compatibility='0074 cache/result/archive immutable. Explicit identical-byte historical implementation snapshots distinguish as-executed source from current reader; numeric/material/physical bindings authenticate actual original paths.',
        owned_view_closure=binding(HERE/'terminal-view-closure.json'),
        api=binding(HERE/'production-api.json'),inventory=binding(HERE/'evidence-inventory.json'),
        evidence_scope='Actual mean/U/S[G0] on all31 genuine nodes0..30/all128/12/fullphase228 with original initial preparation and all64 covariance force. SavedMean all128/fullGamma actual saved-path exercise and explicit original cache reopen performed. No collisionkernel/fullstate/fixedpoint/fullmesh/physical calibration or temporal accuracy claim.',
        limitations=spec['limits'])
    atomic_json(HERE/'complete-handoff.json',handoff)
    print(json.dumps(dict(handoff=binding(HERE/'complete-handoff.json'),actual_total_native_seconds=seconds)))


if __name__ == '__main__':
    main()
