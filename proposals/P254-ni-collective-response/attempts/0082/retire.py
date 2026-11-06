"""Retire only explicitly authorized new raw arrays after terminal committed admission."""
from pathlib import Path
import json
import os
import sys
import time
HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'0077'))
from component_accumulator import reopen
from history_banks import binding, authenticate, save


def main(spec_path):
    spec = json.loads(Path(spec_path).read_bytes())
    if spec.get('cpu',10) not in range(6,11) or sorted(os.sched_getaffinity(0)) != [spec.get('cpu',10)]:
        raise ValueError('retirement requires owned authorized CPU')
    records = []
    for item in spec['items']:
        producer_exit = json.loads(authenticate(item['producer_exit']).read_bytes())
        admission_exit = json.loads(authenticate(item['admission_exit']).read_bytes())
        if producer_exit['returncode'] != 0 or admission_exit['returncode'] != 0:
            raise ValueError('retirement before actual Native0')
        result = json.loads(authenticate(item['producer_result']).read_bytes())
        if result.get('status') != 'complete-owned-production':
            raise ValueError('producer savedconsumer did not settle')
        authenticate(result['catalog'])
        authenticate(result['consumer'])
        consumer = json.loads(authenticate(item['bank_consumer']).read_bytes())
        if consumer['status'] != 'actual-saved-bank-consumer-pass' or consumer['catalog'] != result['catalog']:
            raise ValueError('exact contribution saved-bank consumer required')
        # Actual reopen is deliberately independent of historical raw source files.
        before = reopen(spec['bank'])
        catalog = json.loads(authenticate(result['catalog']).read_bytes())
        arrays = catalog.get('sigma', catalog.get('arrays'))
        if not isinstance(arrays,dict):
            raise ValueError('explicit component array descriptors required')
        descriptors = list(arrays.values())
        paths = []
        for descriptor in descriptors:
            path = authenticate(descriptor).resolve()
            if not any(str(path).startswith(prefix) for prefix in ('/var/tmp/p254-a80-','/var/tmp/p254-a81-')) or path.suffix != '.npy' or path.name == 'catalog.json':
                raise ValueError('unauthorized retirement path: '+str(path))
            paths.append((path,descriptor))
        record = dict(schema='p254-a82-staging-retirement-v1',status='authorized-preunlink',spec=binding(spec_path),catalog=result['catalog'],producer_exit=item['producer_exit'],producer_result=item['producer_result'],producer_consumer=result['consumer'],admission_exit=item['admission_exit'],bank_consumer=item['bank_consumer'],committed_pointer=before['pointer'],retired_arrays=[d for p,d in paths],retired_bytes=sum(d['bytes'] for p,d in paths),recovery_loss='Individual-q raw arrays require recomputation for later inspection; committed summed work and per-cell digests/coverage survive.',actual_time=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
        receipt = Path(item['receipt'])
        if receipt.exists():
            raise ValueError('retirement receipt exists; no automatic retry')
        committed_state = authenticate(before['pointer']['current'])
        snapshot = receipt.with_name(receipt.stem+'-committed-state.json')
        with snapshot.open('xb') as stream:
            stream.write(committed_state.read_bytes())
            stream.flush()
            os.fsync(stream.fileno())
        record['committed_generation_snapshot'] = binding(snapshot)
        if any(record['committed_generation_snapshot'][key] != before['pointer']['current'][key] for key in ('bytes','sha256')):
            raise ValueError('committed generation snapshot differs')
        save(receipt,record)
        for path,descriptor in paths:
            path.unlink()
        after = reopen(spec['bank'])
        if after['pointer'] != before['pointer']:
            raise ValueError('bank pointer changed during sole-writer retirement')
        record.update(status='retired-and-actual-reopen-pass',reopened_pointer=after['pointer'])
        save(receipt,record)
        records.append(binding(receipt))
        print(json.dumps(dict(event='retired-and-reopened',receipt=records[-1],bytes=record['retired_bytes'])),flush=True)
    save(spec['result'],dict(schema='p254-a82-retirement-result-v1',status='complete',receipts=records))


if __name__ == '__main__':
    main(sys.argv[1])
