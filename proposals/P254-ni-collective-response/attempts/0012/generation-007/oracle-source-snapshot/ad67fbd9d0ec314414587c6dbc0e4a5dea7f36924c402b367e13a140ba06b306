"""Resolve registered raw sources and closed package/snapshot graphs by bytes."""
import hashlib
import json
from pathlib import Path
from verified_sources import read_bound_bytes,close_source_records,SourceBindingError

PLAN=Path('/home/administrator/projects/substrate-ni-collective/proposals/P254-ni-collective-response/W1-loaded-force-sensitivity-plan.json')
PLAN_SHA='81d7f2b69ada5cda84900b64a1e833b0be6a87c2e97d1970f7521b90843f67be'
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def load(path):
    return json.loads(Path(path).read_bytes())


def bound_json(path,expected,snapshots=None):
    actual=snapshots.get((str(path),expected),path) if snapshots is not None else path
    return json.loads(read_bound_bytes(actual,expected))


def fingerprint(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def source_binding_valid(binding):
    path,expected=binding
    try:
        read_bound_bytes(path,expected)
    except SourceBindingError:
        return False
    return True


def records(value):
    if isinstance(value,dict):
        if isinstance(value.get('path'),str) and 'sha256' in value:
            yield value['path'],value['sha256']
        for key,item in value.items():
            if key.startswith('/') and isinstance(item,dict) and 'sha256' in item:
                yield key,item['sha256']
            yield from records(item)
    elif isinstance(value,list):
        for item in value:
            yield from records(item)


def closed_graph():
    plan=bound_json(PLAN,PLAN_SHA); access=Path(plan['source_access']['path'])
    inventory=bound_json(access,plan['source_access']['sha256'])
    report_ref=plan['native_wave']['report']
    report=bound_json(report_ref['path'],report_ref['sha256'])
    references={(str(PLAN),PLAN_SHA),(str(access),plan['source_access']['sha256'])}
    references.update(records(inventory)); references.update(records(report))
    snapshots={}; packages=[]
    for package in report['generation_packages']:
        package_data=bound_json(package['package']['path'],package['package']['sha256']); packages.append(package_data)
        references.update(records(package_data))
        snapshot=bound_json(package['source_snapshot']['path'],package['source_snapshot']['sha256'])
        for original,item in snapshot.items():
            saved=Path(item['snapshot'])
            references.add((str(saved),item['sha256']))
            snapshots[(original,item['sha256'])]=str(saved)
    # Every declared source in every completed native spec, not just 71 roots.
    specs=[]
    for row in inventory['raw_registered_sources']:
        if row['path'].endswith('/scf/result.json'):
            result=bound_json(row['path'],row['sha256'],snapshots)
            spec_ref=result['data']['artifacts']['spec.json']
            spec=bound_json(spec_ref['path'],spec_ref['sha256'],snapshots)
            specs.append(spec)
            for original,sha in spec['input_files_sha256'].items():
                references.add((original,sha))
    verified=close_source_records(references,snapshots)
    for row in inventory['current_method_sources'][:2]:
        copied=ROOT/'src/substrate_framework'/Path(row['path']).name
        if not source_binding_valid((copied,row['sha256'])):
            raise ValueError('Imported method differs from captured current source: '+str(copied))
    return {'plan':plan,'inventory':inventory,'native_report':report,'verified':verified,
            'specs':specs,'packages':packages}


def own_sources():
    paths=list(HERE.glob('*.py'))+[ROOT/'src/substrate_framework/material_force_sensitivity.py',
        ROOT/'src/substrate_framework/material_primary_boundary.py',ROOT/'src/substrate_framework/verification.py']
    return {str(path):digest(path) for path in sorted(paths)}


def run(ledger,context):
    graph=closed_graph(); manifest=context['manifest']
    ledger.check('registered full raw source/model/package graph matches prepared exact bytes',
                 graph['verified']==manifest['upstream_verified'] and own_sources()==manifest['oracle_sources'])
    method=graph['inventory']['current_method_sources'][0]
    imported=ROOT/'src/substrate_framework'/Path(method['path']).name
    other=ROOT/'src/substrate_framework/verification.py'
    ledger.mutation_sensitive('load-bearing current-method source substitution and digest corruption',
        source_binding_valid,(imported,method['sha256']),
        [(other,method['sha256']),(imported,'0'*64)])
    context.update(graph)
    return {'verified_sources':graph['verified'],'registered_raw_count':len(graph['inventory']['raw_registered_sources']),
            'registered_current_method_count':len(graph['inventory']['current_method_sources']),
            'nanosecond_metadata_used_as_proof':False,'snapshot_resolution':'Each declared hash resolved to frozen raw object when available; direct bytes otherwise.'}
