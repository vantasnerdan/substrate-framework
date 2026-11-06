"""Audit every public artifact, canonical spec fingerprint, and source binding."""
import copy
import hashlib
import json
from pathlib import Path


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def links_valid(records, root, inventory):
    return all(Path(r['path']).is_relative_to(root) and inventory.get(r['path']) == r['sha256']
               and digest(r['path']) == r['sha256'] for r in records.values())


def run(ledger, context):
    plan, supplement = context['plan'], context['supplement']
    inventory = {**plan['input_sources_sha256'], **supplement['supplemental_sources_sha256']}
    ledger.check('raw-byte native inventory matches frozen revision2 and supplement',
                 all(digest(p) == h for p,h in inventory.items()))
    published, scf, spec, relax = (context[k] for k in ('published','scf','spec','relax'))
    for name,result,source in [('relax',published,relax),('scf',scf,spec)]:
        ledger.check(name+' complete envelope bound to exact canonical specification',
                     result['status']=='complete' and result['schema_version']==1 and result['job_id']==source['id']
                     and result['input_sha256']==fingerprint(source))
        records = result['data']['artifacts']
        predicate = lambda r: links_valid(r, Path(source['directory']), inventory)
        altered = copy.deepcopy(records)
        next(iter(altered.values()))['sha256'] = '0'*64
        outside = copy.deepcopy(records)
        next(iter(outside.values()))['path'] = '/tmp/not-this-stage'
        ledger.mutation_sensitive(name+' artifact path/hash publication', predicate, records, [altered,outside])
        declared = result['data']['input_files_sha256']
        ledger.mutation_sensitive(name+' complete declared source provenance',
            lambda d: d == source['input_files_sha256'] and all(inventory.get(p)==h for p,h in d.items()),
            declared, [{**declared,next(iter(declared)):'0'*64}, {}])
    ledger.check('SCF legacy public artifact set has no XML', set(scf['data']['artifacts']) ==
                 {'command.json','final.extxyz','pw.in','pw.stdout','pw.stderr','spec.json'})
    ledger.check('attempt spec equals force spec', context['attempt_spec']==spec)
    ledger.check('force settings/execution inherit registered primary settings',
        spec['settings']==relax['scf_settings']==context['settings']['scf_settings'] and
        spec['execution']==relax['execution']==context['settings']['execution'])
    command = context['command']
    ledger.check('native command and environment publication is exact', command['argv']==scf['data']['command']
        and command['environment']==scf['data']['environment'] and command['cwd']==scf['data']['attempt_directory']
        and command['mpi_launcher_sha256']==inventory[str(Path(command['argv'][0]).resolve())])
    executable = spec['execution']['pw_executable']
    ledger.check('native executable aliases bind raw inventory', executable in command['argv'] and digest(executable)==inventory[str(Path(executable).resolve())])
    geometry = published['data']['geometry']
    ledger.check('force directory is physical geometry fingerprint', Path(spec['directory']).parent.name == fingerprint(geometry))
    optimizer = context['optimizer']
    ledger.check('optimizer checkpoint bound and uncorrupted', optimizer['input_sha256']==published['input_sha256']
                 and fingerprint(optimizer['state'])==optimizer['state_sha256'])
    return {'native_inventory': inventory, 'count':len(inventory), 'xml_consumed':False,
            'artifact_links': {'relax':published['data']['artifacts'],'scf':scf['data']['artifacts']},
            'spec_fingerprints': {'relax':fingerprint(relax),'scf':fingerprint(spec)},
            'optimizer_state_keys':list(optimizer['state'])}
