"""Expose snapshot-graph to real scientific-consumer byte substitution failures.

All changed originals are isolated copies inside this execution generation.
Neither live registered inputs nor their actual frozen snapshots are mutated.
"""
import copy
import hashlib
from pathlib import Path
from verified_sources import VerifiedSources,SourceBindingError,close_source_records


def save_bytes(path,raw):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('xb') as stream:
        stream.write(raw)


def isolated_source(directory,role,original,expected,raw):
    shadow=directory/'substituted-originals'/role
    snapshot=directory/'snapshot-objects'/expected
    save_bytes(shadow,raw+b'\n# isolated raw-binding substitution\n')
    if not snapshot.exists():
        save_bytes(snapshot,raw)
    # Exactly the production graph closure: the snapshot, not the shadow, is verified.
    records=close_source_records({(str(shadow),expected)},{(str(shadow),expected):str(snapshot)})
    return {'path':str(shadow),'sha256':expected},records[0],{
        'role':role,'registered_original':original,'registered_sha256':expected,
        'isolated_original':str(shadow),'isolated_original_sha256':hashlib.sha256(shadow.read_bytes()).hexdigest(),
        'snapshot':str(snapshot),'snapshot_sha256':hashlib.sha256(snapshot.read_bytes()).hexdigest()}


def binding_rejection(call):
    try:
        call()
    except SourceBindingError as error:
        return str(error)
    raise AssertionError('Actual scientific source consumer accepted substituted bytes')


def run(ledger,context,units_ref,espresso_ref):
    from check_units import source_functions
    sources=context['sources']; directory=context['directory']/'source-consumption-witness'
    originals=[units_ref,espresso_ref]; fixtures=[]; bindings=[]; proof=[]
    for role,ref in zip(('units.py','espresso.py'),originals):
        raw=sources.read_bytes(ref['path'],ref['sha256'])
        fixture,record,evidence=isolated_source(directory,role,ref['path'],ref['sha256'],raw)
        fixtures.append(fixture);bindings.append(record);proof.append(evidence)
    authentic=VerifiedSources(bindings)
    _,executed=source_functions(authentic,*fixtures,'p254_isolated_authentic_units')
    ledger.check('snapshot-backed graph supplies authentic separately resolved ASE consumers',
                 all(row['isolated_original_sha256']!=row['registered_sha256']==row['snapshot_sha256']
                     for row in proof) and
                 {row['consumed_path'] for row in authentic.consumed}=={row['snapshot'] for row in proof})
    # Exercise the actual frozen source energy function on a complete native line,
    # not a literal stand-in or mocked function result.
    baseline=next(row for row in context['inventory']['raw_registered_sources']
                  if '/material/loaded-ni31h6/forces/' in row['path'] and row['path'].endswith('/scf/result.json'))
    result=sources.json(baseline['path'],baseline['sha256'])
    stdout=result['data']['artifacts']['pw.stdout']
    lines=sources.text(stdout['path'],stdout['sha256']).splitlines()
    energy_index=next(i for i,line in enumerate(lines) if line.startswith('!'))
    actual_energy=executed['_get_energy'](lines,energy_index)
    ledger.check('authentic snapshot source executes against actual native energy publication',
                 abs(actual_energy-result['data']['energy_ev'])<3e-10)
    for index,role in enumerate(('units','espresso')):
        replaced=copy.deepcopy(bindings)
        replaced[index]['consumed_path']=fixtures[index]['path']
        error=binding_rejection(lambda:source_functions(VerifiedSources(replaced),*fixtures,
                                                       'p254_rejected_'+role+'_units'))
        ledger.check(role+' substitution rejected through real source_functions consumer before unbound execution',
                     bool(error) and hashlib.sha256(Path(proof[index]['snapshot']).read_bytes()).hexdigest()==originals[index]['sha256'])
        proof[index]['consumer_rejection']=error
    other=exercise_data_consumers(ledger,context,directory,result)
    return {'status':'authentic snapshot consumers accepted; substituted original bytes rejected',
            'ase_sources':proof,'native_data_consumers':other,'authentic_native_energy_ev':actual_energy,
            'authentic_consumption':authentic.consumed,
            'boundary':'Graph closure and all real consumers use shared path+expected-SHA byte acceptance. Comment-only originals differ while snapshots remain exact. No live source mutation.'}


def exercise_data_consumers(ledger,context,directory,result):
    sources=context['sources']; artifacts=result['data']['artifacts']
    spec=sources.json(artifacts['spec.json']['path'],artifacts['spec.json']['sha256'])
    geometry=sources.record(spec['structure']['path'])
    pseudo=Path(spec['execution']['pseudo_dir'])/spec['execution']['pseudopotentials']['Ni']
    refs=[('result.json',sources.record(next(row['path'] for row in context['inventory']['raw_registered_sources']
           if '/material/loaded-ni31h6/forces/' in row['path'] and row['path'].endswith('/scf/result.json'))),'json'),
          ('spec.json',artifacts['spec.json'],'json'),('command.json',artifacts['command.json'],'json'),
          ('pw.in',artifacts['pw.in'],'text'),('pw.stdout',artifacts['pw.stdout'],'text'),
          ('final.extxyz',artifacts['final.extxyz'],'text'),('geometry.extxyz',geometry,'text'),
          ('ni.UPF',sources.record(pseudo),'read_bytes')]
    proof=[]
    for role,ref,method in refs:
        raw=sources.read_bytes(ref['path'],ref['sha256'])
        fixture,record,evidence=isolated_source(directory,role,ref['path'],ref['sha256'],raw)
        authentic=VerifiedSources([record])
        getattr(authentic,method)(fixture['path'],fixture['sha256'])
        substituted={**record,'consumed_path':fixture['path']}
        error=binding_rejection(lambda:getattr(VerifiedSources([substituted]),method)(fixture['path'],fixture['sha256']))
        ledger.check(role+' actual '+method+' boundary accepts snapshot and rejects unbound original',
                     bool(error) and authentic.consumed[0]['consumed_path']==evidence['snapshot'] and
                     hashlib.sha256(Path(evidence['snapshot']).read_bytes()).hexdigest()==ref['sha256'])
        proof.append({**evidence,'consumer_method':method,'consumer_rejection':error})
    return proof
