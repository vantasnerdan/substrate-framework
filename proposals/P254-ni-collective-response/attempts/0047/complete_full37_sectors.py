"""Complete disjoint original full37 sectors, reusing every paid onsite atom."""
import argparse
import json
from pathlib import Path
import sys
from time import perf_counter

HERE=Path(__file__).resolve().parent;A43=HERE.parent/'0043'
sys.path.insert(0,str(A43))
import numpy as np
from cubic_material_source import restore_original,_new_sector
from material_electron_ion_inputs import load_material
from reduced_projection_storage import ImmutableItems
from reduced_projection_cache import ProjectionCache,atomic_json
from substrate_framework.streamed_material_state import authenticate,file_binding


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cpu',type=int,required=True,choices=(6,7))
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args();started=perf_counter()
    owned=list(range(20,74)) if args.cpu==6 else list(range(74,128))
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=args.resume)
    source_path=Path('/var/tmp/p254-a43-full-state-third-002/source.json')
    source_binding=file_binding(source_path)
    if source_binding['sha256']!='064bc5e4786cde656762fc47830c4ad83819c514d1be44078402b362dc01578d' or source_binding['bytes']!=278906:
        raise ValueError('only original banked A43 source admitted')
    source=json.loads(source_path.read_bytes())
    old_cache=json.loads(Path(authenticate(source['sector_cache'])).read_bytes())
    old_reader=ImmutableItems(old_cache['items'])
    c,manifest,binding=load_material()
    control,inputs,rows,states,original,lower,context=restore_original(manifest)
    reader=ImmutableItems(manifest['items'])
    records={record.identity:record for record in original.slots}
    scientific=[file_binding(A43/name) for name in ('cubic_material_source.py','cubic_electronic_donor.py','cubic_position_moments.py')]
    if scientific!=source['source_bindings']:
        raise ValueError('immutable banked scientific implementation differs')
    specification=dict(kind='disjoint-original-full37-electronic-sector-completion',owned_sector_indices=owned,
        cpu=args.cpu,old_complete_selected_atom_source=source_binding,coefficient_manifest=binding,
        scientific_source_bindings=scientific,executor_source=file_binding(__file__))
    cache=ProjectionCache(output/'sector-completion',specification)
    new_count=reused_count=0
    for index in owned:
        unit_started=perf_counter();key='full/%03d'%index
        old_item=old_cache['items']['sector/%03d'%index]
        paid_atoms=old_item['metadata'].get('electronic_atom_indices',list(range(37)))
        if paid_atoms!=[0,31]:
            raise ValueError('disjoint ownership must exclude all paid full37 sectors')
        missing=[atom for atom in range(37) if atom not in paid_atoms]
        lineage=source['full_sector_lineage'][index]
        metadata=dict(original_sector_lineage=lineage,electronic_atom_indices=list(range(37)),
            reused_atom_indices=paid_atoms,reused_original_archive=old_item['archive'],computed_missing_atoms=missing)
        if cache.has(key):
            reused_count+=1
        else:
            paid=old_reader.get('sector/%03d'%index)
            slot=rows[index][1];record=records[lineage['identity']]
            if record.frame!=lineage['frame'] or record.coefficient!=lineage['coefficient']:
                raise ValueError('original full raw column/state source differs')
            new=_new_sector(index,slot,record,lineage,reader,control,inputs,original,missing)
            if set(new)!={'local','nonlocal_term'} or any(value.shape!=(35,3,3,3) or not np.isfinite(value).all() for value in new.values()):
                raise ValueError('complete original full-column missing35 finite third jets required')
            arrays={}
            for name in ('local','nonlocal_term'):
                merged=np.empty((37,3,3,3),complex)
                merged[missing]=new[name];merged[paid_atoms]=paid[name];arrays[name]=merged
            cache.put(key,arrays,metadata=metadata,timings=dict(seconds=perf_counter()-unit_started))
            new_count+=1
            del new,paid,arrays
        if cache.manifest['items'][key]['metadata']!=metadata:
            raise ValueError('cached full37 actual original attribution differs')
        saved=cache.get(key)
        if set(saved)!={'local','nonlocal_term'} or any(value.shape!=(37,3,3,3) or not np.isfinite(value).all() for value in saved.values()):
            raise ValueError('saved full37 sector payload incomplete')
        print(json.dumps(dict(completed_sector=index,identity=lineage['identity'],cpu=args.cpu,
            full37_archive=cache.manifest['items'][key]['archive'],newly_computed=new_count,reused=reused_count,
            sector_elapsed_seconds=perf_counter()-unit_started,elapsed_seconds=perf_counter()-started,
            actual_saved_norms={name:float(np.linalg.norm(value)) for name,value in saved.items()})),flush=True)
        del saved
    if cache.manifest['status']!='complete':
        cache.manifest['status']='complete';cache.publish()
    result=dict(schema='p254-a47-owned-full37-sector-completion-v1',status='complete-at-owned-sectors',
        cpu=args.cpu,owned_sector_indices=owned,completed_sector_count=len(owned),newly_computed=new_count,reused=reused_count,
        all_raw_columns_complexD_and_original_weights=True,electronic_atom_indices=list(range(37)),
        coefficient_manifest=binding,original_selected_atom_source=source_binding,original_lower_order_authority=lower,
        scientific_source_bindings=scientific,executor_source=file_binding(__file__),completed_cache=file_binding(cache.path),
        full37_archives={key:item['archive'] for key,item in cache.manifest['items'].items()},
        original_A43_mutated=False,whole128_source_complete=False,elapsed_seconds=perf_counter()-started)
    path=output/'worker-result.json'
    if path.exists():
        previous=json.loads(path.read_bytes())
        if not args.resume or previous['completed_cache']!=result['completed_cache']:
            raise ValueError('completed owned output differs')
    else:
        atomic_json(path,result)
    print(json.dumps(dict(result=file_binding(path),newly_computed=new_count,reused=reused_count)),flush=True)

if __name__=='__main__':
    main()
