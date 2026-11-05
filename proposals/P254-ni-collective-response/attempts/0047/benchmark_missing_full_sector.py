"""One actual missing full37 sector, reusing paid atoms0/31; no A43 mutation."""
import cProfile
import hashlib
import json
from pathlib import Path
import pstats
import sys
from time import perf_counter

HERE=Path(__file__).resolve().parent
A43=HERE.parent/'0043'
sys.path.insert(0,str(A43))
import numpy as np
from cubic_material_source import restore_original,_new_sector
from material_electron_ion_inputs import load_material
from reduced_projection_storage import ImmutableItems
from reduced_projection_cache import ProjectionCache,atomic_json
from substrate_framework.streamed_material_state import authenticate,file_binding


def main():
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(); started=perf_counter()
    output=args.output.resolve(); output.mkdir(parents=True,exist_ok=False)
    source_path=Path('/var/tmp/p254-a43-full-state-third-002/source.json')
    source_binding=file_binding(source_path)
    if source_binding['sha256']!='064bc5e4786cde656762fc47830c4ad83819c514d1be44078402b362dc01578d' or source_binding['bytes']!=278906:
        raise ValueError('only original banked A43 source admitted')
    source=json.loads(source_path.read_bytes())
    cache_record=json.loads(Path(authenticate(source['sector_cache'])).read_bytes())
    old_reader=ImmutableItems(cache_record['items'])
    sector=19; key='sector/%03d'%sector
    old_item=cache_record['items'][key]
    old_atoms=old_item['metadata']['electronic_atom_indices']
    if old_atoms!=[0,31]:
        raise ValueError('paid selected-atom coverage differs')
    paid=old_reader.get(key)
    missing=[atom for atom in range(37) if atom not in old_atoms]
    c,manifest,binding=load_material()
    control,inputs,rows,states,original,lower,context=restore_original(manifest)
    reader=ImmutableItems(manifest['items'])
    slot=rows[sector][1]
    record=next(record for record in original.slots if record.identity==source['full_sector_lineage'][sector]['identity'])
    lineage=source['full_sector_lineage'][sector]
    if lineage['identity']!='1:1,0,3' or record.frame!=lineage['frame'] or record.coefficient!=lineage['coefficient']:
        raise ValueError('original full raw state/slot binding differs')
    scientific=[file_binding(A43/name) for name in ('cubic_material_source.py','cubic_electronic_donor.py','cubic_position_moments.py')]
    if scientific!=source['source_bindings']:
        raise ValueError('banked A43 scientific source changed')
    specification=dict(kind='single-actual-missing-full37-sector-completion',sector=sector,
        identity=lineage['identity'],old_complete_source=source_binding,coefficient_manifest=binding,
        original_sector_lineage=lineage,paid_selected_atom_archive=old_item['archive'],
        paid_selected_atoms=old_atoms,new_missing_atoms=missing,scientific_source_bindings=scientific,
        benchmark_source=file_binding(__file__),whole128_source_complete=False)
    cache=ProjectionCache(output/'sector-completion',specification)
    profile=cProfile.Profile(); unit_started=perf_counter()
    new=profile.runcall(_new_sector,sector,slot,record,lineage,reader,control,inputs,original,missing)
    compute_seconds=perf_counter()-unit_started
    if set(new)!={'local','nonlocal_term'} or any(value.shape!=(35,3,3,3) or not np.isfinite(value).all() for value in new.values()):
        raise ValueError('actual finite full-column missing35 atom outputs required')
    cache.put('missing/019',new,metadata=dict(original_sector_lineage=lineage,electronic_atom_indices=missing),timings=dict(seconds=compute_seconds))
    saved=cache.get('missing/019')
    assembled={}
    for name in ('local','nonlocal_term'):
        value=np.empty((37,3,3,3),complex)
        value[missing]=saved[name]; value[old_atoms]=paid[name]
        assembled[name]=value
    metadata=dict(electronic_atom_indices=list(range(37)),original_sector_lineage=lineage,
        reused_atoms=old_atoms,reused_archive=old_item['archive'],new_atoms=missing,
        new_missing_archive=cache.manifest['items']['missing/019']['archive'])
    cache.put('full/019',assembled,metadata=metadata)
    consumed=cache.get('full/019')
    if any(value.shape!=(37,3,3,3) or not np.isfinite(value).all() for value in consumed.values()):
        raise ValueError('actual saved full37 sector consumer incomplete')
    cache.manifest['status']='complete';cache.publish()
    profile_path=output/'sector-019.profile';profile.dump_stats(profile_path)
    profile_text=output/'sector-019-profile.txt'
    with profile_text.open('x') as stream:
        pstats.Stats(profile,stream=stream).sort_stats('cumulative').print_stats(40)
    stages=[]
    for (filename,line,function),(primitive,calls,self_seconds,cumulative,callers) in pstats.Stats(profile).stats.items():
        if function in ('_new_sector','full_state_third','full_column_projector_third','compressed_overlaps','position_factors','position_fft','local_third'):
            stages.append(dict(file=filename,line=line,function=function,calls=calls,self_seconds=self_seconds,cumulative_seconds=cumulative))
    result=dict(schema='p254-a47-single-actual-full37-sector-benchmark-v1',status='one-full-sector-complete',
        whole128_source_complete=False,sector=sector,identity=lineage['identity'],full_sector_atom_coverage=list(range(37)),
        all_original_raw_columns_and_complexD=True,reused_paid_atom_indices=old_atoms,
        source=source_binding,scientific_source_bindings=scientific,coefficient_manifest=binding,
        original_sector_lineage=lineage,paid_selected_archive=old_item['archive'],
        new_missing_archive=cache.manifest['items']['missing/019']['archive'],full37_archive=cache.manifest['items']['full/019']['archive'],
        completed_cache=file_binding(cache.path),profile=file_binding(profile_path),profile_text=file_binding(profile_text),
        actual_missing35_compute_seconds=compute_seconds,stage_profile=stages,
        actual_saved_full37_norms={name:float(np.linalg.norm(value)) for name,value in consumed.items()},
        full37_sectors_now_available=20,remaining_full37_sectors=108,
        elapsed_seconds=perf_counter()-started,original_A43_mutated=False)
    atomic_json(output/'sector-result.json',result)
    print(json.dumps(result),flush=True)

if __name__=='__main__':
    main()
