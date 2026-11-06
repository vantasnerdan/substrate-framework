"""Qualify exact declared source-star maps against already paid cubic sectors."""
import json
from pathlib import Path
import numpy as np
from cubic_material_source import restore_original,HERE,atomic_json,restore_basis
from reduced_projection_storage import ImmutableItems
from substrate_framework.streamed_material_state import file_binding,authenticate

manifest_path=Path('/var/tmp/p254-activeband-loaded-fullq-union-001/manifest.json')
manifest_binding=file_binding(manifest_path)
manifest=json.loads(manifest_path.read_bytes())
control,inputs,rows,states,original,lower,context=restore_original(manifest)
records={record.identity:record for record in original.slots}
reader=ImmutableItems(manifest['items'])
cache_path=Path('/var/tmp/p254-a43-full-state-third-002/sector-thirds/manifest.json')
cache=json.loads(cache_path.read_bytes())
sectors=ImmutableItems(cache['items'])
representatives={}
for index,(_,slot,_) in enumerate(rows):
    if not slot.orbital_conjugation and np.array_equal(slot.operation.rotation,np.eye(3)) and not np.any(slot.operation.fractional_translation):
        representatives[(slot.ispin,slot.source_ik)]=index
maps=[]; comparisons=[]
for index,(_,slot,_) in enumerate(rows):
    representative=representatives[(slot.ispin,slot.source_ik)]
    representative_slot=rows[representative][1]
    identity=str(slot.ispin)+':'+','.join(map(str,slot.mesh_index))
    representative_identity=str(representative_slot.ispin)+':'+','.join(map(str,representative_slot.mesh_index))
    record=records[identity]; rep_record=records[representative_identity]
    for binding in (record.frame,record.coefficient):
        authenticate(binding)
    if (record.frame['sha256']!=rep_record.frame['sha256'] or record.coefficient['sha256']!=rep_record.coefficient['sha256']) and not slot.orbital_conjugation:
        raise ValueError('source-star full raw frame/D not byte identical')
    if slot.orbital_conjugation:
        frame,d=original.arrays(record); rep_frame,rep_d=original.arrays(rep_record)
        if not np.array_equal(frame,rep_frame.conj()) or not np.array_equal(d,rep_d.conj()):
            raise ValueError('source-star actual raw conjugation differs')
    geometry_names=('miller','k_cart_bohr_inv','reciprocal_bohr_inv','row_map','source_miller')
    basis=restore_basis(reader.get('preparation/%03d'%index,geometry_names),reader.items['preparation/%03d'%index]['metadata'],slot,control.grid)
    s=slot.operation.rotation
    x=np.linalg.solve(inputs.reference.cell_bohr,inputs.reference.positions_bohr)
    delta=s.T@x-slot.operation.fractional_translation[:,None]-x[:,slot.operation.equivalent_atoms-1]
    residual=delta-np.rint(delta)
    geometry_error=float(np.max(np.abs(inputs.reference.cell_bohr@residual)))
    maps.append(dict(index=index,representative=representative,identity=identity,source_ik=slot.source_ik,spin=slot.ispin,
        rotation=s.tolist(),equivalent_atoms=slot.operation.equivalent_atoms.tolist(),orbital_conjugation=slot.orbital_conjugation,
        source_frame=record.frame,source_coefficient=record.coefficient,source_state=states[(slot.ispin,slot.source_ik)],
        basis_identity=basis.identity,geometry_max_abs_bohr=geometry_error))
    key='sector/%03d'%index; rep_key='sector/%03d'%representative
    if index!=representative and key in cache['items'] and rep_key in cache['items']:
        actual=sectors.get(key); parent=sectors.get(rep_key)
        errors={}
        for name in ('local','nonlocal_term'):
            value=parent[name].conj() if slot.orbital_conjugation else parent[name]
            predicted=np.einsum('ai,bj,ck,nijk->nabc',s,s,s,value[slot.operation.equivalent_atoms-1],optimize=True)
            errors[name]=dict(relative=float(np.linalg.norm(predicted-actual[name])/max(1.,np.linalg.norm(actual[name]))),max_abs=float(np.max(np.abs(predicted-actual[name]))))
        comparisons.append(dict(index=index,representative=representative,errors=errors,geometry_max_abs_bohr=geometry_error))
record=dict(schema='p254-a43-paid-source-star-qualification-v1',coefficient_manifest=manifest_binding,partial_cache=file_binding(cache_path),
    maps=maps,paid_comparisons=comparisons,representative_count=len(representatives),used_operation_count=len({rows[x['index']][1].operation.index for x in maps}),
    tolerance=2e-12,maximum_geometry_error_bohr=max(x['geometry_max_abs_bohr'] for x in maps),oracle=file_binding(__file__))
record['status']='passed' if comparisons and all(error['relative']<2e-12 and error['max_abs']<2e-12 for comparison in comparisons for error in comparison['errors'].values()) else 'failed'
atomic_json(HERE/'actual-paid-source-star-qualification-001.json',record)
print(json.dumps(dict(status=record['status'],representative_count=record['representative_count'],used_operation_count=record['used_operation_count'],maximum_geometry_error_bohr=record['maximum_geometry_error_bohr'],paid_comparisons=comparisons)),flush=True)
if record['status']!='passed':
    raise SystemExit(1)
