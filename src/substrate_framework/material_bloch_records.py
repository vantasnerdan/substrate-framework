"""Complete source/target/frame descriptors for authenticated material coverage.

Original Gamma functional basis identity is a separately named historical
join, never substituted for the new transported basis/gauge identity.
"""
from dataclasses import asdict, dataclass
import hashlib
import json
import numpy as np
from material_pilot_io import array_record
from material_run_contract import canonical
from .qe_bloch_frame import BlochBasis, _basis_identity, _digest
from .qe_bloch_grid import readonly
from .qe_material_grid import PhysicalRhoFFT


@dataclass(frozen=True)
class BlochSlotRecord:
    ispin: int
    mesh_index: tuple
    mesh_slot: int
    source_ik: int
    source_binding: dict
    source_frame_identity: str
    source_coefficients_record: dict
    source_miller_record: dict
    target_miller_record: dict
    row_map_record: dict
    basis_identity: str
    gauge_identity: str
    n_pw: int
    n_bands: int
    operation_index: int
    reciprocal_shift: tuple
    source_k_bohr_inv: tuple
    target_k_bohr_inv: tuple
    raw_band_order: str
    row_order: str
    original_representative: bool
    gamma_basis_join_identity: str | None


def describe_slot(ref, slot, source_miller, source_k, reciprocal, admission, coefficient_record):
    """Rebuild checked target metadata without loading or inventing coefficients."""
    target = readonly(slot.operation.rotation @ source_miller-slot.reciprocal_shift[:, None])
    identity = _basis_identity(target, reciprocal, slot)
    basis = BlochBasis(target, slot.k_cart_bohr_inv, reciprocal, readonly(np.arange(slot.expected_rows)),
                       source_miller, slot, ref.source, admission, identity).validate(ref)
    PhysicalRhoFFT.create(ref.fft_smooth, target, reciprocal, abs(np.linalg.det(ref.cell_bohr)))
    source_p = reciprocal @ source_miller + source_k[:, None]
    target_p = reciprocal @ target + slot.k_cart_bohr_inv[:, None]
    if not np.allclose(np.sum(source_p**2,axis=0), np.sum(target_p**2,axis=0), atol=2e-9, rtol=0):
        raise ValueError('target transport does not preserve source cutoff support')
    original = (np.array_equal(slot.operation.rotation, np.eye(3,dtype=int))
                and not np.any(slot.operation.fractional_translation) and not np.any(slot.reciprocal_shift)
                and np.array_equal(source_miller,target) and np.array_equal(source_k,slot.k_cart_bohr_inv))
    join = None
    if original and not np.any(source_k):
        # A24 material_pilot_inputs.coordinate_arrays stored XML B, not header B.
        # Preserve the latter separately; never claim the two digests are equal.
        join = hashlib.sha256(canonical(dict(pw_miller=array_record(source_miller),
                  cell_bohr=array_record(ref.cell_bohr),
                  reciprocal_bohr_inv=array_record(ref.reciprocal_bohr_inv)))).hexdigest()
    return BlochSlotRecord(slot.ispin, slot.mesh_index, slot.mesh_slot, slot.source_ik,
        dict(slot.source_binding), slot.source_binding['sha256'], coefficient_record,
        array_record(source_miller), array_record(target), array_record(basis.row_map), identity,
        identity+':raw-saved-columns', slot.expected_rows, len(slot.eigenvalues_ry),
        slot.operation.index, tuple(int(x) for x in slot.reciprocal_shift), tuple(source_k.tolist()),
        tuple(slot.k_cart_bohr_inv.tolist()), 'original-all-saved-columns', 'inherited-source-order',
        bool(original), join)


def encode_records(records):
    return np.asarray(json.dumps([asdict(record) for record in records], sort_keys=True, allow_nan=False))


def decode_records(value):
    raw = json.loads(str(value))
    if not isinstance(raw,list):
        raise ValueError('complete typed slot descriptor list required')
    records = []
    for item in raw:
        if set(item) != set(BlochSlotRecord.__dataclass_fields__):
            raise ValueError('incomplete or unknown slot descriptor fields')
        for name in ('ispin','mesh_slot','source_ik','n_pw','n_bands','operation_index'):
            if type(item[name]) is not int:
                raise ValueError('integer slot dimension/index required: '+name)
        if type(item['original_representative']) is not bool:
            raise ValueError('typed original-representative flag required')
        for name in ('mesh_index','reciprocal_shift'):
            if len(item[name]) != 3 or any(type(x) is not int for x in item[name]):
                raise ValueError('integer3 mesh/fold descriptor required')
        for name in ('source_k_bohr_inv','target_k_bohr_inv'):
            if len(item[name]) != 3 or not np.isfinite(item[name]).all():
                raise ValueError('finite physical k descriptor required')
        for name in ('mesh_index','reciprocal_shift','source_k_bohr_inv','target_k_bohr_inv'):
            item[name] = tuple(item[name])
        records.append(BlochSlotRecord(**item))
    return tuple(records)


def source_arrays(source, reference, admission, records, coefficient_record):
    header = {key: getattr(source,key) for key in ('ik','ispin','gamma_only','scalef','ngw','igwx','npol','nbnd','abi')}
    header.update(source=dict(source.source), source_coefficients_record=coefficient_record)
    return dict(source_miller=source.miller, source_k_bohr_inv=source.k_cart_bohr_inv,
                reciprocal_bohr_inv=source.reciprocal_bohr_inv, cell_bohr=reference.cell_bohr,
                reference_reciprocal_bohr_inv=reference.reciprocal_bohr_inv,
                energies_ry=reference.select(source.ik,source.ispin)[1],
                raw_occupations=reference.select(source.ik,source.ispin)[2],
                slot_records_json=encode_records(records), source_header_json=np.asarray(json.dumps(header,sort_keys=True)),
                source_admission_json=np.asarray(json.dumps(dict(admission),sort_keys=True)))


def validate_source_arrays(ref, slots, arrays):
    expected = {'source_miller','source_k_bohr_inv','reciprocal_bohr_inv','reference_reciprocal_bohr_inv',
                'cell_bohr','energies_ry','raw_occupations','slot_records_json','source_header_json','source_admission_json'}
    if set(arrays) != expected or not slots:
        raise ValueError('complete original source archive layout required')
    first = slots[0]
    point, energies, occupations = ref.select(first.source_ik, first.ispin)
    header, admission = (json.loads(str(arrays[name])) for name in ('source_header_json','source_admission_json'))
    if (set(header) != {'ik','ispin','gamma_only','scalef','ngw','igwx','npol','nbnd','abi','source','source_coefficients_record'}
            or header['ik'] != first.source_ik or header['ispin'] != first.ispin
            or header['source'] != dict(first.source_binding) or header['gamma_only'] is not False
            or header['npol'] != 1 or header['scalef'] != 1 or header['ngw'] < point.npw
            or header['igwx'] != point.npw or header['nbnd'] != len(energies)
            or header['abi'] != 'little-endian:marker4:int32:logical32:float64:complex128'):
        raise ValueError('saved original header/source identity mismatch')
    coefficient = header['source_coefficients_record']
    if (set(coefficient) != {'dtype','shape','sha256_C_bytes'} or coefficient['dtype'] != '<c16'
            or coefficient['shape'] != [point.npw,len(energies)] or len(coefficient['sha256_C_bytes']) != 64):
        raise ValueError('all-column coefficient identity required')
    m = arrays['source_miller']
    if m.shape != (3,point.npw) or m.dtype != np.dtype('int32'):
        raise ValueError('all original integer32 source rows required')
    for key,value in (('cell_bohr',ref.cell_bohr),('reference_reciprocal_bohr_inv',ref.reciprocal_bohr_inv),
                      ('energies_ry',energies),('raw_occupations',occupations)):
        if not np.array_equal(arrays[key],value):
            raise ValueError('saved source coordinate/band metadata differs: '+key)
    for key,value in (('source_k_bohr_inv',point.k_cart_bohr_inv),('reciprocal_bohr_inv',ref.reciprocal_bohr_inv)):
        if arrays[key].shape != value.shape or not np.isfinite(arrays[key]).all() or not np.allclose(arrays[key],value,atol=2e-13,rtol=0):
            raise ValueError('saved source physical reciprocal identity differs')
    if (admission['source_rows_sha256'] != _digest(m) or admission['cutoff_ry'] != 2*ref.ecutwfc_ha
            or admission['boundary_tolerance_ry'] != 1e-8
            or admission['row_policy'] != 'complete strict interior; original native cutoff boundary retained'
            or min(admission['included_boundary_rows'],admission['omitted_boundary_rows']) < 0):
        raise ValueError('original source admission/cutoff policy differs')
    actual = decode_records(arrays['slot_records_json'])
    rebuilt = tuple(describe_slot(ref,slot,m,arrays['source_k_bohr_inv'],arrays['reciprocal_bohr_inv'],
                                  admission,coefficient) for slot in slots)
    if actual != rebuilt:
        raise ValueError('complete target/frame descriptors differ from admitted source reconstruction')
    return rebuilt
