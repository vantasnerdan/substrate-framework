"""Source-by-source actual admission and selected-probe checkpoint orchestration."""
import json
from pathlib import Path
import signal
import time
import numpy as np
from material_pilot_io import array_record
from material_run_contract import binding,implementation,verify_consumers
from .qe_saved_reference import load_saved_reference
from .qe_saved_native import load_wavefunctions
from .qe_bloch_grid import build_bloch_grid
from .qe_bloch_frame import admit_source
from .material_bloch_admission import (wave_bindings,source_metadata,final_metadata,
                                       recover_authorized,restore_accepted)
from .material_bloch_records import describe_slot,source_arrays,validate_source_arrays,encode_records
from .material_bloch_probe import selected_slot,probe_arrays,validate_probe_arrays
from .material_bloch_output import load_saved_bloch_coverage,_background_receipt


def arm_boundary(run, seconds):
    run.check()
    remaining = run.grant['remaining_seconds']-run.elapsed()
    if run.summary['first_useful_block_seconds'] is None:
        remaining = min(remaining,run.grant['first_block_seconds']-run.elapsed())
    signal.setitimer(signal.ITIMER_REAL,max(.000001,min(seconds,remaining)))


def produce(run, recipe, paths):
    grant, summary = run.grant, run.summary
    verify_consumers(grant,implementation(paths))
    ref = load_saved_reference(recipe['inputs']['xml'])
    if (ref.nat != 37 or ref.atom_names != ('Ni',)*31+('H',)*6 or ref.nbnd_by_spin != (346,346)
            or ref.fft_grid != (72,72,72)):
        raise ValueError('original actual Ni31H6 material required')
    grid = build_bloch_grid(ref,wave_bindings(recipe['inputs']),symmetry_policy=recipe['settings']['symmetry_policy'])
    if len(grid.slots) != 128:
        raise ValueError('complete128 explicit spin/mesh slots required')
    records, receipt_bindings = [], {}
    for spin in (1,2):
        for ik in range(1,11):
            arm_boundary(run,30)
            started = time.monotonic()
            key = f'source-{spin}-{ik}'
            slots = [slot for slot in grid.slots if (slot.ispin,slot.source_ik)==(spin,ik)]
            metadata = source_metadata(recipe,spin,ik)
            cached = recover_authorized(run,key,metadata)
            if cached is None:
                source = load_wavefunctions(slots[0].source_binding)
                admission = admit_source(ref,slots[0],source)
                coefficient_record = array_record(source.coefficients)
                group = tuple(describe_slot(ref,slot,source.miller,source.k_cart_bohr_inv,
                              source.reciprocal_bohr_inv,admission,coefficient_record) for slot in slots)
                arrays = source_arrays(source,ref,admission,group,coefficient_record)
                run.check()
                # Includes exact dtype/shape/value readback before publication.
                run.commit(key,arrays,metadata)
                del source
            else:
                arrays,_ = cached
            records.extend(validate_source_arrays(ref,slots,arrays))
            receipt_bindings[key] = binding(run.target/(key+'.json'))
            del arrays
            summary.setdefault('source_boundary_seconds',{})[key] = time.monotonic()-started
            run.check()
    records.sort(key=lambda record:(record.mesh_slot,record.ispin))
    if len(records) != 128:
        raise ValueError('incomplete derived slot descriptors')
    verify_consumers(grant,implementation(paths))
    arm_boundary(run,90)
    background_receipt = _background_receipt(recipe['inputs'])
    metadata = final_metadata(recipe,background_receipt)
    cached = recover_authorized(run,'coverage-probe',metadata)
    if cached is None:
        background, restored_receipt = restore_accepted(recipe['inputs'])
        if restored_receipt != background_receipt:
            raise ValueError('restored background differs from admitted original archive')
        slot = selected_slot(grid)
        descriptor = next(record for record in records if (record.ispin,record.mesh_slot)==(slot.ispin,slot.mesh_slot))
        arrays = probe_arrays(background,grid,descriptor,run.check)
        summary['selected_probe_comparisons'] = validate_probe_arrays(arrays,grid,records)
        arrays.update(slot_records_json=encode_records(records),
                      source_receipt_bindings_json=np.asarray(json.dumps(receipt_bindings,sort_keys=True)))
        verify_consumers(grant,implementation(paths))
        run.commit('coverage-probe',arrays,metadata)
        del background,arrays
    else:
        arrays,_ = cached
        stored_records = arrays.pop('slot_records_json')
        if str(stored_records) != str(encode_records(records)):
            raise ValueError('recovered coverage descriptors differ')
        if json.loads(str(arrays.pop('source_receipt_bindings_json'))) != receipt_bindings:
            raise ValueError('recovered coverage receipt dependencies differ')
        summary['selected_probe_comparisons'] = validate_probe_arrays(arrays,grid,records)
        del arrays
    summary['accepted_background'] = background_receipt
    arm_boundary(run,45)
    final_sources = implementation(paths)
    verify_consumers(grant,final_sources)
    summary.update(status='complete',exit_code=0,final_implementation_sources=final_sources)
    # Immutable pre-publication snapshot invokes the SAME public saved reader.
    # Root's post-exit manager outcome remains mandatory and separately recorded.
    snapshot = run.target/('saved-consumption-input-'+summary['run_id']+'.json')
    run.publisher.write(snapshot,summary)
    snapshot_binding = binding(snapshot)
    started = time.monotonic()
    saved = load_saved_bloch_coverage(snapshot_binding)
    summary['saved_consumption'] = dict(input=snapshot_binding,coverage_identity=saved.coverage_identity,
        source_count=len(saved.source_receipts),slot_count=len(saved.slot_records),
        probe_comparisons=dict(saved.probe_comparisons),elapsed_seconds=time.monotonic()-started,
        qualification='authenticated saved reader; post-exit manager success remains separate authority')
    if len(saved.source_receipts) != 20 or len(saved.slot_records) != 128:
        raise ValueError('saved consumption lost complete coverage')
    final_sources = implementation(paths)
    verify_consumers(grant,final_sources)
    summary['final_implementation_sources'] = final_sources
    run.check()
    run.arm()
