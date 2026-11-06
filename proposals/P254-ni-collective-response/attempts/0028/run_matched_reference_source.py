#!/usr/bin/env python3
"""One original source, explicit target selection, direct all-column coefficients.

--target all evaluates every original target of this source; repeated --target
x,y,z evaluates only named slots for a bounded pilot. Rerunning the identical
command restores completed source/slot receipts. No transported moments.
"""
import argparse
from pathlib import Path
import time
from matched_reference_io import (binding, save_json, save_arrays, runtime, registration,
                                  load_case, receipt_arrays, loaded_field_cache, cache_for_slot)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True)
    parser.add_argument('--spin', type=int, choices=(1, 2), required=True)
    parser.add_argument('--source-ik', type=int, required=True)
    parser.add_argument('--target', action='append', required=True)
    parser.add_argument('--array-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-array-bytes', type=int, default=1 << 30)
    parser.add_argument('--row-block', type=int, default=256)
    args = parser.parse_args()
    runtime(args.array_root)
    import numpy as np
    from substrate_framework.qe_saved_native import load_wavefunctions
    from substrate_framework.qe_bloch_frame import admit_source
    from substrate_framework.material_reference_state import restore_source
    from substrate_framework.material_matched_reference import source_state, direct_slot
    reg_binding, reg = registration()
    control, background, ledger, rows, context = load_case(args.case, reg)
    orbit = [r for r in rows if (r[1].ispin, r[1].source_ik) == (args.spin, args.source_ik)]
    if not orbit:
        raise ValueError('original source absent')
    if any(not np.array_equal(orbit[0][2],row[2]) for row in orbit):
        raise ValueError('original source orbit populations differ')
    if args.target == ['all']:
        selected = orbit
    else:
        targets = [tuple(map(int, t.split(','))) for t in args.target]
        if len(set(targets)) != len(targets) or any(len(t) != 3 for t in targets):
            raise ValueError('nonduplicated three-component target indices required')
        selected = [r for r in orbit if r[1].mesh_index in targets]
        if len(selected) != len(targets):
            raise ValueError('target does not belong to original source/spin')
    args.output.mkdir(parents=True, exist_ok=True)
    array_dir = args.array_root.resolve()/args.case/f'source-{args.spin}-{args.source_ik}'
    array_dir.mkdir(parents=True, exist_ok=True)
    context.update(registration=reg_binding, original_source=dict(orbit[0][1].source_binding),
                   spin=args.spin, source_ik=args.source_ik, row_block=args.row_block)
    source = load_wavefunctions(orbit[0][1].source_binding)
    admission = admit_source(control.reference, orbit[0][1], source)
    state_path = args.output/'source-state.json'
    cache = None
    if state_path.exists():
        old, arrays = receipt_arrays(binding(state_path))
        if old['context'] != context or not np.array_equal(arrays['populations'], orbit[0][2]):
            raise ValueError('source-state resume identity differs')
        state = restore_source(arrays, source_identity=source.source['sha256'])
    else:
        if control.case == 'loaded-ni31h6':
            state, cache = loaded_field_cache(control,None,source,reg,expected_populations=orbit[0][2])
        else:
            state = source_state(source, admission, orbit, authority='original control ledger '+context['population_ledger']['sha256'])
        arrays = dict(gram=state.gram, density_coefficients=state.density_coefficients,
                      populations=state.populations,
                      density_eigenvalues=np.asarray(state.admissibility['density_eigenvalues']))
        archive = save_arrays(array_dir/'source-state.npz', arrays, args.max_array_bytes)
        save_json(state_path, dict(status='complete', context=context, archive=archive))
    if control.case == 'loaded-ni31h6' and cache is None:
        state, cache = loaded_field_cache(control,state,source,reg,expected_populations=orbit[0][2])
    receipts = []
    for index, slot, population in sorted(selected, key=lambda r: r[1].mesh_index):
        name = 'slot-'+ '-'.join(map(str, slot.mesh_index))
        path = args.output/(name+'.json')
        if path.exists():
            old, _ = receipt_arrays(binding(path))
            if old['context'] != context or old['metadata']['basis_identity'] != str(ledger['basis_identity'][index]):
                raise ValueError('target resume identity differs')
        else:
            began = time.monotonic()
            arrays, meta = direct_slot(control, background, slot, source, state, admission,
                                      row_block=args.row_block, cached_field=cache_for_slot(cache,slot))
            if meta['basis_identity'] != str(ledger['basis_identity'][index]):
                raise ValueError('original target basis ledger differs')
            archive = save_arrays(array_dir/(name+'.npz'), arrays, args.max_array_bytes)
            save_json(path, dict(status='complete', context=context, metadata=meta,
                                 archive=archive, elapsed_seconds=time.monotonic()-began))
        receipts.append(binding(path))
    # Each invocation has its own manifest; slot receipts remain immutable.
    name = 'selection-'+__import__('hashlib').sha256(
        __import__('json').dumps(receipts, sort_keys=True).encode()).hexdigest()+'.json'
    result = dict(status='complete-selected-targets', context=context, slots=receipts,
                  complete_case=False, response_comparison=False)
    if not (args.output/name).exists():
        save_json(args.output/name, result)
    print(__import__('json').dumps(dict(selection=binding(args.output/name)), sort_keys=True))


if __name__ == '__main__':
    main()
