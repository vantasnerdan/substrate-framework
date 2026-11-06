"""New transport receipt protocol; original matched source/aggregate bindings stay intact."""
from pathlib import Path

from matched_reference_io import HERE, binding, receipt_arrays


def transport_registration():
    from substrate_framework.material_control_reference import bound_json
    from substrate_framework.qe_saved_native import read_bound_bytes
    record = binding(HERE/'matched-transport-registration.json')
    data = bound_json(record)
    for source in data['implementation_sources']:
        read_bound_bytes(source, max_bytes=32 << 20)
    return record, data


def slot_catalog(directories, context, direct_registration, *, transport_registration=None):
    from substrate_framework.material_control_reference import bound_json
    from substrate_framework.material_matched_transport import GEOMETRY_POLICY
    catalog = {}
    for directory in directories:
        for path in sorted(Path(directory).glob('slot-*.json')):
            record = binding(path)
            value = bound_json(record)
            c, meta = value['context'], value['metadata']
            if (value['status'] != 'complete' or any(c.get(k) != v for k, v in context.items())
                    or c['registration'] != direct_registration):
                raise ValueError('different original case/preparation/implementation receipt')
            policy = meta['geometry_policy']
            if policy != 'direct-original-target':
                if (policy != GEOMETRY_POLICY or transport_registration is None
                        or c.get('transport_registration') != transport_registration
                        or meta.get('coefficient_origin') != 'transformed-original-representative'):
                    raise ValueError('unauthenticated transport receipt or unexpected approximation')
            key = (meta['slot']['spin'], tuple(meta['slot']['mesh_index']))
            if key in catalog:
                raise ValueError('duplicate target receipt: '+str(key))
            catalog[key] = record, value
    return catalog


def authenticated_arrays(entry, row, ledger):
    import numpy as np
    from substrate_framework.material_matched_reference import slot_identity
    from substrate_framework.material_matched_transport import NORMALIZATION, require_arrays
    i, slot, population = row
    record, value = entry
    c, meta = value['context'], value['metadata']
    if (meta['slot'] != slot_identity(slot) or meta['basis_identity'] != str(ledger['basis_identity'][i])
            or meta['normalization'] != NORMALIZATION or c['original_source'] != dict(slot.source_binding)
            or c['spin'] != slot.ispin or c['source_ik'] != slot.source_ik
            or meta['rows'] != slot.expected_rows or meta['columns'] != len(population)):
        raise ValueError('original source/basis/normalization/all-column contract differs')
    _, arrays = receipt_arrays(record)
    if max(abs(arrays['number']-population.sum()), np.max(abs(arrays['number_pair']-population.sum()))) > 1e-8*(1+population.sum()):
        raise ValueError('original source population trace differs')
    # Exact nat/FFT geometry checks are performed by the caller's grid.
    require_arrays(arrays, len(arrays['local_onsite_ry_bohr2']), len(arrays['density_bohr3']))
    return arrays


def source_state_receipt(directories, context, registration, spin, source_ik, populations, source):
    import numpy as np
    from substrate_framework.material_control_reference import bound_json
    from substrate_framework.material_reference_state import restore_source
    found = []
    for directory in directories:
        path = Path(directory)/'source-state.json'
        if not path.exists():
            continue
        record = binding(path)
        value = bound_json(record)
        c = value['context']
        if (c.get('spin'), c.get('source_ik')) != (spin, source_ik):
            continue
        if (any(c.get(k) != v for k, v in context.items()) or c['registration'] != registration
                or c['original_source'] != dict(source.source)):
            raise ValueError('source-state original preparation differs')
        _, arrays = receipt_arrays(record)
        if not np.array_equal(arrays['populations'], populations):
            raise ValueError('full source orbit populations differ')
        found.append((record, arrays))
    if not found:
        raise ValueError('completed original source-state.json required; no silent state rebuild')
    record, arrays = found[0]
    for _, other in found[1:]:
        if any(not np.array_equal(arrays[k], other[k]) for k in arrays):
            raise ValueError('multiple differing original source states')
    return record, restore_source(arrays, source_identity=source.source['sha256'])
