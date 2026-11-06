"""Prospective same-coordinate recovery; original generation is never re-pinned."""
from __future__ import annotations
from dataclasses import fields
import hashlib
import json
import numpy as np

from time_pilot_checkpoints import save_json


def saved_reference(target):
    """Summary is advisory; a killed publish may leave only its previous copy."""
    for name in ('result.json', 'result.json.prev'):
        try:
            return json.loads((target/name).read_bytes()).get('reference')
        except FileNotFoundError:
            continue
    return {'availability': 'unavailable', 'reason': 'both advisory summaries absent; committed states remain authoritative'}


def _array_identity(value):
    value = np.asarray(value)
    return {'dtype': value.dtype.str, 'shape': list(value.shape),
            'sha256_C_bytes': hashlib.sha256(value.tobytes(order='C')).hexdigest()}


def coordinate_identity(electronic):
    """Bind actual bases/maps/operators, not merely input/source names or dimensions.

    An eigenvector phase or degenerate-subspace rotation changes state coordinates
    even when H0 is unchanged. Exact match is sufficient for this recovery path;
    mismatch requires a separately justified coordinate transform, never guessing.
    """
    model = electronic.model
    channels = []
    for channel in model.channels:
        record = {}
        for field in fields(channel):
            value = getattr(channel, field.name)
            record[field.name] = _array_identity(value) if isinstance(value, np.ndarray) else value
        channels.append(record)
    fft = {field.name: (_array_identity(value) if isinstance(value, np.ndarray) else value)
           for field in fields(model.fft) for value in [getattr(model.fft, field.name)]}
    return {'channels': channels, 'fft': fft, 'q': _array_identity(model.q_cart_2pi_alat),
            'state_layout': [vars(block) for block in electronic.state_layout]}


def admit_coordinates(target, point, fingerprint, electronic, *, committed):
    path = target/(point+'-coordinates.json')
    actual = {'fingerprint': fingerprint, 'coordinates': coordinate_identity(electronic)}
    # Normalize tuples exactly as their saved JSON representation.
    actual = json.loads(json.dumps(actual, sort_keys=True, allow_nan=False))
    if path.exists():
        if json.loads(path.read_bytes()) != actual:
            raise ValueError('saved spectral coordinates differ; no amplitude reinterpretation')
    elif committed:
        raise ValueError('original spectral coordinates were not recorded; numerical resume is unproven')
    else:
        save_json(path, actual)


def reserved_actions(path, *, required):
    """Recover atomic reservations; conservatively charge one incomplete publish.

    Root policy: a reservation is charged before the derivative. A malformed
    temporary counter may be the next reservation, so retain its raw bytes in
    the returned audit and charge it once. Corrupt committed counters refuse.
    This is accounting recovery, not an automatic computation retry.
    """
    def decode(raw):
        value = json.loads(raw)['action']
        if not isinstance(value, str) or not value.isdecimal():
            raise ValueError('invalid committed action reservation')
        return int(value)
    count, found = 0, False
    for candidate in (path, path.with_suffix(path.suffix+'.prev')):
        if candidate.exists():
            count, found = decode(candidate.read_bytes()), True
            break
    temporary = path.with_suffix(path.suffix+'.tmp')
    audit = None
    if temporary.exists():
        raw = temporary.read_bytes()
        try:
            count = max(count, decode(raw))
        except (ValueError, KeyError, TypeError):
            count += 1
            audit = {'policy': 'charge one ambiguous unpublished reservation',
                     'preserved_raw_hex': raw.hex(), 'sha256': hashlib.sha256(raw).hexdigest(),
                     'charged_actions': str(count)}
        found = True
    if required and not found:
        raise ValueError('action accounting absent for nonzero-time checkpoint; cap cannot be reset')
    return count, audit


def reserve_action(path, count):
    save_json(path, {'action': str(count)})
