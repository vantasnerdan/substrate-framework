"""Committed full-transition output states; resume at completed physical times only."""
from __future__ import annotations
import json
import os
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / '0020'))
from pilot_inputs import binding
from pilot_checkpoints import save_json


def implementation():
    records = {}
    for name, module in list(sys.modules.items()):
        if name.startswith(('substrate_framework', 'catdft')) or name in ('pilot_inputs', 'pilot_checkpoints', 'time_pilot_resume'):
            filename = getattr(module, '__file__', None)
            if filename and Path(filename).suffix == '.py':
                path = Path(filename).resolve()
                records[str(path)] = binding(path)
    for path in (Path(__file__).resolve(), HERE/'run_local_kick.py'):
        records[str(path)] = binding(path)
    return list(records.values())


def save_sample(target, point, index, tau, state, observable, mechanical, fingerprint, event_time):
    name = f'{point}-t{index:02d}'
    archive = target/(name+'.npz')
    receipt_path = target/(name+'.json')
    if receipt_path.exists():
        raise ValueError('committed physical-time sample cannot be overwritten')
    temporary = archive.with_suffix('.npz.tmp')
    with temporary.open('wb') as stream:
        np.savez_compressed(stream, state=state, density=observable.density, ionic_gradient=observable.ionic_gradient)
    os.replace(temporary, archive)
    mechanical = dict(mechanical)
    value = mechanical['ionic_frozen_quadratic_ry']
    mechanical['ionic_frozen_quadratic_ry'] = [value.real, value.imag]
    record = {'name': name, 'point': point, 'index': str(index), 'tau_ry': float(tau),
              'fingerprint': fingerprint, 'state_size': str(state.size), 'archive': binding(archive),
              'mechanical': mechanical, 'committed_at': event_time}
    # Immutable sample receipt commits the state; rotating pointer is convenience only.
    save_json(receipt_path, record)
    save_json(target/(point+'-state.json'), record)
    return record


def recover_samples(target, fingerprint, times):
    result = {}
    for point in ('gamma', 'half'):
        for index, tau in enumerate(times):
            name = f'{point}-t{index:02d}'
            receipt = target/(name+'.json')
            if not receipt.exists():
                continue
            record = json.loads(receipt.read_bytes())
            if (record['name'] != name or record['point'] != point or int(record['index']) != index
                    or record['tau_ry'] != float(tau) or record['fingerprint'] != fingerprint
                    or record['archive'] != binding(target/(name+'.npz'))):
                raise ValueError('physical-time checkpoint source/bytes differ: '+name)
            with np.load(record['archive']['path'], allow_pickle=False) as stored:
                state = stored['state']
                if (state.shape != (int(record['state_size']),) or state.size <= 12
                        or not np.iscomplexobj(state) or not np.all(np.isfinite(state))):
                    raise ValueError('invalid complete transition state: '+name)
            result[name] = record
        present = [index for index in range(len(times)) if f'{point}-t{index:02d}' in result]
        if present != list(range(len(present))):
            raise ValueError('checkpoint physical-time sequence has a hole')
    return result


def cell_trajectory(target, records, masses, times):
    displacement, momentum = [], []
    for point in ('gamma', 'half'):
        us, ps = [], []
        for index in range(len(times)):
            with np.load(records[f'{point}-t{index:02d}']['archive']['path'], allow_pickle=False) as stored:
                us.append(stored['state'][:6].copy())
                ps.append(stored['state'][6:12].copy())
        displacement.append(us)
        momentum.append(ps)
    # q indices 0 and 1 are exactly the complete 2x1x1 grid. No atom-position phase.
    uq, pq = np.asarray(displacement), np.asarray(momentum)
    ucell, pcell = np.fft.ifft(uq, axis=0), np.fft.ifft(pq, axis=0)
    archive = target/'two-cell-trajectory.npz'
    with archive.open('wb') as stream:
        np.savez_compressed(stream, tau_ry=times, displacement_q=uq, momentum_q=pq,
                            displacement_cell=ucell, velocity_cell=pcell/masses, masses_ry=masses)
    return binding(archive)
