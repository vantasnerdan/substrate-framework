"""CompletePi-only publication of committed joint + embedded111 cubic histories.

No Sigma file is opened or required. The bank is read under its existing serial
lock; all output arrays are durable before the atomic catalog publication.
"""
import argparse
from contextlib import contextmanager
import fcntl
import json
import os
from pathlib import Path
import sys
import uuid

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent / n) for n in ('0087', '0062')]
from history_banks import ORDERS, REFERENCE, PI_CONVENTION, binding, authenticate
from bank_storage import load, _open, _new_map, _persist, _sync_dir, _save
from functional_bank import source_identity, complete_pairs

SCHEMA = 'p254-complete-gamma-pi-kernels-v1'
AXES = ['time', 'time', 'coordinate', 'coordinate']
SHAPE = (31, 31, 114, 114)
COMPOSITION = 'zero complex128; add joint Pi; add cubic only [:111,:111]'


@contextmanager
def _snapshot(bank_root):
    root = Path(bank_root).resolve()
    # Existing lock opened read-only: publication never initializes/mutates a bank.
    with (root / 'serial.lock').open('rb') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        pointer, state = load(root)
        if state is None:
            raise ValueError('actual committed dressed bank required')
        yield root, pointer, state


def _read(descriptor, bindings):
    """Authenticate consumed catalog bytes, not its historical dependency tree."""
    authenticate(descriptor)
    old = bindings.get(descriptor['path'])
    if old is not None and any(old[k] != descriptor[k] for k in ('path', 'bytes', 'sha256')):
        raise ValueError('conflicting Pi source catalog binding')
    bindings[descriptor['path']] = descriptor
    return json.loads(Path(descriptor['path']).read_bytes())


def _identity(state):
    identity = state['identity']
    actual = source_identity(identity['input_dataset']['path'],
                             identity['preparation_origin']['path'])
    if actual != identity:
        raise ValueError('changed committed target/original reference/preparation')
    return identity


def _coverage(record):
    coverage = record['coverage']
    if not all(complete_pairs(coverage[key]) for key in ('owned_pairs', 'completed_pairs')):
        raise ValueError('complete owned/completed all961 time pairs required')


def _raw(values, a, b):
    for value in values.values():
        if not np.isfinite(value).all():
            raise ValueError('nonfinite persisted raw Pi')
    expected = values['greater'] - values['lesser'] if a >= b else np.zeros_like(values['retarded'])
    defect = float(np.max(np.abs(values['retarded'] - expected)))
    if defect > 2e-12 * max(1., float(np.max(np.abs(expected)))):
        raise ValueError('raw theta(0)=1 causal identity differs; no refill permitted')
    return defect


def _orders(descriptors, tail, opened):
    if set(descriptors) != set(ORDERS):
        raise ValueError('all three Pi contour orders required')
    result = {}
    for order in ORDERS:
        descriptor = descriptors[order]
        if descriptor['shape'] != [31, 31, *tail] or descriptor['dtype'] != 'complex128':
            raise ValueError('canonical descriptor shape/precision required')
        array = _open(descriptor, (31, 31, *tail), 'complex128')
        opened.append(array)
        if not array.flags.c_contiguous:
            raise ValueError('canonical time-major C-order arrays required')
        result[order] = array
    return result


@contextmanager
def _sources(state, identity, bindings):
    """Consume admitted COW Pi arrays; carry upstream admission without replay."""
    opened = []
    try:
        records, selected = {}, {}
        for name, kind in (('joint', 'pi'), ('cubic', 'cubic')):
            family = state['families'].get(name, {})
            if set(family.get('entries', {})) != {'0'}:
                raise ValueError('complete once-global joint and all64 cubic families required')
            catalog_binding = family['entries']['0']['catalog']
            record = _read(catalog_binding, bindings)
            records[name] = record
            for key in ('input_dataset', 'numerical_reference', 'preparation_origin', 'times_inverse_ry'):
                if record[key] != identity[key]:
                    raise ValueError('Pi source identity differs: ' + name + '/' + key)
            if record.get('failure') is not None or record.get('full_time_square') is not True:
                raise ValueError('incomplete/failed admitted Pi source')
            _coverage(record)
            selected[name] = dict(catalog=catalog_binding, arrays=family['arrays'][kind])
        joint, cubic = records['joint'], records['cubic']
        jc = joint['coverage']
        if (joint['schema'] != 'p254-additive-history-component-v1' or
                joint['status'] != 'complete-global-joint-component' or
                joint['component'] != 'joint' or joint['q_additive'] is not False or
                joint.get('owned_q') or joint.get('completed_q') or
                joint['full_solve_eligible'] is not False or joint['full_mesh'] is not False or
                joint['pi_convention'] != 'full_joint_Gamma_once' or
                joint['sigma_convention'] != 'joint_minus_ionic_Gamma' or
                joint['axes']['pi'] != AXES or
                (joint['slots'], joint['bands'], joint['spins'], joint['gamma_coordinates']) != (128, 12, 2, 114) or
                jc['all114_coordinates'] is not True or jc['external_targets'] != list(range(128))):
            raise ValueError('once-global all128/full114 joint Pi contract differs')
        cc = cubic['coverage']
        if (cubic['schema'] != 'p254-dressed-cubic-family-v1' or cubic['status'] != 'complete' or
                cubic['full_solve_eligible'] is not False or
                cubic['axes'] != ['time', 'time', 'ionic_coordinate', 'ionic_coordinate'] or
                cubic['inherited_q'] != list(range(1, 64)) or cubic['replaced_q'] != [0] or
                cc['owned_q'] != list(range(64)) or cc['completed_q'] != list(range(64)) or
                cc['all111_coordinates'] is not True):
            raise ValueError('complete all64/full111 inherited63/replaced0 cubic required')
        maps = dict(joint=_orders(selected['joint']['arrays'], (114, 114), opened),
                    cubic=_orders(selected['cubic']['arrays'], (111, 111), opened))
        for a in range(31):
            for b in range(31):
                for bank in maps.values():
                    _raw({o: bank[o][a, b] for o in ORDERS}, a, b)
        yield maps, selected
    finally:
        for array in opened:
            array._mmap.close()


def _coverage_record():
    return dict(joint=dict(owned_pairs='all', completed_pairs='all', q_additive=False,
                           all114_coordinates=True, external_targets=list(range(128))),
                cubic=dict(owned_pairs='all', completed_pairs='all', owned_q=list(range(64)),
                           completed_q=list(range(64)), all111_coordinates=True))


def publish_pi(bank_root, output_catalog):
    """Publish three genuine raw Pi arrays only; return the durable catalog Path."""
    with _snapshot(bank_root) as (root, pointer, state):
        identity = _identity(state)
        path = Path(output_catalog).resolve()
        # No publication, temporary metadata, or arrays are written into old storage.
        if path.is_relative_to(root) or path.exists():
            raise ValueError('Pi publication requires a fresh catalog outside the bank')
        protected = set(state['sources']) | {v['path'] for k, v in identity.items() if isinstance(v, dict)}
        if str(path) in protected:
            raise ValueError('Pi publication may not overwrite an upstream source')
        bindings = {}
        with _sources(state, identity, bindings) as (maps, selected):
            path.parent.mkdir(parents=True, exist_ok=True)
            destination = path.parent / (path.stem + '-pi-arrays-' + uuid.uuid4().hex)
            destination.mkdir()
            # A77 collects old state JSON during later Sigma admissions. Retain
            # these exact committed bytes, not a fabricated complete-bank state.
            snapshot_path = destination / 'committed-state.json'
            with snapshot_path.open('xb') as snapshot:
                snapshot.write(authenticate(pointer['current']).read_bytes())
                snapshot.flush()
                os.fsync(snapshot.fileno())
            snapshot_binding = binding(snapshot_path)
            arrays = {}
            for order in ORDERS:
                output_path = destination / ('pi-' + order + '.npy')
                _new_map(output_path, SHAPE, 'complex128')
                output = np.load(output_path, mmap_mode='r+', allow_pickle=False)
                try:
                    # Literal A87 functional_bank.publish Pi operation and loop order.
                    for component in ('joint', 'cubic'):
                        source = maps[component][order]
                        for a in range(31):
                            for b in range(31):
                                if component == 'cubic':
                                    output[a, b, :111, :111] += source[a, b]
                                else:
                                    output[a, b] += source[a, b]
                    for a in range(31):
                        if not np.isfinite(output[a]).all():
                            raise ValueError('nonfinite composed Pi')
                    output.flush()
                finally:
                    output._mmap.close()
                _persist(output_path)
                arrays[order] = dict(binding(output_path), shape=list(SHAPE), dtype='complex128')
            _sync_dir(destination)
            _sync_dir(destination.parent)
            opened = []
            try:
                persisted = _orders(arrays, (114, 114), opened)
                for a in range(31):
                    for b in range(31):
                        _raw({o: persisted[o][a, b] for o in ORDERS}, a, b)
            finally:
                for array in opened:
                    array._mmap.close()
            record = dict(schema=SCHEMA, status='complete-gamma-pi', capability='CompletePi',
                gamma_solve_eligible=True, full_solve_eligible=False, **identity,
                reference_source_kind=REFERENCE, pi_convention=PI_CONVENTION,
                axes=AXES, gamma_coordinates=114, global_q=64, slots=128, bands=12, spins=2,
                coverage=_coverage_record(), composition=COMPOSITION, pi=arrays,
                committed_state=pointer['current'], committed_snapshot=snapshot_binding,
                consumed_sources=selected,
                upstream_costs='Retained unchanged in the bound source catalogs, including timings and boundary_loss; not reset or reclassified.',
                scope='Complete raw Pi for one independent Gamma Picard block; not a full collision bank, electron solve, full G2 or calibrated loss.')
            _save(path, record)
            return path


class PiKernels:
    """Authenticated CompletePi capability with owned read-only mmap lifetime."""
    def __init__(self, catalog_path):
        self.pi = {}
        self.catalog_path = Path(catalog_path).resolve()
        self.catalog = json.loads(self.catalog_path.read_bytes())
        r = self.catalog
        if (r['schema'] != SCHEMA or r['status'] != 'complete-gamma-pi' or
                r['capability'] != 'CompletePi' or r['gamma_solve_eligible'] is not True or
                r['full_solve_eligible'] is not False or 'sigma' in r or
                r['times_inverse_ry'] != list(range(31)) or r['axes'] != AXES or
                r['reference_source_kind'] != REFERENCE or r['pi_convention'] != PI_CONVENTION or
                r['composition'] != COMPOSITION or r['coverage'] != _coverage_record() or
                (r['gamma_coordinates'], r['global_q'], r['slots'], r['bands'], r['spins']) != (114, 64, 128, 12, 2)):
            raise ValueError('complete Gamma-only Pi information contract differs')
        snapshot = r['committed_snapshot']
        if any(snapshot[key] != r['committed_state'][key] for key in ('bytes', 'sha256')):
            raise ValueError('owned snapshot is not the exact original committed state')
        state = json.loads(authenticate(snapshot).read_bytes())
        if state['schema'] != 'p254-dressed-functional-bank-v1':
            raise ValueError('actual committed dressed state required')
        identity = _identity(state)
        if any(r[key] != value for key, value in identity.items()):
            raise ValueError('saved Pi target/original identity differs')
        self.times = np.asarray(r['times_inverse_ry'], dtype=float)
        self.bindings = dict(catalog=binding(self.catalog_path), **identity,
                             committed_snapshot=snapshot, sources=r['consumed_sources'])
        opened = []
        try:
            self.pi = _orders(r['pi'], (114, 114), opened)
            with _sources(state, identity, {}) as (sources, selected):
                if selected != r['consumed_sources']:
                    raise ValueError('saved Pi consumed sources differ from committed state')
                totals = {o: 0. for o in ORDERS}
                maxima = {o: 0. for o in ORDERS}
                defect = 0.
                for a in range(31):
                    for b in range(31):
                        values = {o: self.pi[o][a, b] for o in ORDERS}
                        defect = max(defect, _raw(values, a, b))
                        for order, value in values.items():
                            expected = np.zeros((114, 114), dtype=np.complex128)
                            expected += sources['joint'][order][a, b]
                            expected[:111, :111] += sources['cubic'][order][a, b]
                            if not np.array_equal(value, expected):
                                raise ValueError('saved Pi differs from literal committed joint+embed111cubic')
                            totals[order] += float(np.vdot(value, value).real)
                            maxima[order] = max(maxima[order], float(np.max(np.abs(value))))
                self.report = dict(schema='p254-a89-saved-complete-gamma-pi-consumer-v1',
                    catalog=self.bindings['catalog'], status=r['status'], capability='CompletePi',
                    gamma_solve_eligible=True, full_solve_eligible=False,
                    completed_pairs=961, scalar_elements_consumed_per_order=961 * 114 * 114,
                    all111_plus3_coordinates_consumed=True, all_descriptors_authenticated=True,
                    exact_committed_composition=True, times_inverse_ry=self.times.tolist(),
                    frobenius_norm={o: float(np.sqrt(v)) for o, v in totals.items()},
                    max_abs=maxima, causal_max_abs_defect=defect,
                    sample_30_0={o: [float(v[30, 0, 111, 113].real), float(v[30, 0, 111, 113].imag)]
                                 for o, v in self.pi.items()}, bindings=self.bindings, scope=r['scope'])
        except BaseException:
            for array in opened:
                array._mmap.close()
            self.pi.clear()
            raise

    def close(self):
        for array in self.pi.values():
            array._mmap.close()
        self.pi.clear()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def consume(catalog_path):
    """Read every persisted Pi pair and return the actual authenticated report."""
    with PiKernels(catalog_path) as saved:
        return saved.report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    publish = commands.add_parser('publish')
    publish.add_argument('--bank', type=Path, required=True)
    publish.add_argument('--output', type=Path, required=True)
    reader = commands.add_parser('consume')
    reader.add_argument('--catalog', type=Path, required=True)
    args = parser.parse_args()
    catalog = publish_pi(args.bank, args.output) if args.command == 'publish' else args.catalog
    print(json.dumps(consume(catalog), allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
