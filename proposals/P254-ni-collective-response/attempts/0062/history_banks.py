"""Authenticated time-major NPY kernels; explicit ownership and mmap lifetime."""
from pathlib import Path
import hashlib
import json
import os
import numpy as np

SCHEMA = 'p254-common-history-kernels-v1'
ORDERS = ('greater', 'lesser', 'retarded')
REFERENCE = 'bare-mechanical-with-contacts'
PI_CONVENTION = 'full-common-joint-plus-cubic; represented-only-control-separate'


def binding(path):
    p = Path(path).resolve()
    h = hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda: f.read(1048576), b''):
            h.update(chunk)
    return dict(path=str(p), bytes=p.stat().st_size, sha256=h.hexdigest())


def authenticate(entry):
    actual = binding(entry['path'])
    if any(actual[k] != entry[k] for k in ('path', 'bytes', 'sha256')):
        raise ValueError('changed source: ' + entry['path'])
    return Path(entry['path'])


def save(path, record):
    p = Path(path)
    tmp = p.with_suffix(p.suffix + '.tmp')
    with tmp.open('w') as f:
        json.dump(record, f, indent=2, sort_keys=True, allow_nan=False)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())
    tmp.replace(p)


def close_maps(banks):
    for bank in banks.values():
        for array in bank.values():
            array._mmap.close()
    banks.clear()


def complete_coverage(record):
    nt = len(record['times_inverse_ry'])
    all_pairs = {(a, b) for a in range(nt) for b in range(nt)}
    coverage = record['coverage']
    if set(coverage) != {'ee', 'mixed_ei', 'joint', 'cubic'}:
        return False
    for name, item in coverage.items():
        pairs = all_pairs if item['owned_pairs'] == 'all' else set(map(tuple, item['owned_pairs']))
        done = all_pairs if item['completed_pairs'] == 'all' else set(map(tuple, item['completed_pairs']))
        if pairs != all_pairs or done != all_pairs:
            return False
        if name != 'joint' and (set(item['owned_q']) != set(range(64)) or
                                set(item['completed_q']) != set(range(64))):
            return False
    return True


class FinalBanks:
    """Fresh owned output only; exceptions leave no complete catalog."""
    def __init__(self, directory, times):
        self.directory = Path(directory).resolve()
        self.directory.mkdir(parents=True, exist_ok=False)
        nt = len(times)
        self.banks = {}
        try:
            for kind, tail in (('sigma', (128, 12, 12)), ('pi', (114, 114))):
                self.banks[kind] = {}
                for order in ORDERS:
                    array = np.lib.format.open_memmap(self.directory / (kind + '-' + order + '.npy'),
                        mode='w+', dtype=np.complex128, shape=(nt, nt, *tail))
                    self.banks[kind][order] = array
                    # Bounded writes, including explicit unowned zeros. These are
                    # never admitted as physical values without coverage metadata.
                    for a in range(nt):
                        for b in range(nt):
                            array[a, b] = 0
            self.closed = False
        except BaseException:
            close_maps(self.banks)
            raise

    def publish(self, record):
        for kind, bank in self.banks.items():
            record[kind] = {}
            for order, array in bank.items():
                array.flush()
                record[kind][order] = dict(binding(self.directory / (kind + '-' + order + '.npy')),
                    shape=list(array.shape), dtype=str(array.dtype))
        record['schema'] = SCHEMA
        record['reference_source_kind'] = REFERENCE
        record['pi_convention'] = PI_CONVENTION
        complete = complete_coverage(record)
        record['status'] = 'complete' if complete else 'partial'
        record['full_solve_eligible'] = complete
        save(self.directory / 'catalog.json', record)
        return record

    def close(self):
        if not self.closed:
            close_maps(self.banks)
            self.closed = True

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


class OwnedKernels:
    """Read-only mappings, owned until close/context exit. Partial is opt-in."""
    def __init__(self, catalog, *, require_complete=True):
        self.catalog = json.loads(Path(catalog).read_bytes())
        record = self.catalog
        if (record['schema'] != SCHEMA or record['reference_source_kind'] != REFERENCE or
                record['pi_convention'] != PI_CONVENTION):
            raise ValueError('kernel information contract differs')
        times = np.asarray(record['times_inverse_ry'], dtype=float)
        if times.ndim != 1 or not len(times) or not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
            raise ValueError('exact increasing finite times required')
        complete = complete_coverage(record)
        if record['status'] != ('complete' if complete else 'partial') or record['full_solve_eligible'] != complete:
            raise ValueError('inconsistent coverage/status')
        if require_complete and not complete:
            raise ValueError('partial kernels are ineligible for a full solve')
        self.times = times
        self.banks = {}
        try:
            for kind, tail in (('sigma', (128, 12, 12)), ('pi', (114, 114))):
                if set(record[kind]) != set(ORDERS):
                    raise ValueError('all three contour orders required')
                self.banks[kind] = {}
                for order in ORDERS:
                    item = record[kind][order]
                    array = np.load(authenticate(item), mmap_mode='r', allow_pickle=False)
                    self.banks[kind][order] = array
                    if (list(array.shape) != item['shape'] or str(array.dtype) != item['dtype'] or
                            array.shape != (len(times), len(times), *tail) or array.dtype != np.complex128):
                        raise ValueError('NPY axes/dtype differ')
            self.sigma, self.pi = self.banks['sigma'], self.banks['pi']
        except BaseException:
            close_maps(self.banks)
            raise

    def close(self):
        close_maps(self.banks)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
