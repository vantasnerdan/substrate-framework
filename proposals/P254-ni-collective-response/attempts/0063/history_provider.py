"""Finalized dressed common history provider consumed by A62's functional."""
from pathlib import Path
import json
import numpy as np
from history_storage import StoredNPZ, AUXILIARY, authenticate, binding, mapped
from production_transport import load_solution
from saved_mean import REFERENCE, PI_CONVENTION

SCHEMA = 'p254-dressed-history-provider-v1'


def close_solution(bank):
    for key in ('ret', 'less'):
        bank[key]._mmap.close()


class SavedHistory:
    """Own the two heterogeneous A61 banks and immutable auxiliary mappings."""
    def __init__(self, provider_path):
        self.path = Path(provider_path).resolve()
        self.spec = json.loads(self.path.read_bytes())
        s = self.spec
        if (s['schema'] != SCHEMA or s['status'] != 'complete' or
                s['reference_source_kind'] != REFERENCE or s['pi_convention'] != PI_CONVENTION or
                s['sourcekind'] != 'single-physical-Hermitian-contour'):
            raise ValueError('complete physical common-history provider required')
        self.times = np.asarray(s['times_inverse_ry'], dtype='<f8')
        self.banks, self.aux = {}, {}
        self.bindings = [binding(self.path)]
        try:
            for name, shape, statistics in (
                    ('electron', (128, len(self.times), len(self.times), 12, 12), -1),
                    ('gamma', (1, len(self.times), len(self.times), 114, 114), 1)):
                entry = s[name]
                complete = authenticate(entry['complete'])
                bank = load_solution(complete.parent)
                self.banks[name] = bank
                if (bank['ret'].shape != shape or bank['less'].shape != shape or
                        not np.array_equal(bank['times'], self.times) or bank['statistics'] != statistics or
                        bank['order'] != s['integration_order'] or bank['beta'] != s['beta']):
                    raise ValueError('canonical final-bank dimensions/times/conventions differ')
                self.bindings.append(entry['complete'])
            for name in AUXILIARY:
                self.aux[name] = mapped(s['auxiliary'][name])
                self.bindings.append(s['auxiliary'][name])
            if (self.aux['source_electron_less0'].shape != (128, 12, 12) or
                    any(self.aux[k].shape != (len(self.times), len(self.times), 114, 114)
                        for k in ('control_ret', 'control_less', 'source_gamma_ret', 'source_gamma_less')) or
                    self.aux['n'].shape != () or int(self.aux['n']) <= 0):
                raise ValueError('original auxiliary dimensions differ')
        except BaseException:
            self.close()
            raise

    def electron_slot(self, slot):
        bank = self.banks['electron']
        return bank['ret'][slot], bank['less'][slot]

    def gamma_slot(self):
        bank = self.banks['gamma']
        return bank['ret'][0], bank['less'][0]

    def close(self):
        for bank in self.banks.values():
            close_solution(bank)
        self.banks.clear()
        for a in self.aux.values():
            a._mmap.close()
        self.aux.clear()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


class PreviousHistory:
    """Mapped legacy original A58 preparation/controls or canonical provider."""
    def __init__(self, path):
        self.path = Path(path).resolve()
        self.binding = binding(self.path)
        self.provider = None
        self.legacy = None
        if self.path.suffix == '.npz':
            self.legacy = StoredNPZ(self.path)
            try:
                self.times = self.legacy['times']
                self.aux = {name: self.legacy[name] for name in AUXILIARY}
                nt = len(self.times)
                if (self.legacy['electron_ret'].shape != (nt, nt, 128, 12, 12) or
                        self.legacy['electron_less'].shape != (nt, nt, 128, 12, 12) or
                        self.legacy['gamma_ret'].shape != (nt, nt, 114, 114) or
                        self.legacy['gamma_less'].shape != (nt, nt, 114, 114)):
                    raise ValueError('full original previous history required')
            except BaseException:
                self.legacy.close()
                raise
        else:
            self.provider = SavedHistory(self.path)
            self.times, self.aux = self.provider.times, self.provider.aux

    def electron_slot(self, slot):
        if self.provider is not None:
            return self.provider.electron_slot(slot)
        return self.legacy['electron_ret'][:, :, slot], self.legacy['electron_less'][:, :, slot]

    def gamma_slot(self):
        if self.provider is not None:
            return self.provider.gamma_slot()
        return self.legacy['gamma_ret'], self.legacy['gamma_less']

    def close(self):
        if self.provider is not None:
            self.provider.close()
        if self.legacy is not None:
            self.legacy.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


class DressedHistory(SavedHistory):
    def __init__(self, provider_path, c, bare):
        super().__init__(provider_path)
        self.c, self.bare = c, bare
        if ((c.ns, c.m, c.d, c.f) != (128, 12, 111, 3) or
                not np.array_equal(self.times, bare.times)):
            self.close()
            raise ValueError('original material and exact bare bath times required')
        self.bindings.extend(getattr(bare, 'bindings', ()))

    def electron(self, a, b):
        bank = self.banks['electron']
        less = bank['less'][:, a, b]
        # A56 physical spectral convention. theta(0)=1 keeps raw R diagonal.
        spectral = bank['ret'][:, a, b] if a >= b else -bank['ret'][:, b, a].swapaxes(-1, -2).conj()
        return less+spectral, less

    def ordered(self, q, a, b, selection='positions'):
        if int(q) != 0:
            return self.bare.ordered(q, a, b, selection)
        if selection not in ('positions', 'ionic', 'field'):
            raise ValueError('dressed bank has positions only, never reconstructed phase')
        ret, less = self.gamma_slot()
        spectral = ret[a, b] if a >= b else -ret[b, a].T.conj()
        indices = np.arange(111) if selection == 'ionic' else np.arange(111, 114) if selection == 'field' else np.arange(114)
        selection_ix = np.ix_(indices, indices)
        return (less[a, b]+spectral)[selection_ix], less[a, b][selection_ix]

    def boson(self, q, a, b):
        return self.ordered(q, a, b, 'ionic')

    def close(self):
        super().close()
        bare = getattr(self, 'bare', None)
        if bare is not None:
            if hasattr(bare, 'close'):
                bare.close()
            else:
                # GaussianTransitions is the original raw prescribed bath and
                # owns only its bounded cache, not a file or U history.
                bare.cache.clear()
                bare.cache_size = 0
                bare._transitions = None
            self.bare = None


def open_history(provider_path, c, bare):
    """Take ownership of bare and all mappings until provider.close()."""
    return DressedHistory(provider_path, c, bare)
