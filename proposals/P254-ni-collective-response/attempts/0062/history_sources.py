"""Genuine U and complete saved-G adapters, without whole-history operands."""
from pathlib import Path
import json
import shutil
import tempfile
import zipfile
import numpy as np
from history_banks import binding
from general_moving_reference import SavedMovingReference, GaussianTransitions


def checked_times(times):
    times = np.asarray(times, float)
    if times.ndim != 1 or not len(times) or not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
        raise ValueError('exact increasing finite times required')
    return times


class BareHistory(SavedMovingReference):
    """Own the saved U mapping in addition to the source's bounded S cache."""
    def close(self):
        super().close()
        if hasattr(self.u, '_mmap'):
            self.u._mmap.close()


class SavedDressed:
    """Legacy NPZ byte extraction to temporary NPY mappings, never npz[history].

    Every saved retarded/lesser pair is required. Greater uses the A56 causal
    spectral convention, retaining the raw equal-time jump. Non-Gamma remains
    the prescribed original Gaussian bath, not reconstructed from Gamma.
    """
    def __init__(self, c, dataset, coefficient_binding, times=None):
        self.c = c
        self.dataset = binding(dataset)
        self.bindings = [self.dataset, coefficient_binding]
        self.temp = tempfile.TemporaryDirectory(prefix='p254-a62-dressed-')
        self.arrays = {}
        try:
            with zipfile.ZipFile(dataset) as archive:
                for name in ('times', 'electron_ret', 'electron_less', 'gamma_ret', 'gamma_less'):
                    path = Path(self.temp.name) / (name + '.npy')
                    with archive.open(name + '.npy') as source, path.open('wb') as target:
                        shutil.copyfileobj(source, target, length=1048576)
                    self.arrays[name] = np.load(path, mmap_mode='r', allow_pickle=False)
            saved = checked_times(self.arrays['times'])
            self.times = saved if times is None else checked_times(times)
            self.indices = []
            for t in self.times:
                found = np.flatnonzero(saved == t)
                if len(found) != 1:
                    raise ValueError('requested genuine saved dressed time absent')
                self.indices.append(int(found[0]))
            nt = len(saved)
            for name, tail in (('electron_ret', (128,12,12)), ('electron_less', (128,12,12)),
                               ('gamma_ret', (114,114)), ('gamma_less', (114,114))):
                if self.arrays[name].shape != (nt,nt,*tail) or self.arrays[name].dtype != np.complex128:
                    raise ValueError('complete saved dressed axes/dtype required: ' + name)
                for a in range(nt):
                    for b in range(nt):
                        if not np.isfinite(self.arrays[name][a,b]).all():
                            raise ValueError('nonfinite saved dressed pair')
            self.bare = GaussianTransitions(c, self.times)
        except BaseException:
            self.close()
            raise

    def pair(self, prefix, a, b):
        aa, bb = self.indices[a], self.indices[b]
        ret, less = self.arrays[prefix + '_ret'], self.arrays[prefix + '_less']
        spectral = ret[aa,bb] if self.times[a] >= self.times[b] else -ret[bb,aa].swapaxes(-1,-2).conj()
        return less[aa,bb] + spectral, less[aa,bb]

    def electron(self, a, b):
        return self.pair('electron', a, b)

    def ordered(self, q, a, b, selection='positions'):
        if int(q) != 0:
            return self.bare.ordered(q, a, b, selection)
        if selection == 'phase':
            raise ValueError('saved dressed positions do not contain full phase space')
        indices = {'positions': np.arange(114), 'ionic': np.arange(111), 'field': np.arange(111,114)}[selection]
        return tuple(x[np.ix_(indices, indices)] for x in self.pair('gamma', a, b))

    def boson(self, q, a, b):
        return self.ordered(q, a, b, 'ionic')

    def close(self):
        if hasattr(self, 'bare'):
            self.bare.cache.clear()
        for array in self.arrays.values():
            array._mmap.close()
        self.arrays.clear()
        self.temp.cleanup()


def open_reference(c, coefficient_binding, *, bundle=None, dataset=None, provider=None, times=None):
    if sum(x is not None for x in (bundle, dataset, provider)) != 1:
        raise ValueError('exactly one genuine bare bundle, saved dataset or saved provider required')
    if dataset is not None:
        return SavedDressed(c, dataset, coefficient_binding, times)
    if provider is not None:
        from history_provider import open_history
        record = json.loads(Path(provider).read_bytes())
        exact = checked_times(record['times_inverse_ry'])
        if times is not None and not np.array_equal(exact, checked_times(times)):
            raise ValueError('provider grid must match exact requested times')
        bare = GaussianTransitions(c, exact)
        return open_history(provider, c, bare)
    ref = BareHistory(c, bundle, coefficient_binding, times=times)
    checked_times(ref.times)
    return ref
