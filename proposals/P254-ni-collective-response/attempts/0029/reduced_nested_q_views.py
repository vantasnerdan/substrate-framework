"""One authenticated maximum q bank, with bounded in-memory nested selections.

Loader cutover (inside ``load_projected_coefficients``)::

    from reduced_nested_q_views import NestedQBanks
    q_provider = NestedQBanks(path, manifest, bands, features,
                             QCoefficients, q_feature_geometry)

The existing loader still authenticates the specification/implementation/original
boundary and loads fixed/passive/mechanical data for the *requested* bands. This
provider only replaces its q-provider construction; the numeric QCoefficients
contract, original per-edge NPZs and their bindings do not change.

Disk: only MappedQBanks' manifest-maximum b12-f27 directory is assembled, shared
by all four nested requests. Its original receipt authenticates the maximum
arrays and binds all 128 original edge archives plus the original kernel. Every
source is authenticated once per provider lifetime even when that receipt was
already present. No reduced q arrays or extra scientific receipts are written;
old reduced scratch is neither read nor deleted here.

Memory: ionic band truncations are basic-slice views. With the maximum feature
family, density/potential band truncations and the entire kernel are also views;
the maximum bands/family request returns the original mapped record itself. A
non-strided feature subset allocates only its final requested density, potential
and kernel arrays, filled from basic-slice runs (no full-size gather temporary).
The provider retains Gamma and one requested non-Gamma, in addition to the
underlying mapper's Gamma and one maximum non-Gamma mapping. It drops the old
requested non-Gamma before constructing another subset. Returned records may
outlive eviction if callers retain them; the caller then owns that extra
lifetime. Mappings share the same disk payload across provider instances and
processes, but separate reduced providers have separate subset buffers.

Concurrency: a per-provider RLock serializes cache mutation; a per-q Linux/POSIX
flock serializes maximum assembly, receipt publication and opening/authentication
across cooperating readers. Locks are released by the OS on process exit; an
unreceipted interrupted assembly is rebuilt by MappedQBanks. All assemblers of
this directory must use this factory: the unchanged legacy MappedQBanks factory
does not take this lock. There is no worker/monitor/background framework.

Selections use the actual q and minus-q geometry, including both spin halves and
variable self-inverse Nyquist closure. This module does no projection, kernel
application, rescaling, conjugation, thresholding, stiffness change or equation
work: all physical Omega/NFFT and Fourier/spin conventions remain in the source
coefficients and their numeric consumers.
"""

import fcntl
from threading import RLock

import numpy as np

from reduced_projection_storage import MappedQBanks


def _feature_selection(full_g, requested_g):
    """Locate actual per-q Gs; a nominal family size is not a feature stride."""
    lookup = {tuple(g): i for i, g in enumerate(full_g)}
    try:
        within = np.asarray([lookup[tuple(g)] for g in requested_g], dtype=int)
    except KeyError as error:
        raise ValueError('requested q geometry is not nested in the maximum bank') from error
    return np.concatenate((within, len(full_g) + within))


def _regular_slice(selection):
    """Return a true basic slice only when every selected index follows it."""
    if len(selection) == 1:
        return slice(int(selection[0]), int(selection[0]) + 1)
    step = int(selection[1] - selection[0])
    if step > 0 and np.all(np.diff(selection) == step):
        return slice(int(selection[0]), int(selection[-1]) + 1, step)
    return None


def _runs(selection):
    """Destination/source slices for consecutive runs, preserving exact order."""
    start = 0
    for end in range(1, len(selection) + 1):
        if end == len(selection) or selection[end] != selection[end - 1] + 1:
            yield (slice(start, end),
                   slice(int(selection[start]), int(selection[end - 1]) + 1))
            start = end


def _subset_axis(array, selection, axis):
    regular = _regular_slice(selection)
    if regular is not None:
        key = [slice(None)] * array.ndim
        key[axis] = regular
        return array[tuple(key)]
    shape = list(array.shape)
    shape[axis] = len(selection)
    result = np.empty(shape, dtype=array.dtype)
    for destination, source in _runs(selection):
        dst = [slice(None)] * array.ndim
        src = [slice(None)] * array.ndim
        dst[axis] = destination
        src[axis] = source
        np.copyto(result[tuple(dst)], array[tuple(src)])
    result.flags.writeable = False
    return result


def _subset_kernel(kernel, selection):
    regular = _regular_slice(selection)
    if regular is not None:
        return kernel[regular, regular]
    result = np.empty((len(selection), len(selection)), dtype=kernel.dtype)
    runs = tuple(_runs(selection))
    for dst_row, src_row in runs:
        for dst_column, src_column in runs:
            np.copyto(result[dst_row, dst_column], kernel[src_row, src_column])
    result.flags.writeable = False
    return result


class NestedQBanks:
    """Callable QCoefficients provider over one manifest-maximum MappedQBanks.

    Only the original nested 8/12-band, 7/27-family contract is accepted. The
    supplied manifest must contain the canonical maximum 12/27 export, although
    it may be Gamma-complete: absent q data retains the mapper's explicit error.
    ``geometry`` and ``q_type`` are the existing loader's physical geometry and
    numeric QCoefficients constructor, not alternate physics implementations.
    """

    def __init__(self, manifest_path, manifest, bands, family, q_type, geometry):
        spec = manifest['specification']
        if spec['bands'] != 12 or spec['features'] != 27:
            raise ValueError('nested q views require the canonical maximum 12-band/27-family export')
        if bands not in (8, 12) or family not in (7, 27):
            raise ValueError('requested nested basis/features not present')
        self.bands = bands
        self.family = family
        self.q_type = q_type
        self.geometry = geometry
        self._minus = spec['indexing']['minus']
        self._maximum = MappedQBanks(manifest_path, manifest, spec['bands'],
                                     spec['features'], q_type, geometry)
        self._lock = RLock()
        self._gamma = None
        self._recent_index = None
        self._recent = None

    def __call__(self, iq):
        if not 0 <= iq < 64:
            raise IndexError(iq)
        with self._lock:
            if iq == 0 and self._gamma is not None:
                return self._gamma
            if iq == self._recent_index:
                return self._recent
            if iq != 0:
                # Do not retain a previous allocated subset while filling its successor.
                self._recent_index = None
                self._recent = None
            with (self._maximum.root / ('q%02d.lock' % iq)).open('a+b') as lock:
                fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
                try:
                    # MappedQBanks checks receipt source bindings and authenticates
                    # assembled arrays. It only opens sources during assembly, so
                    # authenticate their original bytes explicitly on reuse too.
                    keys = ['edge/%02d/%03d' % (iq, s) for s in range(128)]
                    keys.append('kernel/%02d' % iq)
                    items = self._maximum.manifest['items']
                    if any(key not in items for key in keys):
                        raise ValueError('physical q bank is incomplete: q=' + str(iq))
                    for key in keys:
                        self._maximum.items.authenticate(items[key]['archive'])
                    maximum = self._maximum(iq)
                finally:
                    fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
            value = self._nested(iq, maximum)
            if iq == 0:
                self._gamma = value
            else:
                self._recent_index = iq
                self._recent = value
            return value

    def _nested(self, iq, maximum):
        m = self.bands
        if m == self._maximum.bands and self.family == self._maximum.family:
            return maximum
        ionic = maximum.ionic[:, :, :m, :m]
        reverse = maximum.ionic_reverse[:, :, :m, :m]
        density = maximum.density[:, :, :m, :m]
        potential = maximum.potential[:, :m, :m, :]
        if self.family == self._maximum.family:
            return self.q_type(ionic=ionic, ionic_reverse=reverse, density=density,
                               potential=potential, kernel=maximum.kernel,
                               feature_partner=maximum.feature_partner,
                               feature_g=maximum.feature_g,
                               feature_spin=maximum.feature_spin)
        gs, partner = self.geometry(iq, self.family)
        selection = _feature_selection(maximum.feature_g, gs)
        # The partner indexes the *minus-q requested list*, never this q's max
        # list. Check exact restriction without opening/storing the minus-q bank.
        opposite = int(self._minus[iq])
        opposite_full_g, _ = self.geometry(opposite, self._maximum.family)
        opposite_g, _ = self.geometry(opposite, self.family)
        opposite_selection = _feature_selection(opposite_full_g, opposite_g)
        if not np.array_equal(maximum.feature_partner[selection],
                              opposite_selection[partner]):
            raise ValueError('nested q selection does not preserve the physical partner geometry')
        spin = maximum.feature_spin[selection]
        for array in (gs, partner, spin):
            array.flags.writeable = False
        return self.q_type(ionic=ionic, ionic_reverse=reverse,
                           density=_subset_axis(density, selection, 1),
                           potential=_subset_axis(potential, selection, 3),
                           kernel=_subset_kernel(maximum.kernel, selection),
                           feature_partner=partner, feature_g=gs, feature_spin=spin)
