"""Literal saved common mean/bare source; no material reload or integration."""
from pathlib import Path
import json
import numpy as np
from history_storage import StoredNPZ, authenticate, binding, finite_blocks

REFERENCE = 'bare-mechanical-with-contacts'
PI_CONVENTION = 'full-common-joint-plus-cubic; represented-only-control-separate'
POSITIONS = np.r_[np.arange(111), np.arange(222, 225)]


class SavedMean:
    def __init__(self, result_path):
        self.path = Path(result_path).resolve()
        self.spec = json.loads(self.path.read_bytes())
        s = self.spec
        if (s['schema'] != 'p254-common-mean-feedback-reference-v1' or
                s['status'] != 'complete' or s['reference_source_kind'] != REFERENCE):
            raise ValueError('complete genuine common mean/bare reference required')
        inputs = [s['input_dataset'], s['baseline_catalog'], s['coefficient_manifest'],
                  *s['physical_history_bindings'], *s['source_bindings'], s['archive']]
        for item in inputs:
            authenticate(item)
        self.bindings = [binding(self.path), *inputs]
        self.archive = StoredNPZ(s['archive']['path'])
        try:
            self.times = self.archive['times']
            self.u = self.archive['U']
            self.rho0 = self.archive['rho0']
            self.s = self.archive['S_gamma']
            self.w0 = self.archive['W0_gamma']
            mean = self.archive['mean_gamma']
            nt = len(self.times)
            if (self.times.tolist() != s['times_inverse_ry'] or
                    self.times.dtype != np.dtype('<f8') or nt < 2 or self.times[0] != 0 or
                    self.times[1] <= 0 or not np.array_equal(self.times, np.arange(nt)*self.times[1]) or
                    self.u.shape != (nt, 128, 12, 12) or self.rho0.shape != (128, 12, 12) or
                    self.s.shape != (nt, 228, 228) or self.w0.shape != (228, 228) or
                    mean.shape != (nt, 228)):
                raise ValueError('exact full original128/12/full-phase228 saved nodes required')
            if (any(a.dtype != np.dtype('<c16') for a in (self.u, self.rho0, self.w0)) or
                    any(a.dtype not in (np.dtype('<f8'), np.dtype('<c16')) for a in (self.s, mean))):
                raise ValueError('original binary64 saved mean/reference precision required')
            if not all(finite_blocks(a) for a in (self.times, self.u, self.rho0, self.s, self.w0, mean)):
                raise ValueError('nonfinite genuine saved reference')
            if (not np.array_equal(self.u[0], np.broadcast_to(np.eye(12), (128, 12, 12))) or
                    not np.array_equal(self.s[0], np.eye(228)) or not np.array_equal(mean[0], np.zeros(228))):
                raise ValueError('original literal U0/S0/mean0 required')
        except BaseException:
            self.archive.close()
            raise

    def electron_slot(self, slot):
        """Only one128-sector member history exists in active Python storage."""
        if not 0 <= slot < 128:
            raise ValueError('original electron slot outside0..127')
        nt = len(self.times)
        ret = np.zeros((nt, nt, 12, 12), dtype='<c16')
        less = np.empty_like(ret)
        rho0 = self.rho0[slot]
        empty0 = np.eye(12)-rho0
        for a in range(nt):
            for b in range(nt):
                end = self.u[b, slot].conj().T
                less[a, b] = 1j*(self.u[a, slot]@rho0@end)
                if a >= b:
                    greater = -1j*(self.u[a, slot]@empty0@end)
                    ret[a, b] = greater-less[a, b]
        return ret, less

    def gamma_slot(self):
        nt = len(self.times)
        ret = np.zeros((nt, nt, 114, 114), dtype='<c16')
        less = np.empty_like(ret)
        for a in range(nt):
            for b in range(nt):
                # Both ordered contractions are evaluated independently, exactly
                # as A50 ordered_coordinates; no dagger/reverse-q inference.
                w = self.s[a]@self.w0@self.s[b].T
                reverse = self.s[b]@self.w0@self.s[a].T
                less[a, b] = -1j*reverse[np.ix_(POSITIONS, POSITIONS)].T
                if a >= b:
                    greater = -1j*w[np.ix_(POSITIONS, POSITIONS)]
                    ret[a, b] = greater-less[a, b]
        return ret, less

    def close(self):
        self.archive.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
