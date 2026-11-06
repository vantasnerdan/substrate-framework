"""Literal saved common mean/bare source; no material reload or integration."""
from pathlib import Path
import json
import numpy as np
from history_storage import StoredNPZ, authenticate, binding, finite_blocks

REFERENCE = 'bare-mechanical-with-contacts'
PI_CONVENTION = 'full-common-joint-plus-cubic; represented-only-control-separate'
POSITIONS = np.r_[np.arange(111), np.arange(222, 225)]
BARE_SCHEMA = 'p254-prescribed-bare-history-source-v1'
FULL_MEAN_SCHEMA = 'p254-fullhistory-common-mean-reference-v1'
FULL_MEAN_REFERENCE = 'common-cubic-mean-curvature-at-held-input-rho-and-covariance'
LEGACY_MEAN_REFERENCE = 'bare-mechanical-with-original-contacts-plus-common-cubic-mean-curvature-at-held-dressed-rho-and-covariance; no-coupled-electronic-bubble-in-S'


class SavedMean:
    def __init__(self, result_path):
        self.path = Path(result_path).resolve()
        self.spec = json.loads(self.path.read_bytes())
        s = self.spec
        self.is_bare = s['schema'] == BARE_SCHEMA
        self.is_full_mean = s['schema'] == FULL_MEAN_SCHEMA
        expected_kinds = ((REFERENCE,) if self.is_bare else
            (FULL_MEAN_REFERENCE,) if self.is_full_mean else (REFERENCE, LEGACY_MEAN_REFERENCE))
        if (s['schema'] not in ('p254-common-mean-feedback-reference-v1', BARE_SCHEMA, FULL_MEAN_SCHEMA) or
                s['status'] != 'complete' or s['reference_source_kind'] not in expected_kinds):
            raise ValueError('complete genuine common mean or prescribed bare source required')
        self.binding = binding(self.path)
        if self.is_bare:
            inputs = [s['reference_bundle'], s['factor_result'], s['archive']]
            for item in inputs:
                authenticate(item)
            original = json.loads(Path(s['reference_bundle']['path']).read_bytes())
            paid = json.loads(Path(s['factor_result']['path']).read_bytes())
            if (original['schema'] != 'p254-a52-requested-time-reference-v1' or
                    paid['schema'] != 'p254-paid-cubic-long-history-column-v1' or
                    paid['status'] != 'complete-first-column-q1' or
                    paid['bare_source_factors'] != s['archive'] or
                    paid['times_inverse_ry'] != s['times_inverse_ry'] or
                    original['times_inverse_ry'] != s['times_inverse_ry']):
                raise ValueError('literal original bundle and paid raw source factors required')
            actual_inputs = json.loads(authenticate(paid['inputs']).read_bytes())
            if s['reference_bundle'] not in actual_inputs['reference']:
                raise ValueError('paid factors do not bind this original reference')
            inputs += [paid['inputs'], original['U'], original['coefficient_manifest'],
                       *original['physical_history_bindings'], *original['source_bindings']]
            self.reference_identity = s['reference_bundle']
            self.input_identity = dict(actual_dataset=self.binding, numerical_reference=self.reference_identity)
        elif self.is_full_mean:
            self.input_identity = s['input_identity']
            self.reference_identity = self.input_identity['numerical_reference']
            inputs = [self.input_identity['actual_dataset'], self.reference_identity,
                s['coefficient_manifest'], s['force_cache'], s['compact_source'],
                s['compact_source_arrays'], *s['physical_history_bindings'],
                *s['source_bindings'], s['archive'], s['producer_reader']['preserved']]
            for item in inputs:
                authenticate(item)
            if 'producer_result' in s:
                inputs += [s['producer_result'], s['publication_implementation']]
                executed = json.loads(authenticate(s['producer_result']).read_bytes())
                for key in ('schema', 'status', 'archive', 'input_identity', 'coefficient_manifest',
                        'physical_history_bindings', 'times_inverse_ry', 'producer_reader'):
                    if executed[key] != s[key]:
                        raise ValueError('published mean differs from executed producer')
                restored = {}
                for record in s['source_restoration']:
                    old = record['original_binding']
                    archived = record['archived_snapshot']
                    if (old not in executed['source_bindings'] or old['path'] in restored or
                            any(old[k] != archived[k] for k in ('bytes', 'sha256'))):
                        raise ValueError('published source snapshot differs from executed bytes')
                    authenticate(archived)
                    restored[old['path']] = archived
                if [restored.get(item['path'], item) for item in executed['source_bindings']] != s['source_bindings']:
                    raise ValueError('published sources are not exact executed restorations')
            original = json.loads(Path(self.reference_identity['path']).read_bytes())
            force = json.loads(Path(s['force_cache']['path']).read_bytes())
            old_reader = s['producer_reader']['original']
            preserved = s['producer_reader']['preserved']
            if (original['schema'] != 'p254-a52-requested-time-reference-v1' or
                    original['coefficient_manifest'] != s['coefficient_manifest'] or
                    original['times_inverse_ry'] != s['times_inverse_ry'] or
                    force['schema'] != 'p254-fullwindow-bare-covariance-force-cache-v1' or
                    force['status'] != 'complete' or
                    force['original_input_dataset'] != self.input_identity['actual_dataset'] or
                    force['coefficient_manifest'] != s['coefficient_manifest'] or
                    force['compact_source'] != s['compact_source'] or
                    force['compact_source_arrays'] != s['compact_source_arrays'] or
                    force['times_inverse_ry'] != s['times_inverse_ry'] or
                    any(old_reader[k] != preserved[k] for k in ('bytes', 'sha256'))):
                raise ValueError('fullhistory numerical/material/force/source identity differs')
            with SavedMean(self.input_identity['actual_dataset']['path']) as bare:
                if not bare.is_bare or bare.reference_identity != self.reference_identity:
                    raise ValueError('fullhistory original bare input ancestry differs')
                self._prepared_rho0 = np.array(bare.rho0)
                self._prepared_w0 = np.array(bare.w0)
                self._prepared_n = int(bare.archive['n'])
        else:
            inputs = [s['input_dataset'], s['baseline_catalog'], s['coefficient_manifest'],
                      *s['physical_history_bindings'], *s['source_bindings'], s['archive']]
            self.reference_identity = s['input_dataset']
            self.input_identity = dict(actual_dataset=s['input_dataset'], numerical_reference=self.reference_identity)
        for item in inputs:
            authenticate(item)
        self.bindings = [self.binding, *inputs]
        self.archive = StoredNPZ(s['archive']['path'])
        try:
            self.times = self.archive['times']
            self.u = self.archive['U']
            self.rho0 = self.archive['rho0']
            self.s = self.archive['S_gamma']
            self.w0 = self.archive['W0_gamma']
            mean = None if self.is_bare else self.archive['mean_gamma']
            self.mean = mean
            nt = len(self.times)
            if (self.times.tolist() != s['times_inverse_ry'] or
                    self.times.dtype != np.dtype('<f8') or nt < 2 or self.times[0] != 0 or
                    self.times[1] <= 0 or not np.array_equal(self.times, np.arange(nt)*self.times[1]) or
                    self.u.shape != (nt, 128, 12, 12) or self.rho0.shape != (128, 12, 12) or
                    self.s.shape != (nt, 228, 228) or self.w0.shape != (228, 228) or
                    (mean is not None and mean.shape != (nt, 228))):
                raise ValueError('exact full original128/12/full-phase228 saved nodes required')
            real_or_complex = (self.s,) if mean is None else (self.s, mean)
            values = (self.times, self.u, self.rho0, self.s, self.w0, *(() if mean is None else (mean,)))
            if (any(a.dtype != np.dtype('<c16') for a in (self.u, self.rho0, self.w0)) or
                    any(a.dtype not in (np.dtype('<f8'), np.dtype('<c16')) for a in real_or_complex)):
                raise ValueError('original binary64 saved mean/reference precision required')
            if not all(finite_blocks(a) for a in values):
                raise ValueError('nonfinite genuine saved reference')
            if (not np.array_equal(self.u[0], np.broadcast_to(np.eye(12), (128, 12, 12))) or
                    not np.array_equal(self.s[0], np.eye(228)) or
                    (mean is not None and not np.array_equal(mean[0], np.zeros(228)))):
                raise ValueError('original literal U0/S0/mean0 required')
            if self.is_bare:
                saved_u = np.load(original['U']['path'], mmap_mode='r', allow_pickle=False)
                try:
                    if not np.array_equal(self.u, saved_u):
                        raise ValueError('paid raw U differs from original reference')
                finally:
                    saved_u._mmap.close()
                n = self.archive['n']
                if n.shape != () or n.dtype.kind not in 'iu' or int(n) <= 0:
                    raise ValueError('original exact integer normalization required')
            if self.is_full_mean:
                n = self.archive['n']
                if (not np.array_equal(self.times, np.arange(31, dtype=float)) or
                        not np.array_equal(self.rho0, self._prepared_rho0) or
                        not np.array_equal(self.w0, self._prepared_w0) or
                        n.shape != () or n.dtype.kind not in 'iu' or int(n) != self._prepared_n):
                    raise ValueError('fullhistory original all31 preparation/normalization differs')
                del self._prepared_rho0, self._prepared_w0, self._prepared_n
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
