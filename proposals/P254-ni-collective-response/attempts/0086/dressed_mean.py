"""Authenticated all31 mean/U/full-phase S driven by the actual G1 history.

Only the new dressed descriptor is admitted. Original bare preparation is a
separate identity, never an alias for the actual dressed forcing input. Saved
operations are inherited unchanged from A63; no material reload or integration.
"""
from pathlib import Path
import json
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent / n) for n in (
    '0074', '0063', '0062', '0060', '0057', '0055', '0052', '0051',
    '0050', '0047', '0046', '0044', '0045', '0041', '0030', '0029')]
from history_storage import StoredNPZ, authenticate, binding, finite_blocks
from saved_mean import SavedMean
from history_provider import SCHEMA as PROVIDER_SCHEMA
from fullwindow_force_cache import actual_nodes, provider_origin, SCHEMA as CACHE_SCHEMA

SCHEMA = 'p254-dressed-fullhistory-common-mean-reference-v1'
REFERENCE = 'common-cubic-mean-curvature-at-held-input-rho-and-covariance'


class DressedMean(SavedMean):
    """Context-owned saved reference with A63's electron_slot/gamma_slot API.

    The four node_* arrays retain raw G1 forcing. In particular an unchanged
    CommonMeanFunctional can consume scalar(node_rho[a], mean[a],
    node_cubic_covariance_force[a]), curvature(mean[a]), and
    mean_rhs(mean[a], node_electronic_gradient[a],
             node_cubic_covariance_force[a]). No baseline catalog is invented.
    """

    def __init__(self, result_path):
        self.path = Path(result_path).resolve()
        self.spec = json.loads(self.path.read_bytes())
        s = self.spec
        if (s['schema'] != SCHEMA or s['status'] != 'complete' or
                s['reference_source_kind'] != REFERENCE):
            raise ValueError('complete genuine dressed fullhistory common mean required')
        self.is_bare = False
        self.is_full_mean = True
        self.is_dressed_mean = True
        self.binding = binding(self.path)
        self.input_identity = s['input_identity']
        self.reference_identity = self.input_identity['numerical_reference']
        self.preparation_origin = s['preparation_origin']
        reader = s['producer_reader']
        if (any(reader['original'][k] != reader['preserved'][k]
                for k in ('bytes', 'sha256')) or
                reader['preserved'] not in s['source_bindings']):
            raise ValueError('executed dressed reader snapshot binding differs')
        inputs = [self.input_identity['actual_dataset'], self.reference_identity,
            self.preparation_origin, s['coefficient_manifest'], s['force_cache'],
            s['compact_source'], s['compact_source_arrays'],
            *s['physical_history_bindings'], *s['source_bindings'], s['archive'],
            reader['preserved']]
        for item in inputs:
            authenticate(item)
        self.bindings = [self.binding, *inputs]
        original = json.loads(Path(self.reference_identity['path']).read_bytes())
        force = json.loads(Path(s['force_cache']['path']).read_bytes())
        actual = json.loads(Path(self.input_identity['actual_dataset']['path']).read_bytes())
        if (original['schema'] != 'p254-a52-requested-time-reference-v1' or
                original['coefficient_manifest'] != s['coefficient_manifest'] or
                original['times_inverse_ry'] != s['times_inverse_ry'] or
                force['schema'] != CACHE_SCHEMA or force['status'] != 'complete' or
                force['original_input_dataset'] != self.preparation_origin or
                force['coefficient_manifest'] != s['coefficient_manifest'] or
                force['compact_source'] != s['compact_source'] or
                force['compact_source_arrays'] != s['compact_source_arrays'] or
                force['times_inverse_ry'] != s['times_inverse_ry'] or
                force['completed_original_internal_q'] != list(range(64)) or
                actual['schema'] != PROVIDER_SCHEMA or actual['status'] != 'complete' or
                actual['times_inverse_ry'] != s['times_inverse_ry'] or
                actual['input_identity']['numerical_reference'] != self.reference_identity or
                actual['input_identity']['source_artifact'] != actual['input_mean'] or
                actual['input_identity']['previous_artifact'] != actual['input_previous']):
            raise ValueError('dressed numerical/material/cache/G1 ancestry differs')
        provider_origin(self.input_identity['actual_dataset']['path'], self.preparation_origin)
        ancestry = [actual['input_previous'], actual['input_mean'], actual['input_kernels'],
            force['archive'], *force['input_bindings'], *original['physical_history_bindings']]
        for item in ancestry:
            authenticate(item)
            if item not in s['physical_history_bindings']:
                raise ValueError('actual ancestry or original physical history binding omitted')
        self.bindings.extend(ancestry)
        n_descriptor = actual['auxiliary']['n']
        if (n_descriptor['shape'] != [] or np.dtype(n_descriptor['dtype']).kind not in 'iu'):
            raise ValueError('G1 normalization must retain its exact integer representation')
        work = s['force_work']
        if (work['nonGamma_contractions'] != 0 or work['Gamma_contractions'] != 31 or
                work['completed_cache_internal_q'] != list(range(64))):
            raise ValueError('unchanged all64 nonGamma cache and all31 actual Gamma nodes required')
        self.archive = StoredNPZ(s['archive']['path'])
        try:
            self.times = self.archive['times']
            self.u = self.archive['U']
            self.s = self.archive['S_gamma']
            self.mean = self.archive['mean_gamma']
            self.rho0 = self.archive['rho0']
            self.w0 = self.archive['W0_gamma']
            self.n = self.archive['n']
            self.node_rho = self.archive['node_rho']
            self.node_ordered_gamma_covariance = self.archive['node_ordered_gamma_covariance']
            self.node_cubic_covariance_force = self.archive['node_cubic_covariance_force']
            self.node_electronic_gradient = self.archive['node_electronic_gradient']
            layouts = ((self.u, (31, 128, 12, 12)),
                (self.s, (31, 228, 228)), (self.mean, (31, 228)),
                (self.rho0, (128, 12, 12)), (self.w0, (228, 228)),
                (self.node_rho, (31, 128, 12, 12)),
                (self.node_ordered_gamma_covariance, (31, 114, 114)),
                (self.node_cubic_covariance_force, (31, 111)),
                (self.node_electronic_gradient, (31, 228)))
            if (self.times.dtype != np.dtype('<f8') or
                    not np.array_equal(self.times, np.arange(31, dtype='<f8')) or
                    self.times.tolist() != s['times_inverse_ry'] or
                    any(a.shape != shape for a, shape in layouts) or
                    self.n.shape != () or self.n.dtype.kind not in 'iu' or int(self.n) <= 0):
                raise ValueError('exact full31/all128/12/full-phase228 saved axes and integer N required')
            complex_arrays = (self.u, self.rho0, self.w0, self.node_rho,
                self.node_ordered_gamma_covariance)
            flexible_arrays = (self.s, self.mean, self.node_cubic_covariance_force,
                self.node_electronic_gradient)
            if (any(a.dtype != np.dtype('<c16') for a in complex_arrays) or
                    any(a.dtype not in (np.dtype('<f8'), np.dtype('<c16'))
                        for a in flexible_arrays) or
                    not all(finite_blocks(a) for a in (self.times, *(a for a, _ in layouts)))):
                raise ValueError('finite original binary64 dressed mean/reference arrays required')
            if (not np.array_equal(self.u[0], np.broadcast_to(np.eye(12), (128, 12, 12))) or
                    not np.array_equal(self.s[0], np.eye(228)) or
                    not np.array_equal(self.mean[0], np.zeros(228))):
                raise ValueError('literal original U0/S0/mean0 required')
            with SavedMean(self.preparation_origin['path']) as bare:
                if (not bare.is_bare or bare.reference_identity != self.reference_identity or
                        not np.array_equal(self.times, bare.times) or
                        not np.array_equal(self.rho0, bare.rho0) or
                        not np.array_equal(self.w0, bare.w0) or
                        int(self.n) != int(bare.archive['n'])):
                    raise ValueError('original bare rho0/W0/exact integer N preparation differs')
            rho, covariance, actual_bindings = actual_nodes(
                self.input_identity['actual_dataset']['path'], self.preparation_origin, self.times)
            if (not np.array_equal(self.node_rho, rho) or
                    not np.array_equal(self.node_ordered_gamma_covariance, covariance)):
                raise ValueError('persisted forcing differs from actual raw G1 density/covariance')
            for item in actual_bindings:
                if item not in s['physical_history_bindings']:
                    raise ValueError('actual G1 history/array binding omitted')
            self.bindings.extend(actual_bindings)
        except BaseException:
            self.archive.close()
            raise
