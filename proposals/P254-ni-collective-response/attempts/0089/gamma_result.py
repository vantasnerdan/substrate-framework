"""Authenticated Gamma-only result and exact-input later full-history reuse.

This module does not import either continuation driver. A Gamma result is not
an electronic history or a full G2 provider; paths in Pi publications may differ
but their actual NPY bytes, axes and dtype must be identical.
"""
from pathlib import Path
import json
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent / n) for n in ('0086', '0063', '0062', '0061')]
from history_storage import authenticate, binding, finite_blocks
from saved_mean import POSITIONS
from production_transport import load_solution, array_digest

SCHEMA = 'p254-independent-gamma-result-v1'
ORDERS = ('greater', 'lesser', 'retarded')


def source_bindings(value):
    """Collect actual bindings without changing their paths or identities."""
    found = {}
    def visit(item):
        if isinstance(item, dict):
            if {'path', 'bytes', 'sha256'} <= item.keys():
                descriptor = {k: item[k] for k in ('path', 'bytes', 'sha256')}
                if descriptor['path'] in found and found[descriptor['path']] != descriptor:
                    raise ValueError('conflicting source binding: ' + descriptor['path'])
                found[descriptor['path']] = descriptor
            for child in item.values():
                visit(child)
        elif isinstance(item, (list, tuple)):
            for child in item:
                visit(child)
    visit(value)
    return list(found.values())


def validate_inputs(catalog, mean, previous):
    if previous.provider is None:
        raise ValueError('actual previous G1 provider required')
    numerical = previous.provider.spec['input_identity']['numerical_reference']
    if (catalog['input_dataset'] != previous.binding or
            mean.input_identity['actual_dataset'] != previous.binding or
            catalog.get('numerical_reference') != numerical or
            mean.reference_identity != numerical or
            catalog.get('preparation_origin') != mean.preparation_origin):
        raise ValueError('actual G1, original numerical reference or preparation identity differs')
    authenticate(numerical)
    authenticate(mean.preparation_origin)
    if (not np.array_equal(mean.times, np.arange(31, dtype='<f8')) or
            not np.array_equal(previous.times, mean.times) or
            catalog['times_inverse_ry'] != mean.times.tolist() or
            not np.array_equal(1j * mean.rho0, previous.aux['source_electron_less0']) or
            not np.array_equal(-1j * mean.w0[np.ix_(POSITIONS, POSITIONS)].T,
                               previous.aux['source_gamma_less'][0, 0]) or
            mean.n.shape != () or mean.n.dtype.kind not in 'iu' or
            previous.aux['n'].shape != () or previous.aux['n'].dtype.kind not in 'iu' or
            int(mean.n) != 64 or int(previous.aux['n']) != 64):
        raise ValueError('literal original preparation, all31 grid and integer N64 required')
    return dict(kernel_input=catalog['input_dataset'], actual_dataset=previous.binding,
        numerical_reference=numerical, preparation_origin=mean.preparation_origin,
        source_artifact=mean.binding, previous_artifact=previous.binding)


def admitted_binary(binary):
    verification_path = HERE.parent / '0085' / 'verification-001.json'
    verification = json.loads(verification_path.read_bytes())
    actual = binding(binary)
    if (verification['verdict'] != 'finite-sensitive-same-discrete-operator-pass' or
            any(actual[k] != verification['production_binary'][k] for k in ('bytes', 'sha256'))):
        raise ValueError('admitted A85 real-axis binary required')
    return actual, binding(verification_path)


def pi_identity(kernels):
    if set(kernels.pi) != set(ORDERS) or set(kernels.catalog['pi']) != set(ORDERS):
        raise ValueError('all three actual Pi contour arrays required')
    result = {}
    for name in ORDERS:
        item, array = kernels.catalog['pi'][name], kernels.pi[name]
        authenticate(item)
        if (array.shape != (31, 31, 114, 114) or array.dtype != np.dtype('<c16') or
                item['shape'] != list(array.shape) or item['dtype'] != str(array.dtype) or
                not finite_blocks(array)):
            raise ValueError('finite all31/full114 complex128 Pi required')
        result[name] = {key: item[key] for key in ('path', 'bytes', 'sha256', 'shape', 'dtype')}
    return result


class GammaResult:
    """Context-owned, read-only actual complete Gamma bank plus typed receipts."""
    def __init__(self, result_path):
        self.path = Path(result_path).resolve()
        self.binding = binding(self.path)
        self.spec = json.loads(self.path.read_bytes())
        self.solution = None
        try:
            s = self.spec
            if (s['schema'] != SCHEMA or s['status'] != 'complete-gamma-block' or
                    s['full_G2_complete'] is not False or s['gamma_slots'] != 1 or
                    s['gamma_coordinates'] != 114 or s['n'] != 64 or
                    s['times_inverse_ry'] != list(range(31))):
                raise ValueError('complete independent Gamma block required, not full G2')
            for item in source_bindings(s):
                authenticate(item)
            actual_binary, verification = admitted_binary(s['binary']['path'])
            if actual_binary != s['binary'] or verification != s['binary_verification']:
                raise ValueError('Gamma result does not retain admitted actual binary identity')
            invocation = json.loads(authenticate(s['invocation']).read_bytes())
            if (invocation['input_identity'] != s['input_identity'] or
                    invocation['input_pi'] != s['input_pi'] or
                    invocation['pi_arrays'] != s['pi_arrays'] or
                    invocation['binary'] != s['binary'] or
                    invocation['binary_verification'] != s['binary_verification'] or
                    invocation['source_bindings'] != s['source_bindings'] or
                    invocation['order'] != s['order'] or invocation['beta'] != s['beta']):
                raise ValueError('Gamma invocation/result identities differ')
            self.complete = s['gamma']['complete']
            bank_path = authenticate(self.complete).parent
            self.solution = load_solution(bank_path)
            solution = self.solution
            if (solution['ret'].shape != (1, 31, 31, 114, 114) or
                    solution['less'].shape != solution['ret'].shape or
                    solution['ret'].dtype != np.dtype('<c16') or
                    solution['less'].dtype != np.dtype('<c16') or
                    not np.array_equal(solution['times'], np.arange(31, dtype='<f8')) or
                    solution['statistics'] != 1 or solution['order'] != s['order'] or
                    solution['beta'] != s['beta'] or s['order'] != 1 or s['beta'] != 1. or
                    not all(finite_blocks(solution[k]) for k in ('ret', 'less'))):
                raise ValueError('actual complete Gamma bank layout differs')
            complete = json.loads(authenticate(self.complete).read_bytes())
            if len(complete['slots']) != 1:
                raise ValueError('exactly one committed Gamma slot required')
            receipt = json.loads(authenticate(complete['slots'][0]).read_bytes())
            provenance = receipt['provenance']
            native = json.loads(authenticate(provenance['native']).read_bytes())
            if (receipt['slot'] != 0 or native.get('exit') != 0 or
                    provenance['sources'] != s['source_bindings'] or
                    provenance['binary'] != s['binary'] or
                    receipt['ret'] != array_digest(solution['ret'][0]) or
                    receipt['less'] != array_digest(solution['less'][0])):
                raise ValueError('Gamma committed slot/native provenance differs')
            self.provenance = provenance
            self.bindings = source_bindings([self.binding, s, complete, provenance['native']])
        except BaseException:
            self.close()
            raise

    def validate(self, kernels, mean, previous, *, binary, verification, order, beta):
        """Fail before electronic Native if any exact numerical input changed."""
        identity = validate_inputs(kernels.catalog, mean, previous)
        s = self.spec
        if (identity != s['input_identity'] or mean.binding != s['input_mean'] or
                previous.binding != s['input_previous'] or
                mean.preparation_origin != s['preparation_origin'] or
                verification != s['binary_verification'] or
                any(binary[k] != s['binary'][k] for k in ('bytes', 'sha256')) or
                order != s['order'] or beta != s['beta']):
            raise ValueError('precomputed Gamma M2/G1/preparation/binary/order/beta differs')
        current = pi_identity(kernels)
        for name in ORDERS:
            old = s['pi_arrays'][name]
            if any(current[name][key] != old[key] for key in ('bytes', 'sha256', 'shape', 'dtype')):
                raise ValueError('precomputed Gamma exact Pi bytes/axes/dtype changed: ' + name)
        source_ret, source_less = mean.gamma_slot()
        if self.provenance['components'] != [array_digest(a) for a in
                (source_ret, source_less, kernels.pi['retarded'], kernels.pi['lesser'])]:
            raise ValueError('precomputed Gamma original solve payload differs')
        if not np.array_equal(self.solution['less'][0, 0, 0],
                              previous.aux['source_gamma_less'][0, 0]):
            raise ValueError('actual Gamma bank changed original Gamma0')
        return self

    def close(self):
        if self.solution is not None:
            for key in ('ret', 'less'):
                self.solution[key]._mmap.close()
            self.solution = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
