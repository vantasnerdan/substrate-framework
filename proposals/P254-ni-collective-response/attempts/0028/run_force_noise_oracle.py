#!/usr/bin/env python3
"""Independent new force-sign/order/state composition oracle, not old-suite replay.

SymPy differentiates a literal finite Hermitian H(x,y). Exact CAR traces and
independent dense Fock traces confront the production consumer, including an
initially empty third mode. Only the admitted A25 occupation-bit lift and
exterior-minor Fock-state primitives are reused; its oracle suite is not run.
All fixture parameters are dimensionless synthetic algebra, not Ni/H inputs.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
for directory in (ROOT/'src', HERE.parent/'0025', HERE.parent/'0024'):
    sys.path.insert(0, str(directory))


def binding(path):
    path = Path(path).resolve()
    raw = path.read_bytes()
    return dict(path=str(path), bytes=str(len(raw)), sha256=hashlib.sha256(raw).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--result', type=Path, required=True)
    args = parser.parse_args()
    began = time.monotonic()
    report = dict(status='failed', checks=[], mutations=[], scope=__doc__)
    try:
        import numpy as np
        import sympy as sp
        import scipy
        import run_physical_factor_oracle as admitted
        from substrate_framework import material_force_fluctuations as production
        from substrate_framework import physical_factor_preparation as factors
        registration = json.loads((HERE/'force-noise-registration.json').read_bytes())
        expected = {item['path']: item for item in registration['read_source_bindings']}
        for module in (admitted, factors):
            record = binding(module.__file__)
            if record != expected[record['path']]:
                raise ValueError('admitted primitive bytes changed')
        report['sources'] = [binding(path) for path in (Path(__file__), production.__file__,
            factors.__file__, admitted.__file__, HERE/'force-noise-registration.json',
            np.__file__, sp.__file__, scipy.__file__)]
        report['versions'] = dict(numpy=np.__version__, sympy=sp.__version__, scipy=scipy.__version__)
        p, q, x, y = sp.symbols('p q x y', real=True)
        gx = sp.Matrix([[2, 1+sp.I, 2], [1-sp.I, -1, sp.I], [2, -sp.I, 3]])
        gy = sp.Matrix([[1, 2*sp.I, 1+sp.I], [-2*sp.I, 3, 2], [1-sp.I, 2, -2]])
        h = sp.diag(1, 4, 7) + x*gx + y*gy
        forces = [-h.diff(x), -h.diff(y), sp.eye(3)]
        assert all(matrix == matrix.H for matrix in forces)
        # Exact integer CAR matrices are supplied by the admitted primitive,
        # not by the production Wick/factor expression being tested.
        units = {}
        for i in range(3):
            for j in range(3):
                eij = np.zeros((3, 3))
                eij[i, j] = 1
                units[i, j] = sp.Matrix(admitted.second_quantized(eij).real.astype(int))
        def lift(matrix):
            return sum((matrix[i, j]*units[i, j] for i in range(3) for j in range(3)), sp.zeros(8))
        weights = [((p if bits & 1 else 1-p)*(q if bits & 2 else 1-q)
                    if not bits & 4 else 0) for bits in range(8)]
        many_state = sp.diag(*weights)
        lifted = [lift(matrix) for matrix in forces]
        symbolic_means = sp.Matrix([sp.simplify(sp.trace(many_state*a)) for a in lifted])
        symbolic_c = sp.Matrix(3, 3, lambda i, j: sp.simplify(
            sp.trace(many_state*lifted[i]*lifted[j])-symbolic_means[i]*symbolic_means[j]))
        r = sp.Matrix([[2, 1+sp.I], [0, 3]])
        psi_s = r.col_join(sp.zeros(1, 2))
        d_s = r.inv()*sp.diag(p, q)*r.inv().H
        s_s = psi_s.H*psi_s
        m_s = [psi_s.H*a*psi_s for a in forces]
        t_s = [[(a.H*psi_s).H*(b*psi_s) for b in forces] for a in forces]
        factor_c = sp.Matrix(3, 3, lambda i, j: sp.simplify(
            sp.trace(d_s*t_s[i][j])-sp.trace(d_s*m_s[i]*d_s*m_s[j])))
        assert all(sp.simplify(v) == 0 for v in factor_c-symbolic_c)
        assert all(sp.simplify(sp.trace(d_s*m_s[i])-symbolic_means[i]) == 0 for i in range(3))
        lost_s = sp.Matrix(3, 3, lambda i, j: sp.simplify(
            sp.trace(d_s*(t_s[i][j]-m_s[i]*s_s.inv()*m_s[j]))))
        assert any(v != 0 for v in lost_s)
        assert sp.simplify(symbolic_c[0, 1]-symbolic_c[1, 0]) != 0
        assert symbolic_means[0] != 0 and symbolic_c[0, 2] != 0
        report['symbolic'] = dict(hamiltonian=str(h), means=str(symbolic_means),
            connected=str(symbolic_c), empty_complement=str(lost_s),
            identity='exact SymPy equality to finite CAR before numerical consumer calls')
        substitution = {p: sp.Rational(1, 3), q: sp.Rational(3, 5)}
        def numeric(value):
            return np.array(value.subs(substitution).evalf(), dtype=np.complex128)
        psi, d = numeric(psi_s), numeric(d_s)
        operators = [numeric(a) for a in forces]
        rho_fock = admitted.fock_state(np.eye(3), np.array([1/3, 3/5, 0.]))
        ops_fock = [admitted.second_quantized(a) for a in operators]
        means = np.array([np.trace(rho_fock @ a) for a in ops_fock])
        covariance = np.array([[np.trace(rho_fock @ a @ b)-means[i]*means[j]
                              for j, b in enumerate(ops_fock)] for i, a in enumerate(ops_fock)])
        second = covariance + means[:, None]*means[None, :]
        def pack(value):
            value = np.asarray(value)
            return dict(real=value.real.tolist(), imag=value.imag.tolist())
        def compare(name, actual, expected_value):
            actual, expected_value = np.asarray(actual), np.asarray(expected_value)
            error = float(np.max(np.abs(actual-expected_value), initial=0))
            allowance = 1e-12 + 1e-10*float(np.max(np.abs(expected_value), initial=0))
            report['checks'].append(dict(name=name, max_abs_error=error, allowance=allowance,
                actual=pack(actual), expected=pack(expected_value)))
            if error > allowance:
                raise AssertionError(name)
        def exposed(name, actual, mutation):
            gap = float(np.max(np.abs(np.asarray(actual)-mutation), initial=0))
            threshold = 1e-8*(1+float(np.max(np.abs(actual), initial=0)))
            report['mutations'].append(dict(name=name, max_abs_gap=gap, threshold=threshold,
                actual=pack(actual), mutated=pack(mutation)))
            if gap <= threshold:
                raise AssertionError('unexposed mutation: '+name)
        compare('independent dense versus exact CAR means', means, numeric(symbolic_means).ravel())
        compare('independent dense versus exact CAR directed covariance', covariance, numeric(symbolic_c))
        from types import SimpleNamespace
        class Coupling:
            npw, ndof = 3, 6
            @property
            def ionic(self):
                return self
            def prepare(self, direction):
                flat = direction.reshape(-1)
                gradient = flat[0]*numeric(gx) + flat[4]*numeric(gy)
                return SimpleNamespace(apply=lambda vectors: gradient @ vectors)
        coupling = Coupling()
        directions = {'H32-x': np.array([1., 0., 0., 0., 0., 0.]),
                      'Ni1-y': np.array([0., 0., 0., 0., 1., 0.])}
        kwargs = dict(factor_identity='synthetic nonorthogonal raw frame',
            pw_basis_identity='synthetic three full modes',
            state_authority='explicit synthetic p=1/3,q=3/5,empty third mode',
            operator_authority='literal SymPy derivative of synthetic Hermitian H(x,y)',
            coordinate_identity='synthetic two Cartesian coordinates', row_block_size=2)
        psi_before, d_before = psi.copy(), d.copy()
        result = production.initial_force_fluctuations(psi, d, coupling, directions, block_size=1, **kwargs)
        arrays, metadata = result.archive_payload()
        compare('production means include negative ionic sign', arrays['means'], means)
        compare('full ordered connected AB', arrays['ordered_connected'], covariance)
        compare('full ordered connected BA', arrays['reverse_ordered_connected'], covariance.T)
        compare('unconnected ordered AB', arrays['ordered_second_moment'], second)
        compare('unconnected ordered BA', arrays['reverse_ordered_second_moment'], second.T)
        compare('symmetrized noise', arrays['symmetrized_connected'], (covariance+covariance.T)/2)
        compare('expected CAR commutator', arrays['commutator'], covariance-covariance.T)
        compare('retained empty complement contribution', arrays['empty_complement_contribution'], numeric(lost_s))
        for i, left in enumerate(metadata['operator_order']):
            compare('raw overlap '+left, result.preparation.functional.overlaps[left], numeric(m_s[i]))
            for j, right in enumerate(metadata['operator_order']):
                compare('full ordered product '+left+' '+right,
                        result.preparation.functional.ordered_products[left, right], numeric(t_s[i][j]))
        if not np.array_equal(psi, psi_before) or not np.array_equal(d, d_before):
            raise AssertionError('consumer mutated its input state/frame')
        restored = production.initial_force_fluctuations(psi, d, coupling, directions,
            gram=numeric(s_s), gram_authority='exact synthetic same-frame Gram', block_size=2, **kwargs)
        restored_arrays, _ = restored.archive_payload()
        compare('honest authenticated-Gram constructor path', restored_arrays['ordered_connected'], covariance)
        compare('different column block partition', restored_arrays['means'], means)
        projected_only = arrays['ordered_connected']-arrays['empty_complement_contribution']
        exposed('discarding empty complement', covariance, projected_only)
        exposed('reversing AB and BA', covariance, covariance.T)
        class WrongSign(Coupling):
            def prepare(self, direction):
                correct = super().prepare(direction)
                return SimpleNamespace(apply=lambda vectors: -correct.apply(vectors))
        wrong, _ = production.initial_force_fluctuations(psi, d, WrongSign(), directions,
                                                         block_size=2, **kwargs).archive_payload()
        exposed('wrong force sign in means', means[:2], wrong['means'][:2])
        exposed('wrong force sign in force-number covariance', covariance[:2, 2], wrong['ordered_connected'][:2, 2])
        compare('force-force alone is blind to global sign', wrong['ordered_connected'][:2, :2], covariance[:2, :2])
        changed, _ = production.initial_force_fluctuations(psi, .8*d, coupling, directions,
                                                          block_size=2, **kwargs).archive_payload()
        exposed('changed density coefficients', covariance, changed['ordered_connected'])
        normalized = psi/np.linalg.norm(psi, axis=0)
        changed_frame, _ = production.initial_force_fluctuations(normalized, d, coupling, directions,
                                                                block_size=2, **kwargs).archive_payload()
        exposed('renormalized raw columns with untransported D', covariance, changed_frame['ordered_connected'])
        try:
            production.initial_force_fluctuations(psi, d, coupling,
                {'invalid-complex-direction': np.array([1., 0., 0., 0., 1j, 0.])}, block_size=1, **kwargs)
        except ValueError:
            report['checks'].append(dict(name='complex non-Hermitian direction refused', passed=True))
        else:
            raise AssertionError('complex force direction was silently accepted')
        report['status'] = 'passed'
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error), traceback=traceback.format_exc())
    report['elapsed_seconds'] = time.monotonic()-began
    with args.result.open('x') as stream:
        json.dump(report, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(status=report['status'], result=str(args.result), elapsed_seconds=report['elapsed_seconds'])))
    return 0 if report['status'] == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
