#!/usr/bin/env python3
"""Independent exact-polynomial and SciPy angular oracle; no material covariance claim."""
import argparse
from hashlib import sha256
import itertools
import json
import math
import os
from pathlib import Path
import sys
from types import SimpleNamespace as NS

for _name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT/'src'))
import numpy as np
import scipy
from scipy.special import lpmv
import sympy as sp
from threadpoolctl import threadpool_limits
from substrate_framework.material_angular_intertwiners import (
    harmonic_intertwiner, original_channel_intertwiner, intertwiner_diagnostics)
from substrate_framework.qe_material_field_derivatives import solid_harmonics, centered_derivatives
from substrate_framework.qe_nc_projectors import qe_real_harmonics, centered_projectors


def binding(path):
    path = Path(path).resolve()
    data = path.read_bytes()
    return dict(path=str(path), bytes=str(len(data)), sha256=sha256(data).hexdigest())


def authenticate(entry):
    actual = binding(entry['path'])
    if actual['sha256'] != entry['sha256'] or int(actual['bytes']) != int(entry['bytes']):
        raise ValueError('original byte binding mismatch: '+str(entry['path']))
    return actual


def registration():
    record = json.loads((HERE/'angular-intertwiner-registration.json').read_text())
    for relative, digest in record['source_authority'].items():
        if binding(ROOT/relative)['sha256'] != digest:
            raise ValueError('registered source changed: '+relative)
    return record


def sources(record):
    return [binding(ROOT/path) for path in sorted(record['source_authority'])] + [
        binding(HERE/'angular-intertwiner-registration.json')]


def publish(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')


def sphere_quadrature():
    """Independent Gauss-Legendre x uniform-azimuth rule, exact through l=4."""
    z, wz = np.polynomial.legendre.leggauss(8)
    phi = (np.arange(32)+.317)*2*np.pi/32
    zz, pp = np.meshgrid(z, phi, indexing='ij')
    radius = np.sqrt(1-zz*zz)
    directions = np.column_stack([(radius*np.cos(pp)).ravel(),
                                  (radius*np.sin(pp)).ravel(), zz.ravel()])
    return directions, np.repeat(wz, len(phi))*2*np.pi/len(phi)


def independent_harmonics(momenta):
    """SciPy associated Legendre definition, not production Cartesian tensors."""
    p = np.asarray(momenta)
    radius = np.linalg.norm(p, axis=1)
    if np.any(radius == 0):
        raise ValueError('independent angular oracle requires nonzero directions')
    z = p[:, 2]/radius
    phi = np.arctan2(p[:, 1], p[:, 0])
    columns = []
    for ell in range(3):
        for m in range(ell+1):
            factor = np.sqrt((2*ell+1)/(4*np.pi)*math.factorial(ell-m)/math.factorial(ell+m))
            base = factor*lpmv(m, ell, z)
            if m == 0:
                columns.append(base)
            else:
                columns.extend((np.sqrt(2)*base*np.cos(m*phi), np.sqrt(2)*base*np.sin(m*phi)))
    return np.column_stack(columns)


def exact_blocks(rotation):
    """Solve literal polynomial coefficients symbolically; verify exact identity."""
    x, y, z = sp.symbols('x y z', real=True)
    root3 = sp.sqrt(3)
    bases = ((sp.Integer(1),), (z, -x, -y),
             ((2*z*z-x*x-y*y)/2, -root3*x*z, -root3*y*z,
              root3*(x*x-y*y)/2, root3*x*y))
    transformed = sp.Matrix(rotation.tolist())*sp.Matrix([x, y, z])
    replacements = dict(zip((x, y, z), transformed))
    result = []
    for ell, basis in enumerate(bases):
        columns = []
        for polynomial in basis:
            target = sp.expand(polynomial.subs(replacements, simultaneous=True))
            p = sp.Poly(target, x, y, z)
            c = p.coeff_monomial
            if ell == 0:
                coefficients = [target]
            elif ell == 1:
                coefficients = [c(z), -c(x), -c(y)]
            else:
                coefficients = [c(z*z), -c(x*z)/root3, -c(y*z)/root3,
                                (c(x*x)-c(y*y))/root3, c(x*y)/root3]
            if sp.expand(target-sum(a*b for a, b in zip(coefficients, basis))) != 0:
                raise AssertionError('independent symbolic polynomial identity failed')
            columns.append(sp.Matrix(coefficients))
        q = sp.Matrix.hstack(*columns)
        if sp.simplify(q.T*q-sp.eye(len(basis))) != sp.zeros(len(basis)):
            raise AssertionError('exact symbolic orthogonality failed')
        result.append(q)
    return result


def rotations():
    for permutation in itertools.permutations(range(3)):
        for signs in itertools.product((-1, 1), repeat=3):
            yield np.diag(signs) @ np.eye(3, dtype=int)[list(permutation)]


def channel_fixture():
    # Deliberately interleaved l, repeated radial channels, unequal radial values.
    angular = np.array([2, 0, 1, 2, 1, 0], dtype=int)
    radial, ell, lm = [], [], []
    for i, l in enumerate(angular):
        for j in range(2*l+1):
            radial.append(i); ell.append(l); lm.append(l*l+j)
    return NS(angular_l=angular, channel_radial=np.array(radial), channel_l=np.array(ell),
              channel_lm=np.array(lm), nh=len(radial))


class AnalyticRadialFixture:
    """Oracle-only unequal analytic radial functions; never material input."""
    def __init__(self, upf):
        self.tables = NS(upf=upf)
        self.scale = np.arange(1, len(upf.angular_l)+1)*.19

    def evaluate(self, s):
        a = self.scale[None, :]
        value = (1+3*a)*np.exp(-s[:, None]*a)
        return value, -a*value, a*a*value

    def beta(self, q):
        return self.evaluate(q*q)[0]*q[:, None]**self.tables.upf.angular_l


def residual(left, right):
    return float(np.linalg.norm(left-right)/max(1., float(np.linalg.norm(left))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cpu', required=True, type=int, choices=range(6, 11))
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    os.sched_setaffinity(0, {args.cpu})
    record = registration()
    if args.output.resolve().is_relative_to(ROOT):
        raise ValueError('oracle output must be outside Git/repository')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with threadpool_limits(limits=1):
        direction, weight = sphere_quadrature()
        y = independent_harmonics(direction)
        upf = channel_fixture()
        radial = AnalyticRadialFixture(upf)
        tables = NS(upf=upf, beta=radial.beta)
        p = direction*np.linspace(.2, 1.7, len(direction))[:, None]
        e, f = np.array([.31, -.27, .81]), np.array([-.4, .73, .19])
        b = centered_derivatives(radial, p, e, f)
        expected_b = (independent_harmonics(p)[:, upf.channel_lm]
                      * np.linalg.norm(p, axis=1)[:, None]**upf.channel_l
                      * radial.evaluate(np.sum(p*p, axis=1))[0][:, upf.channel_radial]
                      * (-1j)**upf.channel_l)
        metrics = dict(polynomial=0., quadrature=0., angular=0., solid=0.,
                       centered=residual(b[0], expected_b), derivatives=0., nc_centered=0.,
                       orthogonality=0., composition=0.)
        sensitivity = dict(transpose=0., improper_parity=0., channel_order=0.,
                           repeated_radial_swap=0., cross_channel=0., identity=0.)
        for r in rotations():
            exact = exact_blocks(r)
            q = original_channel_intertwiner(upf, r)
            expected = np.zeros_like(q)
            first = 0
            for ell in upf.angular_l:
                last = first+2*int(ell)+1
                expected[first:last, first:last] = np.array(exact[int(ell)]).astype(float)
                first = last
            metrics['polynomial'] = max(metrics['polynomial'], residual(expected, q))
            yr = independent_harmonics(direction @ r.T)
            for ell in range(3):
                sl = slice(ell*ell, (ell+1)**2)
                block = harmonic_intertwiner(r, ell)
                integrated = y[:, sl].T @ (weight[:, None]*yr[:, sl])
                metrics['quadrature'] = max(metrics['quadrature'], residual(block, integrated))
                metrics['angular'] = max(metrics['angular'], residual(yr[:, sl], y[:, sl] @ block))
                values = solid_harmonics(p, e, f)[0][:, sl]
                rotated = solid_harmonics(p @ r.T, e, f)[0][:, sl]
                metrics['solid'] = max(metrics['solid'], residual(rotated, values @ block))
            br = centered_derivatives(radial, p @ r.T, r @ e, r @ f)
            metrics['derivatives'] = max(metrics['derivatives'], *(residual(vr, v @ q) for vr, v in zip(br, b)))
            metrics['nc_centered'] = max(metrics['nc_centered'],
                residual(centered_projectors(tables, p @ r.T), centered_projectors(tables, p) @ q))
            metrics['orthogonality'] = max(metrics['orthogonality'], residual(np.eye(upf.nh), q.T @ q))
            # A noncommuting fixed operation checks the right-action order.
            s = np.array([[0, 1, 0], [0, 0, -1], [1, 0, 0]])
            metrics['composition'] = max(metrics['composition'], residual(
                original_channel_intertwiner(upf, s @ r), q @ original_channel_intertwiner(upf, s)))
            wrong = dict(transpose=q.T, identity=np.eye(upf.nh))
            ordering = np.arange(upf.nh)
            ordering[[0, 2]] = ordering[[2, 0]]
            wrong['channel_order'] = q[np.ix_(ordering, ordering)]
            repeated = np.arange(upf.nh)
            repeated[:5], repeated[9:14] = repeated[9:14].copy(), repeated[:5].copy()
            wrong['repeated_radial_swap'] = q[:, repeated]
            cross = q.copy(); cross[0, 5] += .25
            wrong['cross_channel'] = cross
            if round(np.linalg.det(r)) == -1:
                wrong['improper_parity'] = original_channel_intertwiner(upf, -r)
            for name, candidate in wrong.items():
                sensitivity[name] = max(sensitivity[name], residual(br[0], b[0] @ candidate))
        metrics['qe_against_legendre'] = residual(qe_real_harmonics(direction, 2), y)
        # An arbitrary complex non-Hermitian D must retain its real defect.
        r = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
        q = original_channel_intertwiner(upf, r)
        d = np.diag(np.arange(upf.nh)+1).astype(complex)
        d[0, 7], d[6, 2] = .37+.81j, -.29j
        before = d.copy()
        report = intertwiner_diagnostics(q, d)
        d_expected = np.linalg.norm(q @ d @ q.T-d)
        checks = {name: value < 3e-12 for name, value in metrics.items()}
        checks.update({'sensitive_'+name: value > 1e-3 for name, value in sensitivity.items()})
        checks['raw_complex_D_preserved'] = bool(np.array_equal(d, before) and
            abs(report['d_conjugacy_frobenius']-d_expected) < 1e-12 and d_expected > 1)
        checks['raw_D_not_hermitized'] = report['d_source_hermiticity_frobenius'] > .1
        rejected = {}
        for name, call in (
            ('noncubic', lambda: harmonic_intertwiner(np.diag([1, 1, 2]), 1)),
            ('unsupported_l', lambda: harmonic_intertwiner(np.eye(3), 3)),
            ('reordered_ledger', lambda: original_channel_intertwiner(
                NS(**{**vars(upf), 'channel_lm': upf.channel_lm[::-1]}), r))):
            try:
                call()
            except ValueError:
                rejected[name] = True
            else:
                rejected[name] = False
        checks.update({'reject_'+k: v for k, v in rejected.items()})
    result = dict(schema='P254-angular-intertwiner-oracle-v1',
        status='passed' if all(checks.values()) else 'failed', checks=checks,
        relative_residuals=metrics, mutation_residuals=sensitivity,
        exact_polynomial_identities='all 48 signed permutations, l=0,1,2, independent SymPy coefficients',
        quadrature='8-point Gauss-Legendre by 32 shifted azimuth points; SciPy associated Legendre',
        d_diagnostics=report, source_bindings=sources(record),
        versions=dict(numpy=np.__version__, scipy=scipy.__version__, sympy=sp.__version__),
        cpu_affinity=sorted(os.sched_getaffinity(0)),
        scope='small independent angular oracle only; not saved-output material or complete operator covariance')
    publish(args.output, result)
    if not all(checks.values()):
        raise SystemExit('angular oracle failed; inspect saved diagnostics')


if __name__ == '__main__':
    main()
