#!/usr/bin/env python3
"""Execute canonical scalar/Fourier checks of physical coordinate variations.

This is an independent finite mathematical program, not a Native trajectory,
material energy result or a new equilibrium/finite-amplitude licence.
"""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import time
from types import SimpleNamespace as NS

import numpy as np
import sympy as sp

from mixed_coordinate_response import coordinate_rhs


def fixture():
    rat = lambda n, d=17: sp.Rational(n, d)
    matrix = lambda rows, cols: np.array(
        [[rat(1+i*cols+j) for j in range(cols)] for i in range(rows)], object)
    c0 = matrix(6, 6)
    c0 = c0+c0.T+np.eye(6, dtype=int)*3
    contact = matrix(3, 3)
    contact = contact+contact.T
    reference = NS(ncells=5, masses=np.array([rat(i+3, 1) for i in range(6)], object))
    coefficients = NS(delta_curvature=np.eye(6, dtype=int)*rat(2),
                      mixed=matrix(6, 3), contact=contact)
    polarizations = np.array([[0, 1, 0], [0, 0, -1]], int)
    charges = np.array([rat(1, 1), rat(2, 1)], object)
    return reference, coefficients, NS(mass=rat(7, 1)), charges, polarizations, c0


def jets(length, offset):
    mean = sp.Matrix([sp.Rational(offset+i, 13) for i in range(length)])
    plus = sp.Matrix([sp.Rational(offset+2*i, 19)+sp.I*sp.Rational(i+1, 23)
                      for i in range(length)])
    return dict(mean=mean, plus=plus, minus=sp.conjugate(plus),
                mixed=sp.Matrix([sp.Rational(offset+3*i+1, 29) for i in range(length)]))


def canonical_oracle():
    r, c, f, z, e, c0 = fixture()
    symbols = {name: sp.Matrix(sp.symbols(name+'0:'+str(n), real=True))
               for name, n in [('u', 6), ('p', 6), ('a', 2), ('pi', 2), ('g', 6), ('ell', 3)]}
    u, p, a, pi, g, ell = (symbols[k] for k in ('u', 'p', 'a', 'pi', 'g', 'ell'))
    charge = sp.Matrix(np.repeat(z, 3)[:, None]*np.tile(np.eye(3, dtype=int), (2, 1)))
    shift = sp.Matrix(e).T*a
    kinetic = p-charge*shift
    fixed_force = sp.Matrix([sp.Rational(i+1, 11) for i in range(6)])
    energy = r.ncells*((u.T*(g-fixed_force))[0]
        +(u.T*sp.Matrix(c0+c.delta_curvature)*u)[0]/2
        +sum(kinetic[i]**2/(2*r.masses[i]) for i in range(6))
        +(shift.T*ell)[0]+(shift.T*sp.Matrix(c.contact)*shift)[0]/2
        +(u.T*sp.Matrix(c.mixed)*shift)[0])+(pi.T*pi)[0]/(2*f.mass)
    vectors = dict(u=energy.diff(p)/r.ncells, p=-energy.diff(u)/r.ncells,
                   coordinate=energy.diff(pi), momentum=-energy.diff(a))
    families = {name: jets(len(vector), 2+i) for i, (name, vector) in enumerate(symbols.items())}
    lam, eta = sp.symbols('lambda eta')
    substitutions = {}
    for name, vector in symbols.items():
        family = families[name]
        for i, symbol in enumerate(vector):
            substitutions[symbol] = (family['mean'][i]+lam*family['plus'][i]
                                     +eta*family['minus'][i]+lam*eta*family['mixed'][i])
    results = []
    outputs = {}
    for direction in ('plus', 'mixed'):
        small = NS(**{field: np.array(families[name][direction], object).reshape(-1)
                      for field, name in [('u', 'u'), ('p', 'p'), ('coordinate', 'a'), ('momentum', 'pi')]})
        actual = coordinate_rhs(reference=r, coefficients=c, field=f, charges=z,
            polarizations=e, curvature=c0, gamma=True, small=small,
            gradient=np.array(families['g'][direction], object).reshape(-1),
            field_gradient=np.array(families['ell'][direction], object).reshape(-1))
        expected = {}
        for name, vector in vectors.items():
            expanded = vector.subs(substitutions, simultaneous=True)
            derivative = expanded.diff(lam) if direction == 'plus' else expanded.diff(lam).diff(eta)
            expected[name] = np.array(derivative.subs({lam: 0, eta: 0}), object).reshape(-1)
        defect = max(abs(complex(a-b)) for name in expected for a, b in zip(actual[name], expected[name]))
        results.append(dict(name='canonical_scalar_'+direction, passed=defect < 1e-11, max_abs_defect=defect))
        outputs[direction] = actual
    return results, outputs, r


def fourier_oracle():
    r, c, f, z, e, c0 = fixture()
    c0, delta = np.asarray(c0, float), np.asarray(c.delta_curvature, float)
    b = np.arange(36, dtype=float).reshape(6, 6)/47
    mode, n = 1, r.ncells
    phase = np.exp(2j*np.pi*mode*np.arange(n)/n)
    q = 2*np.pi*mode/n
    cq = c0+b*np.exp(-1j*q)+b.T*np.exp(1j*q)
    up = np.arange(1, 7)/17+1j*np.arange(6)/23
    pp, gp = up*1.3, up*(.7-.2j)
    small = NS(u=up, p=pp, coordinate=np.zeros(0, complex), momentum=np.zeros(0, complex))
    numeric_r = NS(ncells=n, masses=np.asarray(r.masses, float))
    numeric_c = NS(delta_curvature=delta, mixed=np.asarray(c.mixed, float), contact=np.asarray(c.contact, float))
    actual = coordinate_rhs(reference=numeric_r, coefficients=numeric_c, field=NS(mass=float(f.mass)),
        charges=np.asarray(z, float), polarizations=e, curvature=cq, gamma=False,
        small=small, gradient=gp, field_gradient=np.zeros(3, complex))
    ur = phase[:, None]*up
    force = -(ur @ (c0+delta).T+np.roll(ur, 1, axis=0) @ b.T
              +np.roll(ur, -1, axis=0) @ b)-phase[:, None]*gp
    expected_force = np.sum(phase.conj()[:, None]*force, axis=0)/n
    defect = float(np.max(np.abs(actual['p']-expected_force)))
    return dict(name='independent_periodic_force', passed=defect < 1e-11,
                max_abs_defect=defect, q_curvature_transpose_gap=float(np.max(np.abs(cq-cq.T)))), actual, expected_force, cq, up


def run():
    checks, outputs, reference = canonical_oracle()
    fourier, output, expected, cq, up = fourier_oracle()
    checks.append(fourier)
    first = outputs['plus']
    gaps = dict(global_photon_missing_N=max(abs(complex(v-v/reference.ncells)) for v in first['momentum']),
                unphysical_matter_N=max(abs(complex(v-v*reference.ncells)) for v in first['u']),
                half_mixed=max(abs(complex(v-v/2)) for v in outputs['mixed']['p']),
                transpose_q_force=float(np.max(np.abs((cq-cq.T) @ up))),
                real_projection=max(abs(complex(v).imag) for v in first['p']))
    mutations = [dict(name=name, detected=bool(np.isfinite(gap) and gap > 1e-9), max_abs_defect=gap)
                 for name, gap in gaps.items()]
    return dict(passed=all(x['passed'] for x in checks) and all(x['detected'] for x in mutations),
                checks=checks, mutations=mutations)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    start, cpu = time.monotonic(), time.process_time()
    root = Path(__file__).resolve().parent
    sources = []
    for name in ('mixed_coordinate_response.py', Path(__file__).name):
        path = root/name
        raw = path.read_bytes()
        sources.append(dict(path=str(path), bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()))
    report = dict(schema='P254-physical-coordinate-response-oracle-v1', source_bindings=sources,
                  scope='Canonical variation block only; not complete coupled propagation or material energy.')
    try:
        report.update(run())
    except Exception as error:
        report.update(passed=False, error=dict(type=type(error).__name__, message=str(error)))
    report['execution'] = dict(elapsed_seconds=time.monotonic()-start, cpu_seconds=time.process_time()-cpu,
                              process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(passed=report['passed'], output=str(args.output))))
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
