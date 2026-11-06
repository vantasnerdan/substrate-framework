"""Build the same-functional Gamma/nonzero-q static bridge for P254."""
import json
from pathlib import Path
import time
from decimal import Decimal
import xml.etree.ElementTree as ET

import numpy as np
import sympy as sp
from catdft.qe_phonon_io import parse_dynamical_xml
from substrate_framework.periodic_curvature import DIRECTIONS, scalar_hessian, force_hessian
from native_cases import HERE, REG, WORK, prepare, run_case, save, setup


def stencil_identity():
    x, y, h = sp.symbols('x y h', real=True, nonzero=True)
    a, b, c, d, e, f = sp.symbols('a b c d e f', real=True)
    functional = a + b*x + c*y + (d*x*x + 2*e*x*y + f*y*y)/2
    pairs = {name: sp.expand(functional.subs({x: h*u, y: h*v}))
             + sp.expand(functional.subs({x: -h*u, y: -h*v}))
             for name, (u, v) in DIRECTIONS.items()}
    actual = sp.Matrix([[(pairs['e1']-2*a)/h**2,
                        (pairs['sum']-pairs['e1']-pairs['e2']+2*a)/(2*h**2)],
                       [(pairs['sum']-pairs['e1']-pairs['e2']+2*a)/(2*h**2),
                        (pairs['e2']-2*a)/h**2]])
    assert sp.simplify(actual - sp.hessian(functional, (x, y))) == sp.zeros(2)
    return {'quadratic_hessian': str(sp.simplify(actual)),
            'scope': 'Exact stencil identity, including nonstationary linear terms; native accuracy is separate.'}


def target(q, primitive):
    parsed = parse_dynamical_xml(REG['inputs'][q.replace('-', '_') + '_dyn'], primitive)
    record = parsed['q_star'][0]['force_matrix_ry_per_bohr2']
    full = np.array(record['real']) + 1j*np.array(record['imag'])
    block = full.reshape(6, 6)[np.ix_([0, 4], [0, 4])]
    if np.max(np.abs(block.imag)) > 1e-10:
        raise RuntimeError('Selected self-conjugate-q target has unresolved imaginary curvature')
    return block.real


def baseline_consistency(gamma, half):
    old_path = Path(REG['inputs']['primitive_pw_in']).parent / 'data-file-schema.xml'
    old = ET.fromstring(old_path.read_bytes())
    old_energy = Decimal(old.findtext('output/total_energy/etot').strip()) * 2
    checks = {
        'new_gamma_minus_saved_ry': float(Decimal(gamma['energy_ry_per_primitive'])-old_energy),
        'scope': 'Gamma reference matching' if half is None else 'Gamma and folded-supercell matching'}
    checks['matched'] = abs(checks['new_gamma_minus_saved_ry']) <= 1e-7
    if half is not None:
        checks.update(
            folded_minus_gamma_ry=float(Decimal(half['energy_ry_per_primitive'])
                                       -Decimal(gamma['energy_ry_per_primitive'])),
            folded_minus_gamma_magnetization=half['magnetization_per_primitive']
                                             -gamma['magnetization_per_primitive'])
        checks['matched'] = (checks['matched'] and abs(checks['folded_minus_gamma_ry']) <= 1e-7
                            and abs(checks['folded_minus_gamma_magnetization']) <= 1e-3)
    save(WORK / 'baseline-consistency.json', checks)
    if not checks['matched']:
        raise RuntimeError('Repair branch/folding before displacement wave: ' + json.dumps(checks))
    return checks


def curvature_rows(template, primitive, bindings, q, center):
    reference = target(q, primitive)
    allowance = max(1e-4, .01 * float(np.max(np.abs(reference))))
    seed = Path(center['xml']).parent / 'scratch/scf.save'
    rows = []
    for step in REG['numerical_design']['steps_bohr']:
        energies = {'center': center['energy_ry_per_primitive']}
        forces, cases = {}, []
        for direction, vector in DIRECTIONS.items():
            for sign, suffix in ((1, '+'), (-1, '-')):
                name = f'h{step:g}-{direction}-' + ('plus' if sign == 1 else 'minus')
                path, pattern = prepare(template, primitive, bindings, q, name,
                                       np.array(vector)*sign*step, seed)
                result = run_case(path, pattern)
                energies[direction+suffix] = result['energy_ry_per_primitive']
                forces[direction+suffix] = result['conjugate_force_ry_per_bohr']
                cases.append(result)
        scalar, force = scalar_hessian(energies, step), force_hessian(forces, step)
        errors = {'scalar_vs_DFPT': float(np.max(np.abs(scalar-reference))),
                  'force_vs_DFPT': float(np.max(np.abs(force-reference))),
                  'scalar_vs_force': float(np.max(np.abs(scalar-force)))}
        row = {'step_bohr': step, 'scalar_hessian': scalar.tolist(),
               'force_hessian_raw': force.tolist(), 'DFPT_hessian': reference.tolist(),
               'allowance_ry_per_bohr2': allowance, 'errors': errors,
               'agrees': all(value <= allowance for value in errors.values()), 'cases': cases}
        rows.append(row)
        save(WORK / q / 'curvature-progress.json', rows)
        print('ACTION14_CURVATURE ' + q + ' ' + json.dumps({k:v for k,v in row.items() if k!='cases'}), flush=True)
    rows[-1]['change_from_coarser_scalar'] = (
        np.array(rows[-1]['scalar_hessian'])-rows[0]['scalar_hessian']).tolist()
    return rows


def main():
    started = time.time()
    analytic = stencil_identity()
    template, primitive, bindings = setup()
    # Independently useful Gamma work does not wait for supercell initialization.
    import os
    phase = os.environ.get('ACTION14_PHASE', 'all')
    if phase not in ('all', 'gamma', 'half-b1'):
        raise ValueError('Unregistered action phase')
    selected = ('gamma', 'half-b1') if phase == 'all' else (phase,)
    centers = {}
    if phase == 'half-b1':
        centers['gamma'] = json.loads((WORK / 'gamma/center/action-observables.json').read_text())
    for q in selected:
        if q == 'gamma':
            seed, center_name = REG['inputs']['existing_primitive_seed'], 'center'
        else:
            initializer = os.environ.get('ACTION14_INITIALIZER_SPEC')
            if not initializer:
                raise ValueError('Nonzero-q continuation requires its registered derived initializer')
            seed = json.loads(Path(initializer).read_text())
            center_name = 'center-folded-density-001'
        path, signs = prepare(template, primitive, bindings, q, center_name, [0, 0], seed)
        centers[q] = run_case(path, signs)
    baseline = baseline_consistency(centers['gamma'], centers.get('half-b1'))
    rows = {q: curvature_rows(template, primitive, bindings, q, centers[q]) for q in selected}
    if phase == 'half-b1':
        completed_gamma = json.loads((WORK / 'gamma-result.json').read_text())
        rows['gamma'] = completed_gamma['curvatures']['gamma']
    complete = set(rows) == {'gamma', 'half-b1'}
    result = {'status': 'native_static_bridge_constructed' if all(r['agrees'] for v in rows.values() for r in v)
              else 'native_bridge_requires_numerical_or_mapping_repair',
              'phase': phase, 'full_registered_bridge_complete': complete,
              'elapsed_seconds': time.time()-started, 'analytic_stencil': analytic,
              'baseline_consistency': baseline, 'centers': centers, 'curvatures': rows,
              'scientific_scope': 'Actual nonstationary small SG15/FD spin-PBE constrained static action '
                'at the q members explicitly present in curvatures. Not Ni31H6 calibration, '
                'projected dynamic HXC, physical loss, coherence domain or nuclear rates.',
              'next_construction': 'Transfer the native static action and actual spatial spin-PBE kernel '
                'into the retarded coupled response; repair discrepancies rather than manufacture rejection.'}
    save(WORK / ('result.json' if phase == 'all' else phase + '-result.json'), result)
    if complete and phase != 'all':
        save(WORK / 'result.json', result)
    print('ACTION14_COMPLETE ' + phase + ' ' + result['status'], flush=True)


if __name__ == '__main__':
    main()
