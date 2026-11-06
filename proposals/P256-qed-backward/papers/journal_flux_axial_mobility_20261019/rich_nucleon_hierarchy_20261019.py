"""Local nucleon-hierarchy experiment, not a spatial baryon/nucleus model.

Compare the captured number-coupled scalar supplier with a NEW finite quantum
colour-adjoint spin-vector cell.  Primitive CAR, SU(3), spin, flavour and the
new coupling are inputs.  No photon, gluon dispersion or measured scale is
claimed.  Run from this directory with the campaign's scientific Python.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import platform
from pathlib import Path
from typing import Any

import numpy as np
import scipy
from scipy import sparse
from scipy.linalg import eigh

from rich_composite_local_flux_binding_20261017_model import Law, generators

HERE = Path(__file__).resolve().parent
INPUT = HERE / 'rich_composite_local_flux_binding_20261017_actual3_input.json'
SIGMA = [np.array([[0, 1], [1, 0]], complex),
         np.array([[0, -1j], [1j, 0]]), np.diag([1., -1.])]
COLOUR = generators()


def lift(matrix: Any, states: list[int], index: dict[int, int]) -> Any:
    """Literal shared-CAR cdag matrix c, including identical-quark signs."""
    rows, cols, values = [], [], []
    for p, q in zip(*np.nonzero(matrix)):
        p, q = int(p), int(q)
        for col, bits in enumerate(states):
            if not (bits >> q) & 1:
                continue
            sign = (-1) ** (bits & ((1 << q) - 1)).bit_count()
            changed = bits ^ (1 << q)
            if (changed >> p) & 1:
                continue
            sign *= (-1) ** (changed & ((1 << p) - 1)).bit_count()
            changed |= 1 << p
            rows.append(index[changed]); cols.append(col)
            values.append(sign * matrix[p, q])
    return sparse.csr_matrix((values, (rows, cols)), shape=(len(states),) * 2)


def casimir(operators: list[Any], size: int) -> Any:
    result = sparse.csr_matrix((size, size), dtype=complex)
    for operator in operators:
        result += operator @ operator
    return result


def field_generators() -> tuple[list[Any], list[Any]]:
    # [T_c,T_a]=i f_cab T_b; the quantum field has the adjoint representation.
    colour = []
    for c in COLOUR:
        f = np.array([[float((-2j * np.trace((c @ a - a @ c) @ b)).real)
                       for b in COLOUR] for a in COLOUR])
        colour.append(np.kron(-1j * f, np.eye(3)))
    spin = []
    for j in range(3):
        rotation = np.zeros((3, 3), complex)
        for k, ell in itertools.permutations(range(3), 2):
            epsilon = ((j - k) * (k - ell) * (ell - j)) // 2
            rotation[k, ell] = -1j * epsilon
        spin.append(np.kron(np.eye(8), rotation))
    return colour, spin


FIELD_COLOUR, FIELD_SPIN = field_generators()


def matter_sector(number: int) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    states = [sum(1 << i for i in occupied)
              for occupied in itertools.combinations(range(12), number)]
    index = {bits: i for i, bits in enumerate(states)}
    colour = [lift(np.kron(np.eye(4), t), states, index) for t in COLOUR]
    spin = [lift(np.kron(np.kron(np.eye(2), s / 2), np.eye(3)), states, index)
            for s in SIGMA]
    isospin = [lift(np.kron(np.kron(s / 2, np.eye(2)), np.eye(3)), states, index)
               for s in SIGMA]
    mixed = [lift(np.kron(np.kron(np.eye(2), s), t), states, index)
             for t in COLOUR for s in SIGMA]
    c2, s2, i2, g2 = [casimir(ops, len(states))
                      for ops in (colour, spin, isospin, mixed)]
    selected = []
    for row, bits in enumerate(states):
        occupied = [k for k in range(12) if (bits >> k) & 1]
        sz2 = sum(1 if (k // 3) % 2 == 0 else -1 for k in occupied)
        iz2 = sum(1 if k // 6 == 0 else -1 for k in occupied)
        if sz2 == number % 2 and iz2 == number % 2:
            selected.append(row)
    values, vectors = eigh(c2[np.ix_(selected, selected)].toarray())
    singlets = vectors[:, abs(values) < 1e-9]
    embedding = np.zeros((len(states), singlets.shape[1]), complex)
    embedding[selected] = singlets
    reduced_spin = embedding.conj().T @ (s2 @ embedding)
    reduced_iso = embedding.conj().T @ (i2 @ embedding)
    _, multiplets = eigh(reduced_spin + np.sqrt(2) * reduced_iso)
    charge = np.array([sum(2 / 3 if k // 6 == 0 else -1 / 3
                           for k in range(12) if (bits >> k) & 1)
                       for bits in states])
    data = {'number': number, 'colour': colour, 'spin': spin, 'isospin': isospin,
            'mixed': mixed, 'g2': g2, 'charge': charge}
    levels = []
    for column in range(multiplets.shape[1]):
        psi = embedding @ multiplets[:, column]
        spin2 = float(np.vdot(psi, s2 @ psi).real)
        iso2 = float(np.vdot(psi, i2 @ psi).real)
        lam = float(np.vdot(psi, g2 @ psi).real)
        levels.append({'S': round((np.sqrt(1 + 4 * spin2) - 1) / 2, 6),
                       'I': round((np.sqrt(1 + 4 * iso2) - 1) / 2, 6),
                       'lambda_G2': lam, 'hyperfine_per_kappa': -(lam - 4 * number) / 2,
                       'lambda_eigen_residual': float(np.linalg.norm(g2 @ psi - lam * psi)),
                       'psi': psi})
    summary = {'N': number, 'full_CAR_dimension': len(states),
               'fixed_weight_dimension': len(selected),
               'singlet_fixed_weight_dimension': singlets.shape[1],
               'levels': [{k: v for k, v in level.items() if k != 'psi'} for level in levels]}
    data['summary'] = summary
    return data, levels


def dressed(data: dict[str, Any], psi: Any, coupling: float, omega: float) -> dict[str, float]:
    """Act with the full 25*d Hamiltonian on an exact two-layer cyclic state.

    Field law is DEFINED as |0> plus 24 adjoint-vector states, not an asserted
    converged harmonic oscillator. H00=0, H11=omega, H_a0=g G_a.
    """
    images = np.array([operator @ psi for operator in data['mixed']])
    lam = float(np.vdot(images, images).real)
    energy = -2 * coupling**2 * lam / (omega + np.sqrt(omega**2 + 4 * coupling**2 * lam))
    alpha = 1 / np.sqrt(1 + coupling**2 * lam / (energy - omega)**2)
    zero = alpha * psi
    one = alpha * coupling / (energy - omega) * images
    hzero = coupling * sum((operator @ one[a] for a, operator in enumerate(data['mixed'])),
                           np.zeros_like(zero))
    hone = omega * one + coupling * np.array([operator @ zero for operator in data['mixed']])
    residual = np.sqrt(np.linalg.norm(hzero - energy * zero)**2
                       + np.linalg.norm(hone - energy * one)**2)
    gauss = []
    for matter, field in zip(data['colour'], FIELD_COLOUR):
        action = np.array([matter @ row for row in one]) + field @ one
        gauss.append(np.linalg.norm(matter @ zero)**2 + np.linalg.norm(action)**2)
    j2zero = np.zeros_like(zero); j2one = np.zeros_like(one)
    for matter, field in zip(data['spin'], FIELD_SPIN):
        jzero = matter @ zero
        jone = np.array([matter @ row for row in one]) + field @ one
        j2zero += matter @ jzero
        j2one += np.array([matter @ row for row in jone]) + field @ jone
    expected_j2 = sum(float(np.vdot(operator @ psi, operator @ psi).real)
                      for operator in data['spin'])
    spin_residual = np.sqrt(np.linalg.norm(j2zero - expected_j2 * zero)**2
                            + np.linalg.norm(j2one - expected_j2 * one)**2)
    probability = float(np.vdot(one, one).real)
    charge = float(np.vdot(zero, data['charge'] * zero).real
                   + np.vdot(one, one * data['charge'][None, :]).real)
    return {'energy': float(energy), 'field_probability': probability,
            'vacuum_source_weight': float(alpha**2), 'charge': charge,
            'total_spin2': float((np.vdot(zero, j2zero) + np.vdot(one, j2one)).real),
            'Hamiltonian_eigen_residual': float(residual),
            'Gauss_max_norm2': float(max(gauss)), 'total_spin_eigen_residual': float(spin_residual)}


def scalar_and_return() -> dict[str, Any]:
    parameters = json.loads(INPUT.read_text())['parameters']
    law = Law(parameters)
    energies = [float(e[0]) for e in law.local_energies]
    returns = []
    for number in (3, 6):
        channels = []
        for a, b in itertools.product(range(2), repeat=2):
            overlap = law.overlap[number - 1, number][a, 0] * law.overlap[1, 0][b, 0]
            gap = float(law.local_energies[number - 1][a] + law.local_energies[1][b]
                        - energies[number] - energies[0] + parameters['colour_electric'] * 4 / 3)
            strength = float(3 * (parameters['hopping']**2 + parameters['spin_orbit']**2)
                             * overlap**2)
            channels.append({'cloud_indices': [a, b], 'gap': gap,
                             'norm2_per_quark': strength})
        returns.append({'N': number, 'channels': channels,
                        'second_order_energy': -number * sum(c['norm2_per_quark'] / c['gap']
                                                            for c in channels)})
    return {'parameters': parameters, 'local_scalar_ground_energies': energies,
            'triplet_vs_three_singles': energies[3] + 2 * energies[0] - 3 * energies[1],
            'scalar_mergers': {str(n): energies[n] + (n // 3 - 1) * energies[0]
                              - (n // 3) * energies[3] for n in (6, 9, 12)},
            'one_hop_returns': returns,
            'one_hop_merger_correction': returns[1]['second_order_energy']
                                        - 2 * returns[0]['second_order_energy'],
            'one_hop_scope': 'H0=local cloud plus electric flux, empty neighbours; no gradients or sea vertex. '
                            'TdagT=(t^2+v^2)I gives no spin splitting at this order only.'}


def run() -> dict[str, Any]:
    scalar = scalar_and_return()
    sectors = {n: matter_sector(n) for n in (3, 6, 9, 12)}
    summaries = [data['summary'] for data, _ in sectors.values()]
    trials = []
    for coupling in (0., .05, .2, .5, 1., 2.):
        levels = []
        for number, (data, multiplets) in sectors.items():
            for level in multiplets:
                result = dressed(data, level['psi'], coupling, 1.)
                levels.append({'N': number, 'S': level['S'], 'I': level['I'], **result})
        nucleon = next(x for x in levels if x['N'] == 3 and x['S'] == .5)
        delta = next(x for x in levels if x['N'] == 3 and x['S'] == 1.5)
        mergers = {}
        for number in (6, 9, 12):
            best = min(x['energy'] for x in levels if x['N'] == number)
            correction = best - (number // 3) * nucleon['energy']
            mergers[str(number)] = {'adjoint_field_only': correction,
                                   'new_local_product_law_scalar_plus_adjoint':
                                   scalar['scalar_mergers'][str(number)] + correction}
        data, multiplets = sectors[3]
        proton = next(x for x in multiplets if x['S'] == .5)
        lowering = data['isospin'][0] - 1j * data['isospin'][1]
        neutron_psi = lowering @ proton['psi']
        neutron_psi /= np.linalg.norm(neutron_psi)
        neutron = dressed(data, neutron_psi, coupling, 1.)
        trials.append({'g': coupling, 'omega': 1., 'levels': levels,
                       'Delta_minus_nucleon': delta['energy'] - nucleon['energy'],
                       'neutron': neutron, 'proton_neutron_energy_difference':
                       nucleon['energy'] - neutron['energy'], 'compact_mergers': mergers})
    controls = {'no_coupling': 'g=0 rows above: no split and no induced merger penalty',
                'colour_only_field': 'C_a kills matter singlets, so induced shift is exactly zero',
                'spin_only_field': {str(s): {'lambda': 4 * s * (s + 1),
                    'energy_g0p2_omega1': (1 - np.sqrt(1 + 4 * .2**2 * 4 * s * (s + 1))) / 2}
                    for s in (.5, 1.5)}}
    sources = [Path(__file__), INPUT, HERE / 'rich_composite_local_flux_binding_20261017_model.py']
    return {'kind': 'NEW_LOCAL_QUANTUM_ADJOINT_VECTOR_HIERARCHY_NOT_SPATIAL_PARTICLES',
            'plan': 'issue222-categorical-backward-20261019',
            'captured_base': '9b0d4efd04fd165e742240d5b911011b761e5a21',
            'runtime': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__},
            'source_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
            'old_scalar_supplier': scalar, 'CAR_colour_spin_spectra': summaries,
            'new_quantum_field_trials': trials, 'adverse_controls': controls,
            'scope': ['One spatial orbital per cell, 12 CAR modes (u,d x spin2 x colour3).',
                      '25-state quantum field is the declared law, not oscillator/gluon convergence.',
                      'Scalar-plus-adjoint comparison is a NEW local product law, not old spatial-H output.',
                      'g values are exploratory model inputs, not fitted or measured couplings.',
                      'No confinement, propagating gluon/photon, spatial binding, nuclear saturation or mass prediction.',
                      'Capped local shell and primitive SU3/CAR/spin/flavour/coupling remain input debt.']}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=HERE / 'rich_nucleon_hierarchy_20261019_result.json')
    args = parser.parse_args()
    result = run()
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'output': str(args.output), 'scalar_mergers': result['old_scalar_supplier']['scalar_mergers'],
                      'trials': [{'g': t['g'], 'Delta_minus_N': t['Delta_minus_nucleon'],
                                  'mergers': t['compact_mergers']} for t in result['new_quantum_field_trials']]}, indent=2))


if __name__ == '__main__':
    main()
