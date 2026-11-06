"""Compact six-quark T1,I=0 carrier versus two independent canonical P2 nucleons.

A source-connected finite-channel merger calculation, not a nucleus or NN
potential. Imports the unchanged canonical bouquet operator and literal CAR.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import platform
import time
from pathlib import Path
from typing import Any

import numpy as np
import scipy
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh
from threadpoolctl import threadpool_limits

from rich_journal_flux_carrier_20261019 import (
    Carrier, canonical_degree_two_loop, canonical_exterior_maps, matter_law,
)
from rich_nucleon_hierarchy_20261019 import matter_sector

HERE = Path(__file__).resolve().parent
PRIOR = HERE/'rich_journal_flux_carrier_20261019_result.json'


def six_quark_matter() -> dict[str, Any]:
    """Full 12-mode six-quark CAR, exact I=0, colour-singlet S=1 seed at Sz=1."""
    data, levels = matter_sector(6)
    states = [sum(1 << mode for mode in occupied)
              for occupied in itertools.combinations(range(12), 6)]
    selected = [j for j, bits in enumerate(states)
                if sum(1 if mode//6 == 0 else -1
                       for mode in range(12) if bits >> mode & 1) == 0]
    iso2 = sum((a@a for a in data['isospin']),
               sparse.csr_matrix((len(states), len(states)), dtype=complex))
    values, vectors = eigh(iso2[np.ix_(selected, selected)].toarray())
    singlets = vectors[:, abs(values) < 1e-9]
    embedding = np.zeros((len(states), singlets.shape[1]), complex)
    embedding[selected] = singlets
    reduce = lambda operator: embedding.conj().T@(operator@embedding)
    source = next(level['psi'] for level in levels if level['S'] == 1. and level['I'] == 0.)
    # The shared even-N supplier returns Sz=0. Raising selects the actual
    # vector polarization Sz=+1; no spinor C4 convention is carried over.
    source = (data['spin'][0]+1j*data['spin'][1])@source
    source /= np.linalg.norm(source)
    seed = embedding.conj().T@source
    g2 = reduce(data['g2'])
    lam = float(np.vdot(seed, g2@seed).real)
    return {'G': [reduce(a) for a in data['mixed']],
            'Q': [reduce(a) for a in data['colour']],
            'S': [reduce(a) for a in data['spin']],
            'charge': reduce(sparse.diags(data['charge'])),
            'seed': seed, 'dimension': embedding.shape[1], 'I': 0.,
            'source_G2': lam,
            'source_G2_eigen_residual': float(np.linalg.norm(g2@seed-lam*seed)),
            'full_CAR_dimension': len(states), 'fixed_I3_dimension': len(selected)}


class CompactCarrier(Carrier):
    """Only the integer-spin polarization and channel observables differ."""
    def polarized(self, x: Any) -> Any:
        part = x
        result = x.copy()
        for turn in range(1, 4):
            part = self.rotate(part)
            result += (1j**turn)*part
        return result/4

    def vector_representation(self, x: Any) -> dict[str, Any]:
        result = super().cubic_representation(x)
        central = np.zeros_like(x)
        for dest, inverse, spin in self.quarter_turns:
            central += (x@spin.T)[inverse]+(x@spin.conj())[dest]
        result.update({'six_quarter_turn_class_sum_target': 2.,
                       'class_sum_eigen_residual': float(np.linalg.norm(central-2*x)),
                       'interpretation': 'T1 vector triplet, not continuous SO3'})
        return result

    def solve_compact(self, mass: float, electric: float) -> dict[str, Any]:
        started = time.monotonic()
        self.calls = 0
        h = self.operator(mass, electric)
        values, vectors = eigsh(h, k=1, which='SA', v0=self.seed.ravel(),
                               tol=2e-10, maxiter=1600, ncv=24)
        energy = float(values[0])
        x = self.polarized(vectors[:, 0].reshape(self.n, self.d))
        x /= np.linalg.norm(x)
        self.last_state, self.last_energy, self.last_hamiltonian = x, energy, h
        weights = np.sum(abs(x)**2, axis=1)
        degree = np.array([sum(self.loops['degree'][v] for v in s) for s in self.states])
        interaction = -self.interaction(x)/(2*mass)
        electric_part = 4*electric*self.electric_C2[:, None]*x
        gauss = []
        for a in range(8):
            image = sum((self.field_operator(self.loops['Q'][a], k)@x for k in range(3)),
                        np.zeros_like(x))+x@self.matter['Q'][a].T
            gauss.append(float(np.linalg.norm(image)))
        source = self.interaction(self.seed)
        source /= np.linalg.norm(source)
        z0 = float(abs(np.vdot(self.seed[0], x[0]))**2)
        w1 = float(abs(np.vdot(source, x))**2)
        gap = 16*electric/3
        g = -1/(2*np.sqrt(3)*mass)
        lam = self.matter['source_G2']
        return {'quark_number': 6, 'source_S': 1., 'I': 0., 'source_Sz': 1.,
                'mass_m': mass, 'electric_kE': electric, 'a': 1.,
                'normalized_source_g': g, 'one_source_gap_Omega': gap,
                'field_regulator': 'canonical total representation degree <=2 per occupied cell',
                'field_dimension': self.n, 'matter_dimension': self.d,
                'tensor_dimension': self.dimension,
                'full_CAR_dimension': self.matter['full_CAR_dimension'],
                'fixed_I3_dimension': self.matter['fixed_I3_dimension'],
                'source_G2': lam,
                'source_G2_eigen_residual': self.matter['source_G2_eigen_residual'],
                'energy': energy,
                'electric_energy_expectation': float(np.vdot(x, electric_part).real),
                'Pauli_energy_expectation': float(np.vdot(x, interaction).real),
                'vacuum_seed_overlap_Z0': z0, 'one_source_projection_weight': w1,
                'excluded_from_vacuum_plus_one_source_weight': 1-z0-w1,
                'probability_by_total_representation_degree':
                    [float(weights[degree == j].sum()) for j in range(3)],
                'probability_by_active_loop_count':
                    [float(weights[self.count == j].sum()) for j in range(3)],
                'charge': float(np.vdot(x, x@self.matter['charge'].T).real),
                'cubic_representation': self.vector_representation(x),
                'root_Gauss_max_norm': max(gauss),
                'unprojected_H_eigen_residual': float(np.linalg.norm(electric_part+interaction-energy*x)),
                'C4z_eigenphase_minus_i_residual': float(np.linalg.norm(self.rotate(x)+1j*x)),
                'one_source_only_energy': float((gap-np.sqrt(gap*gap+4*g*g*lam))/2),
                'leading_energy': float(-g*g*lam/gap),
                'H_matvec_calls': self.calls, 'elapsed_seconds': time.monotonic()-started}


def ground_exterior(carrier: Carrier, maps: dict[str, Any], mass: float) -> dict[str, Any]:
    """Exact Q3 H P2 ground-channel leakage; no exterior diagonalization."""
    out = np.zeros((maps['exterior_field_dimension'], carrier.d), complex)
    # The two tensor factors commute: doing the small matter product before
    # rectangular field expansion avoids a large redundant dense product.
    for b, g in zip(maps['maps'], carrier.matter['G']):
        out -= (b@(carrier.last_state@g.T))/(2*mass)
    weights = np.sum(abs(out)**2, axis=1)
    return {'H_exterior_norm2': float(weights.sum()),
            'by_field_channel': {str(label): float(weights[maps['categories'] == label].sum())
                                 for label in np.unique(maps['categories'])},
            'exterior_field_dimension': maps['exterior_field_dimension'],
            'interpretation': 'Full-H variance on the P2 eigenstate, not an energy-error bound.'}


def run() -> dict[str, Any]:
    prior = json.loads(PRIOR.read_text())
    comparisons = []
    with threadpool_limits(limits=1):
        loops = canonical_degree_two_loop()
        compact = CompactCarrier(loops, six_quark_matter(), 2)
        nucleon = Carrier(loops, matter_law(.5), 2)
        maps = canonical_exterior_maps(loops, compact)
        for mass in (6., 3.):
            six = compact.solve_compact(mass, .27)
            six['canonical_degree3_ground_exterior'] = ground_exterior(compact, maps, mass)
            if mass == 3.:
                old = next(row for row in prior['rows'] if row['I'] == .5 and row['mass_m'] == mass
                           and row['loop_representation'] == 'PW155_degree2')
                three = {'quark_number': 3, 'mass_m': mass, 'energy': old['energy'],
                         'one_source_only_energy': old['source_only_energy'],
                         'degree_probabilities': old['canonical_degree3_exterior'][
                             'ground_probability_by_total_representation_degree'],
                         'H_exterior_norm2': old['canonical_degree3_exterior'][
                             'ground_state_H_exterior_norm2'],
                         'cubic_representation': old['cubic_representation'],
                         'provenance': 'Existing same-law canonical P2 result; not reexecuted.'}
            else:
                fresh = nucleon.solve(mass, .27)
                ext = ground_exterior(nucleon, maps, mass)
                degree = np.array([sum(loops['degree'][v] for v in s) for s in nucleon.states])
                weights = np.sum(abs(nucleon.last_state)**2, axis=1)
                three = {'quark_number': 3, 'mass_m': mass, 'energy': fresh['energy'],
                         'one_source_only_energy': fresh['source_only_energy'],
                         'degree_probabilities': [float(weights[degree == j].sum()) for j in range(3)],
                         'H_exterior_norm2': ext['H_exterior_norm2'],
                         'cubic_representation': fresh['cubic_representation'],
                         'provenance': 'New same-law canonical P2 nucleon calculation.'}
            g, omega = six['normalized_source_g'], six['one_source_gap_Omega']
            merger = six['energy']-2*three['energy']
            comparisons.append({'mass_m': mass, 'compact_N6': six, 'one_N3_cell': three,
                                'separated_pn_I0_T1_energy': 2*three['energy'],
                                'compact_minus_separated': merger,
                                'one_source_only_compact_minus_separated':
                                    six['one_source_only_energy']-2*three['one_source_only_energy'],
                                'leading_compact_minus_separated': (28/3)*g*g/omega,
                                'physical_sign': 'compact channel penalized' if merger > 0 else 'compact channel favoured',
                                'separated_two_cell_H_exterior_norm2': 2*three['H_exterior_norm2']})
            print(json.dumps(comparisons[-1]), flush=True)
    sources = [Path(__file__), HERE/'rich_journal_flux_carrier_20261019.py',
               HERE/'rich_nucleon_hierarchy_20261019.py',
               HERE/'rich_composite_local_flux_binding_20261017_model.py', PRIOR]
    return {'kind': 'CANONICAL_COMPACT_SIX_QUARK_T1_I0_MERGER_CHANNEL',
            'plan': 'issue222-journal-construction-20261019',
            'scoped_item': 'journal-result',
            'plan_revision_supplied_by_parent_prefix': '74246a88',
            'runtime': {'python': platform.python_version(), 'numpy': np.__version__,
                        'scipy': scipy.__version__},
            'source_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
            'comparisons': comparisons,
            'regulator_contract': 'One canonical P2 bouquet per occupied cell; separated cells have independent regulators.',
            'conserved_quantum_number_contract':
                'Compact and separated channels both charge1,I0,T1. Separated pn isosinglet uses exactly degenerate N3 proton/neutron doublets.',
            'scope': ['Source-connected T1 N6 versus two independent Gamma6 N3 cells coupled to I0,T1.',
                      'No scalar supplier, no extra N6 coupling, and no fitted merger term.',
                      'Rest masses cancel: six quark masses versus two times three.',
                      'No absolute global N6 minimum, spatial NN force, nucleus, or saturation theorem.',
                      'Degree3 norm is exterior leakage on P2 states, not full-H energy or an energy-error bound.']}


if __name__ == '__main__':
    result = run()
    output = HERE/'rich_journal_compact_hierarchy_20261019_result.json'
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'output': str(output), 'merger_penalties': [
        {'mass': row['mass_m'], 'penalty': row['compact_minus_separated']}
        for row in result['comparisons']]}), flush=True)
