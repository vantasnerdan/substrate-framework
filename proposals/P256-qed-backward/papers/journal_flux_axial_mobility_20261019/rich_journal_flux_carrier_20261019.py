"""Cubic three-loop SU(3) projected-rotor baryon, with repeated flux actions.

This is a gauge-reduced bouquet law, NOT a reduction of the old shared-edge
nine-link electric Hamiltonian. Each independent four-edge loop costs 4*kE*C2.
PW19 full products and canonical PW155 total-degree-two sources are distinguished.
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
from scipy.linalg import eigh, expm, null_space
from scipy.sparse.linalg import LinearOperator, eigsh
from threadpoolctl import threadpool_limits

from rich_composite_local_flux_binding_20261017_model import link_operators
from rich_nucleon_hierarchy_20261019 import COLOUR, SIGMA, lift, matter_sector

HERE = Path(__file__).resolve().parent


def loop_law() -> dict[str, Any]:
    """P U_ab P, electric Casimir, and root conjugation on exact PW19."""
    transport, left, right = link_operators()
    magnetic = []
    for t in COLOUR:
        raw = sum((t[b, a] * transport[a, b, False]
                   for a, b in itertools.product(range(3), repeat=2)),
                  np.zeros((19, 19), complex))
        magnetic.append(1j * (raw - raw.conj().T))
    conjugation = [l+r for l, r in zip(left, right)]
    inversion = np.zeros((19, 19))
    inversion[0, 0] = 1
    for a, b in itertools.product(range(3), repeat=2):
        inversion[10+3*b+a, 1+3*a+b] = 1
        inversion[1+3*b+a, 10+3*a+b] = 1
    covariance = 0.
    for c, t in enumerate(COLOUR):
        for a, u in enumerate(COLOUR):
            expected = sum((2*np.trace((t@u-u@t)@v)*magnetic[b]
                            for b, v in enumerate(COLOUR)), np.zeros((19, 19), complex))
            # [Q_c,B_a] = +i f_cab B_b; its contraction with G_a is a scalar.
            covariance = max(covariance, float(np.linalg.norm(
                conjugation[c]@magnetic[a]-magnetic[a]@conjugation[c]-expected)))
    gram = np.array([[np.vdot(x[:, 0], y[:, 0]) for y in magnetic] for x in magnetic])
    return {'B': magnetic, 'Q': conjugation, 'inversion': inversion,
            'degree': np.array([0]+[1]*18), 'C2': np.array([0.]+[4/3]*18),
            'name': 'PW19',
            'gram_residual': float(np.linalg.norm(gram-np.eye(8)/3)),
            'conjugation_commutator_residual': covariance,
            'inversion_odd_residual': max(float(np.linalg.norm(inversion@b@inversion+b))
                                          for b in magnetic)}


def canonical_degree_two_loop() -> dict[str, Any]:
    """Exact PW multiplication through every two-fundamental source channel."""
    symmetric = []
    for a in range(3):
        for b in range(a, 3):
            v = np.zeros((3, 3), complex)
            v[a, b] = 1 if a == b else 1/np.sqrt(2)
            v[b, a] = v[a, b]
            symmetric.append(v.ravel())
    e6 = np.array(symmetric).T
    e8 = np.array([np.sqrt(2)*t.ravel() for t in COLOUR]).T
    six = [e6.conj().T@(np.kron(t, np.eye(3))+np.kron(np.eye(3), t))@e6
           for t in COLOUR]
    eight = [e8.conj().T@(np.kron(t, np.eye(3))-np.kron(np.eye(3), t.conj()))@e8
             for t in COLOUR]
    representations = [
        ('1', [np.zeros((1, 1), complex) for _ in COLOUR], 0, 0.),
        ('3', COLOUR, 1, 4/3), ('bar3', [-t.conj() for t in COLOUR], 1, 4/3),
        ('6', six, 2, 10/3), ('bar6', [-t.conj() for t in six], 2, 10/3),
        ('8', eight, 2, 3.)]
    offsets, offset, degree, casimirs = {}, 0, [], []
    for name, generators, order, c2 in representations:
        d = len(generators[0])
        offsets[name] = offset
        offset += d*d
        degree += [order]*(d*d)
        casimirs += [c2]*(d*d)
    transport = [np.zeros((offset, offset), complex) for _ in range(9)]
    allowed = {'1': ('3',), '3': ('6', 'bar3'), 'bar3': ('1', '8'),
               '6': ('8',), 'bar6': ('bar3',), '8': ('3', 'bar6')}
    intertwining_residual = 0.
    orthonormality_residual = 0.
    for r, tr, _, _ in representations:
        dr = len(tr[0])
        parent = [np.kron(t, np.eye(dr))+np.kron(np.eye(3), u)
                  for t, u in zip(COLOUR, tr)]
        for s, ts, _, _ in representations:
            if s not in allowed[r]:
                continue
            ds = len(ts[0])
            equations = np.vstack([np.kron(np.eye(ds), a)-np.kron(b.T, np.eye(3*dr))
                                   for a, b in zip(parent, ts)])
            kernel = null_space(equations, rcond=1e-11)
            if kernel.shape[1] != 1:
                raise RuntimeError(f'Nonunique SU3 CG intertwiner {r} -> {s}: {kernel.shape}')
            cg = kernel[:, 0].reshape(3*dr, ds, order='F')*np.sqrt(ds)
            intertwining_residual = max(intertwining_residual,
                max(float(np.linalg.norm(a@cg-cg@b)) for a, b in zip(parent, ts)))
            orthonormality_residual = max(orthonormality_residual,
                float(np.linalg.norm(cg.conj().T@cg-np.eye(ds))))
            cg = cg.reshape(3, dr, ds)
            for a, b in itertools.product(range(3), repeat=2):
                block = np.sqrt(dr/ds)*np.einsum('mM,nN->MNmn', cg[a], cg[b].conj())
                transport[3*a+b][offsets[s]:offsets[s]+ds*ds,
                                 offsets[r]:offsets[r]+dr*dr] = block.reshape(ds*ds, dr*dr)
    # Remove only roundoff zeros so sparse actions preserve exact source support.
    for u in transport:
        u[abs(u) < 2e-14] = 0
    magnetic = []
    for t in COLOUR:
        raw = sum((t[b, a]*transport[3*a+b] for a, b in itertools.product(range(3), repeat=2)),
                  np.zeros((offset, offset), complex))
        magnetic.append(1j*(raw-raw.conj().T))
    conjugation = []
    for a in range(8):
        blocks = [-np.kron(g[a].T, np.eye(len(g[a])))+
                  np.kron(np.eye(len(g[a])), g[a]) for _, g, _, _ in representations]
        conjugation.append(sparse.block_diag(blocks).toarray())
    partners = {'1': '1', '3': 'bar3', 'bar3': '3', '6': 'bar6', 'bar6': '6', '8': '8'}
    inversion = np.zeros((offset, offset))
    for name, generators, _, _ in representations:
        d = len(generators[0])
        for a, b in itertools.product(range(d), repeat=2):
            inversion[offsets[partners[name]]+d*b+a, offsets[name]+d*a+b] = 1
    gram = np.array([[np.vdot(x[:, 0], y[:, 0]) for y in magnetic] for x in magnetic])
    covariance = max(float(np.linalg.norm(
        conjugation[c]@magnetic[a]-magnetic[a]@conjugation[c]-
        sum((2*np.trace((t@u-u@t)@v)*magnetic[b] for b, v in enumerate(COLOUR)),
            np.zeros((offset, offset), complex))))
        for c, t in enumerate(COLOUR) for a, u in enumerate(COLOUR))
    return {'B': magnetic, 'Q': conjugation, 'inversion': inversion,
            '_representations': representations, '_offsets': offsets,
            'degree': np.array(degree), 'C2': np.array(casimirs), 'name': 'PW155_degree2',
            'gram_residual': float(np.linalg.norm(gram-np.eye(8)/3)),
            'conjugation_commutator_residual': covariance,
            'inversion_odd_residual': max(float(np.linalg.norm(inversion@b@inversion+b))
                                         for b in magnetic),
            'CG_intertwining_residual': intertwining_residual,
            'CG_orthonormality_residual': orthonormality_residual}


def matter_law(isospin: float) -> dict[str, Any]:
    data, levels = matter_sector(3)
    states = [sum(1 << i for i in occupied)
              for occupied in itertools.combinations(range(12), 3)]
    selected = [j for j, bits in enumerate(states)
                if sum((1 if mode//6 == 0 else -1)
                       for mode in range(12) if bits >> mode & 1) == 1]
    iso2 = sum((x@x for x in data['isospin']), sparse.csr_matrix((220, 220), dtype=complex))
    values, vectors = eigh(iso2[np.ix_(selected, selected)].toarray())
    embedding = np.zeros((220, sum(abs(values-isospin*(isospin+1)) < 1e-8)), complex)
    embedding[selected] = vectors[:, abs(values-isospin*(isospin+1)) < 1e-8]
    reduce = lambda op: embedding.conj().T @ (op @ embedding)
    seed = next(x['psi'] for x in levels if x['I'] == isospin and x['S'] == isospin)
    axial = lift(np.kron(np.kron(SIGMA[2], SIGMA[2]), np.eye(3)),
                 states, {bits: j for j, bits in enumerate(states)})
    colour = [reduce(x) for x in data['colour']]
    spin = [reduce(x) for x in data['spin']]
    c2 = sum((x@x for x in colour), np.zeros((embedding.shape[1],)*2, complex))
    s2 = sum((x@x for x in spin), np.zeros_like(c2))
    _, joint = eigh(c2+np.sqrt(2)*s2)
    channels = []
    for casimir in (0., 3., 6.):
        for spin2 in (.75, 3.75):
            chosen = [j for j in range(joint.shape[1])
                      if abs(np.vdot(joint[:, j], c2@joint[:, j]).real-casimir) < 1e-8
                      and abs(np.vdot(joint[:, j], s2@joint[:, j]).real-spin2) < 1e-8]
            if chosen:
                channels.append((casimir, spin2, joint[:, chosen]))
    return {'G': [reduce(x) for x in data['mixed']],
            'Q': colour, 'S': spin, 'A': reduce(axial), 'channels': channels,
            'charge': reduce(sparse.diags(data['charge'])),
            'seed': embedding.conj().T @ seed,
            'embedding': embedding, 'axial_full': axial,
            'dimension': embedding.shape[1], 'I': isospin}


class Carrier:
    def __init__(self, loops: dict[str, Any], matter: dict[str, Any], cutoff: int):
        self.loops, self.matter, self.cutoff = loops, matter, cutoff
        self.states = [(a, b, c) for a in range(len(loops['degree']))
                       for b in range(len(loops['degree']))
                       if loops['degree'][a]+loops['degree'][b] <= cutoff
                       for c in range(len(loops['degree']))
                       if loops['degree'][a]+loops['degree'][b]+loops['degree'][c] <= cutoff]
        self.index = {s: j for j, s in enumerate(self.states)}
        self.count = np.array([sum(x != 0 for x in s) for s in self.states])
        self.electric_C2 = np.array([sum(loops['C2'][x] for x in s) for s in self.states])
        self.d = matter['dimension']
        self.n = len(self.states)
        self.dimension = self.n*self.d
        self.field_B = [self.field_operator(b, k) for b in loops['B'] for k in range(3)]
        self.seed = np.zeros((self.n, self.d), complex)
        self.seed[0] = matter['seed']
        self.calls = 0
        inverse_local = np.argmax(loops['inversion'], axis=0)
        destinations = [self.index[(int(inverse_local[s[1]]), s[0], s[2])]
                        for s in self.states]
        self.rotation_inverse = np.argsort(destinations)
        self.spin_rotation = expm(-.5j*np.pi*matter['S'][2])
        self.c3_inverse = np.argsort([self.index[(s[2], s[0], s[1])] for s in self.states])
        self.c3_spin = expm(-2j*np.pi/3*sum(matter['S'])/np.sqrt(3))
        self.quarter_turns = []
        for axis in range(3):
            changed = []
            for s in self.states:
                if axis == 0:
                    target = (s[0], int(inverse_local[s[2]]), s[1])
                elif axis == 1:
                    target = (s[2], s[1], int(inverse_local[s[0]]))
                else:
                    target = (int(inverse_local[s[1]]), s[0], s[2])
                changed.append(self.index[target])
            self.quarter_turns.append((np.array(changed), np.argsort(changed),
                                       expm(-.5j*np.pi*matter['S'][axis])))

    def rotate(self, x: Any) -> Any:
        """C4z: axial x -> y, y -> -x; negative holonomy is inverted."""
        return (x@self.spin_rotation.T)[self.rotation_inverse]

    def polarized(self, x: Any) -> Any:
        """Exact C4 eigenphase exp(-i*pi/4), shared by both m=+1/2 seeds."""
        part = x
        result = x.copy()
        for n in range(1, 4):
            part = self.rotate(part)
            result += np.exp(.25j*np.pi*n)*part
        return result/4

    def rotate_c3(self, x: Any) -> Any:
        """Proper 120-degree rotation about (1,1,1): x -> y -> z -> x."""
        return (x@self.c3_spin.T)[self.c3_inverse]

    def cubic_representation(self, x: Any) -> dict[str, Any]:
        """Actual finite orbit and binary-cubic class sum, not an SO3 label."""
        c3x = self.rotate_c3(x)
        c3square = self.rotate_c3(c3x)
        orbit = np.array([x, c3x, c3square, self.rotate(c3x), self.rotate(c3square)])
        flat = orbit.reshape(len(orbit), -1)
        gram = flat.conj()@flat.T
        weights, vectors = eigh(gram)
        chosen = weights > 1e-9
        basis = (vectors[:, chosen].T@flat/np.sqrt(weights[chosen])[:, None])
        # Columns of vectors supply state coefficients, with no conjugation:
        # Gram_ij=<orbit_i|orbit_j>.
        basis = basis.reshape(sum(chosen), self.n, self.d)
        matrices, closures = {}, {}
        for name, action in (('C3_111', self.rotate_c3), ('C4z', self.rotate)):
            images = np.array([action(v) for v in basis])
            reduced = basis.reshape(len(basis), -1).conj()@images.reshape(len(basis), -1).T
            reconstructed = np.einsum('ij,ikl->jkl', reduced, basis)
            character = np.trace(reduced)
            matrices[name] = {'real': float(character.real), 'imag': float(character.imag)}
            closures[name] = float(np.linalg.norm(images-reconstructed))
        central = np.zeros_like(x)
        for dest, inverse, spin in self.quarter_turns:
            central += (x@spin.T)[inverse]+(x@spin.conj())[dest]
        target = 3*np.sqrt(2) if self.matter['I'] == .5 else 0.
        return {'orbit_dimension': len(basis), 'characters': matrices,
                'rotation_orbit_closure_norms': closures,
                'six_quarter_turn_class_sum_target': float(target),
                'class_sum_eigen_residual': float(np.linalg.norm(central-target*x)),
                'interpretation': 'Gamma6 doublet' if self.matter['I'] == .5 else 'Gamma8 quartet'}

    def field_operator(self, matrix: Any, axis: int) -> Any:
        rows, cols, vals = [], [], []
        transitions = [list(zip(np.flatnonzero(matrix[:, col]),
                                matrix[np.flatnonzero(matrix[:, col]), col]))
                       for col in range(len(self.loops['degree']))]
        for col, state in enumerate(self.states):
            for changed, value in transitions[state[axis]]:
                target = list(state); target[axis] = int(changed)
                row = self.index.get(tuple(target))
                if row is not None:
                    rows.append(row); cols.append(col); vals.append(value)
        return sparse.csr_matrix((vals, (rows, cols)), shape=(self.n, self.n))

    def interaction(self, vector: Any) -> Any:
        x = np.asarray(vector).reshape(self.n, self.d)
        out = np.zeros_like(x)
        for b, g in zip(self.field_B, self.matter['G']):
            out += (b@x)@g.T
        return out

    def operator(self, mass: float, electric: float) -> Any:
        def action(v: Any) -> Any:
            self.calls += 1
            x = self.polarized(np.asarray(v).reshape(self.n, self.d))
            return self.polarized(4*electric*self.electric_C2[:, None]*x-
                                  self.interaction(x)/(2*mass)).ravel()
        return LinearOperator((self.dimension, self.dimension), matvec=action, dtype=complex)

    def solve(self, mass: float, electric: float) -> dict[str, Any]:
        coefficient = -1/(2*mass)  # q cancels between Pauli coefficient and B.
        gap = 16*electric/3
        start = time.monotonic()
        self.calls = 0
        h = self.operator(mass, electric)
        values, vectors = eigsh(h, k=1, which='SA', v0=self.seed.ravel(),
                               tol=2e-10, maxiter=1400, ncv=24)
        energy = float(values[0])
        x = self.polarized(vectors[:, 0].reshape(self.n, self.d))
        x /= np.linalg.norm(x)
        residual = float(np.linalg.norm(h@x.ravel()-energy*x.ravel()))
        self.last_state, self.last_hamiltonian = x, h
        self.last_energy = energy
        weights = np.sum(abs(x)**2, axis=1)
        probabilities = [float(weights[self.count == n].sum()) for n in range(4)]
        axial = float(np.vdot(x, x@self.matter['A'].T).real)
        gauss = []
        for a in range(8):
            root = sum((self.field_operator(self.loops['Q'][a], k)@x for k in range(3)),
                       np.zeros_like(x)) + x@self.matter['Q'][a].T
            gauss.append(float(np.linalg.norm(root)))
        # Field-conditioned axial expectation, not an assumed probability law.
        conditional = []
        for n, probability in enumerate(probabilities):
            part = x[self.count == n]
            conditional.append(float(np.vdot(part, part@self.matter['A'].T).real/probability)
                               if probability > 1e-14 else None)
        source = self.interaction(self.seed)
        source /= np.linalg.norm(source)
        one_source_weight = float(abs(np.vdot(source, x))**2)
        vacuum_overlap = float(abs(np.vdot(self.seed[0], x[0]))**2)
        p = 1-probabilities[0]
        source_only_energy = (gap-np.sqrt(gap**2+4*coefficient**2*
                                         np.linalg.norm(self.interaction(self.seed))**2))/2
        source_only_p = .5*(1-gap/np.sqrt(gap**2+4*coefficient**2*
                                         np.linalg.norm(self.interaction(self.seed))**2))
        channel_rows = []
        diagonal_axial = 0.
        for casimir, spin2, basis in self.matter['channels']:
            part = x@basis.conj()
            reduced_a = basis.conj().T@self.matter['A']@basis
            contribution = float(np.vdot(part, part@reduced_a.T).real)
            diagonal_axial += contribution
            channel_rows.append({'colour_C2': casimir, 'quark_spin2': spin2,
                                 'dimension': basis.shape[1],
                                 'probability': float(np.vdot(part, part).real),
                                 'axial_diagonal_contribution': contribution})
        source_amplitude = np.vdot(source, x)
        vacuum_amplitude = np.vdot(self.seed, x)
        excluded = x-vacuum_amplitude*self.seed-source_amplitude*source
        included = x-excluded
        excluded_axial = float(np.vdot(excluded, excluded@self.matter['A'].T).real)
        axial_cross = float(2*np.vdot(included, excluded@self.matter['A'].T).real)
        raw_hx = 4*electric*self.electric_C2[:, None]*x+coefficient*self.interaction(x)
        result = {'I': self.matter['I'], 'mass_m': mass, 'electric_kE': electric,
                  'loop_representation': self.loops['name'],
                  'a': 1., 'q_for_B': .7, 'Pauli_coefficient_after_q_cancel': coefficient,
                  'one_loop_gap': gap, 'normalized_source_g': coefficient/np.sqrt(3),
                  'total_representation_degree_cutoff': self.cutoff, 'field_dimension': self.n,
                  'matter_dimension': self.d, 'tensor_dimension': self.dimension,
                  'energy': energy, 'eigen_residual': residual,
                  'root_Gauss_max_norm': max(gauss),
                  'unprojected_H_eigen_residual': float(np.linalg.norm(raw_hx-energy*x)),
                  'C4z_eigenphase_residual': float(np.linalg.norm(
                      self.rotate(x)-np.exp(-.25j*np.pi)*x)),
                  'cubic_representation': self.cubic_representation(x),
                  'probability_by_active_loop_count': probabilities,
                  'field_probability': p, 'vacuum_seed_overlap_Z0': vacuum_overlap,
                  'vacuum_weight_minus_seed_overlap': probabilities[0]-vacuum_overlap,
                  'one_source_projection_weight': one_source_weight,
                  'excluded_from_vacuum_plus_one_source_weight': 1-vacuum_overlap-one_source_weight,
                  'axial_neutral_tau3_sigma3': axial,
                  'axial_conditioned_on_active_loops': conditional,
                  'matter_colour_spin_channels': channel_rows,
                  'axial_cross_quark_spin_channels': axial-diagonal_axial,
                  'excluded_sector_axial_contribution': excluded_axial,
                  'included_excluded_axial_coherence': axial_cross,
                  'charge': float(np.vdot(x, x@self.matter['charge'].T).real),
                  'source_only_energy': float(source_only_energy),
                  'source_only_field_probability': float(source_only_p),
                  'H_matvec_calls': self.calls, 'elapsed_seconds': time.monotonic()-start}
        if self.matter['I'] == .5:
            result['old_probability_law_prediction_using_actual_p'] = 5/3-5*p/9
            result['probability_law_defect'] = axial-(5/3-5*p/9)
            result['axial_mobility_line_defect'] = axial-(10/9+5*vacuum_overlap/9)
            result['exact_three_component_axial_reconstruction'] = (
                5*vacuum_overlap/3+10*one_source_weight/9+excluded_axial+axial_cross)
            result['source_only_axial'] = 5/3-5*source_only_p/9
            result['conditional_empty_cell_three_hop_t0p2_interface_factor'] = -6*.2**3*vacuum_overlap
        return result


def inelastic_axial(n: Carrier, d: Carrier, row: dict[str, Any],
                    mass: float, electric: float) -> dict[str, Any]:
    """Full neutral A_z^3 response at fixed I3, including I=3/2 final states."""
    x, energy = n.last_state, n.last_energy
    a = n.matter['A']
    cross = d.matter['embedding'].conj().T@(
        n.matter['axial_full']@n.matter['embedding'])
    yn = x@a.T-row['axial_neutral_tau3_sigma3']*x
    yd = x@cross.T
    kn = (n.last_hamiltonian@yn.ravel()).reshape(yn.shape)-energy*yn
    kd = (d.last_hamiltonian@yd.ravel()).reshape(yd.shape)-energy*yd
    strength_n = float(np.vdot(yn, yn).real)
    strength_d = float(np.vdot(yd, yd).real)
    mu0 = strength_n+strength_d
    mu1 = float((np.vdot(yn, kn)+np.vdot(yd, kd)).real)
    mu2 = float((np.vdot(kn, kn)+np.vdot(kd, kd)).real)
    result: dict[str, Any] = {
        'operator': 'sum_quarks sigma_z tau3; full I3=+1/2 CAR response',
        'inelastic_strength_mu0': mu0, 'energy_weighted_mu1': mu1,
        'energy_squared_mu2': mu2, 'centroid_mu1_over_mu0': mu1/mu0,
        'strength_final_Ihalf': strength_n, 'strength_final_Ithreehalf': strength_d}
    p = (5/3-row['axial_neutral_tau3_sigma3'])*9/5
    if not 0 < p < .5:
        result['matched_25_state_reference'] = {
            'exists_positive_gap': False, 'required_field_probability': p,
            'reason': 'Positive-gap two-layer N law requires 0<p<1/2.'}
        return result
    omega = -energy*(1-2*p)/p
    coupling = np.sqrt(energy**2*(1-p)/(16*p))
    full_a = np.block([[n.matter['A'], cross.conj().T], [cross, d.matter['A']]])
    full_g = [sparse.block_diag((gn, gd)).toarray()
              for gn, gd in zip(n.matter['G'], d.matter['G'])]
    seed = np.concatenate((n.matter['seed'], np.zeros(d.d)))
    source = np.array([g@seed for g in full_g])/4
    reference = np.zeros((25, n.d+d.d), complex)
    reference[0] = np.sqrt(1-p)*seed
    reference[1:] = -np.sqrt(p)*source
    def href(v: Any) -> Any:
        out = np.zeros_like(v)
        out[0] = coupling*sum((g@v[j+1] for j, g in enumerate(full_g)), np.zeros(n.d+d.d))
        out[1:] = omega*v[1:]+coupling*np.array([g@v[0] for g in full_g])
        return out
    aref = float(np.vdot(reference, reference@full_a.T).real)
    yref = reference@full_a.T-aref*reference
    kref = href(yref)-energy*yref
    ref0 = float(np.vdot(yref, yref).real)
    ref1 = float(np.vdot(yref, kref).real)
    ref2 = float(np.vdot(kref, kref).real)
    result['matched_25_state_reference'] = {
        'exists_positive_gap': True, 'fit_field_probability': p,
        'fit_omega': omega, 'fit_g': coupling,
        'ground_energy': energy, 'ground_axial': aref,
        'ground_eigen_residual': float(np.linalg.norm(href(reference)-energy*reference)),
        'inelastic_strength_mu0': ref0, 'energy_weighted_mu1': ref1,
        'energy_squared_mu2': ref2, 'centroid_mu1_over_mu0': ref1/ref0,
        'actual_minus_reference_mu0': mu0-ref0,
        'actual_minus_reference_mu1': mu1-ref1,
        'actual_minus_reference_mu2': mu2-ref2}
    return result


def canonical_exterior_maps(loops: dict[str, Any], carrier: Carrier) -> dict[str, Any]:
    """Exact Q_degree3 B P_degree2 action, not an interior PW805 Hamiltonian."""
    representations = loops['_representations']
    by_name = {name: generators for name, generators, _, _ in representations}
    cube_basis = []
    for indices in itertools.combinations_with_replacement(range(3), 3):
        v = np.zeros((3, 3, 3), complex)
        permutations = set(itertools.permutations(indices))
        for index in permutations:
            v[index] = 1/np.sqrt(len(permutations))
        cube_basis.append(v.ravel())
    e10 = np.array(cube_basis).T
    ten = [e10.conj().T@(np.kron(np.kron(t, np.eye(3)), np.eye(3))+
            np.kron(np.kron(np.eye(3), t), np.eye(3))+
            np.kron(np.kron(np.eye(3), np.eye(3)), t))@e10 for t in COLOUR]
    parent15 = [np.kron(t, np.eye(8))+np.kron(np.eye(3), u)
                for t, u in zip(COLOUR, by_name['8'])]
    c2 = sum((t@t for t in parent15), np.zeros((24, 24), complex))
    values, vectors = eigh(c2)
    e15 = vectors[:, abs(values-16/3) < 1e-9]
    if e15.shape[1] != 15:
        raise RuntimeError('Incorrect (2,1) SU3 representation dimension')
    fifteen = [e15.conj().T@t@e15 for t in parent15]
    new = [('10', ten), ('bar10', [-t.conj() for t in ten]),
           ('15', fifteen), ('bar15', [-t.conj() for t in fifteen])]
    by_name.update(new)
    offsets = dict(loops['_offsets'])
    degree = list(loops['degree'])
    labels = ['retained']*155
    total = 155
    for name, generators in new:
        offsets[name] = total
        size = len(generators[0])**2
        total += size
        degree += [3]*size
        labels += [name]*size
    magnetic = [np.zeros((total, 155), complex) for _ in COLOUR]
    for b, previous in zip(magnetic, loops['B']):
        b[:155] = previous
    residual = 0.
    # Both U_out,in and U_in,out are required because b contains U-Udag.
    pairs = [('6', '10'), ('8', '15'), ('bar6', 'bar15'),
             ('bar10', 'bar6'), ('15', '6'), ('bar15', '8')]
    for r, s in pairs:
        tr, ts = by_name[r], by_name[s]
        dr, ds = len(tr[0]), len(ts[0])
        parents = [np.kron(t, np.eye(dr))+np.kron(np.eye(3), u)
                   for t, u in zip(COLOUR, tr)]
        equations = np.vstack([np.kron(np.eye(ds), a)-np.kron(b.T, np.eye(3*dr))
                               for a, b in zip(parents, ts)])
        kernel = null_space(equations, rcond=1e-11)
        if kernel.shape[1] != 1:
            raise RuntimeError(f'Nonunique exterior intertwiner {r}->{s}')
        cg = kernel[:, 0].reshape(3*dr, ds, order='F')*np.sqrt(ds)
        residual = max(residual, max(float(np.linalg.norm(a@cg-cg@b))
                                    for a, b in zip(parents, ts)))
        cg = cg.reshape(3, dr, ds)
        for a, b in itertools.product(range(3), repeat=2):
            block = np.sqrt(dr/ds)*np.einsum('mM,nN->MNmn', cg[a], cg[b].conj())
            block = block.reshape(ds*ds, dr*dr)
            for generator, target in zip(COLOUR, magnetic):
                if offsets[r] < 155:
                    target[offsets[s]:offsets[s]+ds*ds, offsets[r]:offsets[r]+dr*dr] += (
                        1j*generator[b, a]*block)
                else:
                    target[offsets[r]:offsets[r]+dr*dr, offsets[s]:offsets[s]+ds*ds] += (
                        -1j*generator[b, a].conjugate()*block.conj().T)
    for b in magnetic:
        b[abs(b) < 2e-14] = 0
    exterior = [(a, b, c) for a in range(total) for b in range(total)
                if degree[a]+degree[b] <= 3 for c in range(total)
                if degree[a]+degree[b]+degree[c] == 3]
    index = {s: j for j, s in enumerate(exterior)}
    maps = []
    for b in magnetic:
        transitions = [list(zip(np.flatnonzero(b[:, col]), b[np.flatnonzero(b[:, col]), col]))
                       for col in range(155)]
        for axis in range(3):
            rows, cols, vals = [], [], []
            for col, state in enumerate(carrier.states):
                for target, value in transitions[state[axis]]:
                    changed = list(state); changed[axis] = int(target)
                    row = index.get(tuple(changed))
                    if row is not None:
                        rows.append(row); cols.append(col); vals.append(value)
            maps.append(sparse.csr_matrix((vals, (rows, cols)),
                                         shape=(len(exterior), carrier.n)))
    categories = []
    for state in exterior:
        new_labels = [labels[x] for x in state if degree[x] == 3]
        if new_labels:
            categories.append('single_loop_10_bar10' if '10' in new_labels[0]
                              else 'single_loop_15_bar15')
        elif any(degree[x] == 2 for x in state):
            categories.append('two_loops_degree2_plus1')
        else:
            categories.append('three_loops_fundamental')
    return {'maps': maps, 'categories': np.array(categories),
            'exterior_field_dimension': len(exterior),
            'local_representation_dimension_for_exterior': total,
            'CG_intertwining_residual': residual}


def exterior_axial_response(exterior: dict[str, Any], n: Carrier, d: Carrier,
                           row: dict[str, Any], mass: float) -> dict[str, Any]:
    """Positive exterior correction on the retained trial state, without rediagonalizing."""
    def action(x: Any, matter: dict[str, Any]) -> Any:
        out = np.zeros((exterior['exterior_field_dimension'], x.shape[1]), complex)
        for b, g in zip(exterior['maps'], matter['G']):
            out -= ((b@x)@g.T)/(2*mass)
        return out
    x = n.last_state
    ga = row['axial_neutral_tau3_sigma3']
    cross = d.matter['embedding'].conj().T@(
        n.matter['axial_full']@n.matter['embedding'])
    yn = x@n.matter['A'].T-ga*x
    yd = x@cross.T
    en, ed = action(yn, n.matter), action(yd, d.matter)
    leakage = action(x, n.matter)
    component = np.sum(abs(en)**2, axis=1)+np.sum(abs(ed)**2, axis=1)
    correction = float(component.sum())
    # For the nonstationary full-H trial state, commutator norm and the
    # centered response-vector norm differ by the ground-state exterior.
    cn = en-(leakage@n.matter['A'].T-ga*leakage)
    cd = ed-leakage@cross.T
    degree = np.array([sum(n.loops['degree'][v] for v in state) for state in n.states])
    weights = np.sum(abs(x)**2, axis=1)
    mu2 = row['inelastic_axial']['energy_squared_mu2']
    return {'meaning': 'Exact canonical H action on retained P2 trial state; not exact full-H ground response.',
            'internal_P2_response_mu2': mu2,
            'positive_exterior_centered_response_norm2': correction,
            'full_H_centered_trial_response_mu2': mu2+correction,
            'exterior_response_by_field_channel': {
                str(label): float(component[exterior['categories'] == label].sum())
                for label in np.unique(exterior['categories'])},
            'ground_state_H_exterior_norm2': float(np.vdot(leakage, leakage).real),
            'full_H_commutator_norm2_on_retained_trial': float(
                mu2+np.vdot(cn, cn).real+np.vdot(cd, cd).real),
            'ground_probability_by_total_representation_degree':
                [float(weights[degree == j].sum()) for j in range(3)],
            'exterior_field_dimension': exterior['exterior_field_dimension'],
            'local_representation_dimension_for_exterior':
                exterior['local_representation_dimension_for_exterior'],
            'CG_intertwining_residual': exterior['CG_intertwining_residual'],
            'weak_coupling_power_counting':
                'P2 top-degree amplitude O(c^2), Q3 H deltaA psi O(c^3), positive norm O(c^6).'}


def run() -> dict[str, Any]:
    rows = []
    with threadpool_limits(limits=1):
        laws = [loop_law(), canonical_degree_two_loop()]
        for loops in laws:
            matter_n, matter_d = matter_law(.5), matter_law(1.5)
            for cutoff in ((1, 2, 3) if loops['name'] == 'PW19' else (2,)):
                n, d = Carrier(loops, matter_n, cutoff), Carrier(loops, matter_d, cutoff)
                exterior = canonical_exterior_maps(loops, n) if loops['name'] == 'PW155_degree2' else None
                for mass in (3., 1., .3):
                    nr, dr = n.solve(mass, .27), d.solve(mass, .27)
                    nr['inelastic_axial'] = inelastic_axial(n, d, nr, mass, .27)
                    if exterior is not None:
                        nr['canonical_degree3_exterior'] = exterior_axial_response(
                            exterior, n, d, nr, mass)
                    rows.extend((nr, dr))
                    print(json.dumps(nr), flush=True)
                    print(json.dumps(dr), flush=True)
    comparison = []
    for loops in laws:
        for cutoff in ((1, 2, 3) if loops['name'] == 'PW19' else (2,)):
            for mass in (3., 1., .3):
                selected = [x for x in rows if x['mass_m'] == mass and
                            x['loop_representation'] == loops['name'] and
                            x['total_representation_degree_cutoff'] == cutoff]
                n = next(x for x in selected if x['I'] == .5)
                d = next(x for x in selected if x['I'] == 1.5)
                comparison.append({'representation': loops['name'], 'cutoff': cutoff, 'mass': mass,
                                   'Delta_minus_nucleon': d['energy']-n['energy'],
                                   'source_only_Delta_minus_nucleon':
                                       d['source_only_energy']-n['source_only_energy']})
    sources = [Path(__file__), HERE/'rich_nucleon_hierarchy_20261019.py',
               HERE/'rich_composite_local_flux_binding_20261017_model.py',
               HERE/'rich_colour_spin_root_20261019.py']
    return {'kind': 'CANONICAL_TWO_SOURCE_AND_REPEATED_PW19_CUBIC_FLUX_BARYON',
            'plan': 'issue222-journal-construction-20261019',
            'runtime': {'python': platform.python_version(), 'numpy': np.__version__,
                        'scipy': scipy.__version__},
            'source_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
            'local_operator_evidence': {l['name']: {k: v for k, v in l.items()
                if k not in ('B', 'Q', 'inversion', 'degree', 'C2',
                             '_representations', '_offsets')} for l in laws},
            'rows': rows, 'ordering': comparison,
            'scope': ['PW19 full three-loop tensor and PW155 total source-degree2 are distinct calculations.',
                      'PW155 degree2 contains every exact B squared vacuum source irrep, not full L2 closure.',
                      'Independent four-edge loops sharing only root: defined electric bouquet law.',
                      'This is not an exact maximal-tree electric reduction of the previous nine-link graph.',
                      'Signed permutations plus loop inversion realize cubic axial symmetry, not continuous SO3.',
                      'Primitive CAR, flavour, mass, lattice spacing, kE, and Wilson-Pauli operator remain inputs.',
                      'm=.3 and m=1 are stress regimes, not controlled heavy-Dirac continuum predictions.',
                      'No colour magnetic plaquette potential, spatial bound state, or sea matching is asserted.']}


if __name__ == '__main__':
    result = run()
    output = HERE/'rich_journal_flux_carrier_20261019_result.json'
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'output': str(output), 'ordering': result['ordering']}), flush=True)
