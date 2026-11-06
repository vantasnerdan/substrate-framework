"""Order g^2 t^3 physical SU(3) baryon hopping, not a three-hop moment.

Two independent 0+adjoint-vector cells and a Haar SU(3) tree link.
The one-field restriction is exact for this coefficient, not a full Hilbert
space truncation claim. Units hbar=c=1; no scalar cloud is appended.
"""
from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse

from rich_nucleon_hierarchy_20261019 import (
    COLOUR, SIGMA, FIELD_COLOUR, matter_sector, lift,
)
from rich_composite_local_flux_binding_20261017_model import fermion, link_operators

HERE = Path(__file__).resolve().parent


def bits_sector(left_number, up_number):
    return [sum(1 << k for k in occupied)
            for occupied in itertools.combinations(range(24), 3)
            if sum(k < 12 for k in occupied) == left_number
            and sum((k % 12) < 6 for k in occupied) == up_number]


class Orbit:
    """CAR x link orbit; field channels are vacuum, 24 left, 24 right."""

    def __init__(self, up_number, angle=0.):
        self.bits = [bits_sector(n, up_number) for n in (3, 2, 1, 0)]
        self.index = [{b: j for j, b in enumerate(s)} for s in self.bits]
        self.transport, self.left, self.right = link_operators()
        self.hops = []
        for stage in range(3):
            rows, cols, vals, charges = [], [], [], []
            for col, bits in enumerate(self.bits[stage]):
                for flavour, spin, a, b in itertools.product(range(2), range(2), range(3), range(3)):
                    change = fermion(bits, [(True, 12+6*flavour+3*spin+a),
                                            (False, 6*flavour+3*spin+b)])
                    if change is None:
                        continue
                    newbits, sign = change
                    q = 2/3 if flavour == 0 else -1/3
                    u = self.transport[a, b, False]
                    for out, old in zip(*np.nonzero(u)):
                        rows.append(int(out)*len(self.bits[stage+1])+self.index[stage+1][newbits])
                        cols.append(int(old)*len(self.bits[stage])+col)
                        vals.append(-sign*u[out, old]*np.exp(-1j*q*angle))
                        charges.append(q)
            shape = (19*len(self.bits[stage+1]), 19*len(self.bits[stage]))
            self.hops.append([sparse.csr_matrix((np.array(vals)*(-1j*np.array(charges))**d,
                                                (rows, cols)), shape=shape) for d in range(3)])
        self.mixed = []
        self.colour = []
        for stage in range(4):
            operators, colour_ops = [], []
            for site in range(2):
                gs, cs = [], []
                for t in COLOUR:
                    mat = np.zeros((24, 24), complex)
                    mat[12*site:12*(site+1), 12*site:12*(site+1)] = np.kron(np.eye(4), t)
                    cs.append(lift(mat, self.bits[stage], self.index[stage]))
                    for sigma in SIGMA:
                        mat = np.zeros((24, 24), complex)
                        mat[12*site:12*(site+1), 12*site:12*(site+1)] = np.kron(np.kron(np.eye(2), sigma), t)
                        gs.append(lift(mat, self.bits[stage], self.index[stage]))
                operators.append(gs); colour_ops.append(cs)
            self.mixed.append(operators); self.colour.append(colour_ops)

    def zero(self, stage):
        return np.zeros((49, 19, len(self.bits[stage])), complex)

    def hop(self, wave, stage, derivative=0):
        flat = wave.reshape(49, -1)
        return (self.hops[stage][derivative] @ flat.T).T.reshape(49, 19, -1)

    def vertex(self, wave, stage):
        out = self.zero(stage)
        for site in range(2):
            for a, op in enumerate(self.mixed[stage][site]):
                f = 1+24*site+a
                out[f] = (op @ wave[0].T).T
                out[0] += (op @ wave[f].T).T
        return out

    def inverse(self, wave, delta, omega, power=1):
        gap = np.full((49, 1, 1), delta+omega)
        gap[0] = delta
        return wave/gap**power

    def resolvent(self, wave, stage, order, delta, omega, lam, part=None):
        # R=(E_N-H)^-1; E_N=-lam*g^2/omega+O(g^4).
        inv = lambda x, power=1: self.inverse(x, delta, omega, power)
        if order == 0:
            return -inv(wave)
        if order == 1:
            return inv(self.vertex(inv(wave), stage))
        shift = lam/omega*inv(wave, 2)
        vv = -inv(self.vertex(inv(self.vertex(inv(wave), stage)), stage))
        return shift if part == 'energy' else vv if part == 'VV' else shift+vv

    def gauss_norm(self, wave, stage):
        worst = 0.
        for site, bond in [(0, self.right), (1, self.left)]:
            for a in range(8):
                out = np.empty_like(wave)
                op = self.colour[stage][site][a]
                for f in range(49):
                    out[f] = (op @ wave[f].T).T + bond[a] @ wave[f]
                base = 1+24*site
                out[base:base+24] += np.einsum('ab,blm->alm', FIELD_COLOUR[a], wave[base:base+24])
                worst = max(worst, float(np.linalg.norm(out)))
        return worst


def endpoint(orbit, psi, local_bits, stage, omega, lam):
    zero = orbit.zero(stage)
    for bits, value in zip(local_bits, psi):
        global_bits = bits if stage == 0 else bits << 12
        if global_bits in orbit.index[stage]:
            zero[0, 0, orbit.index[stage][global_bits]] = value
    return [zero, -orbit.vertex(zero, stage)/omega, -lam/(2*omega**2)*zero]


def matrix_element(orbit, initial, final, r1, r2, delta, omega, lam,
                   part1=None, part2=None, jets=False):
    waves = [initial, np.zeros_like(initial), np.zeros_like(initial)] if jets else [initial]
    gauss = []
    for stage in range(3):
        new = []
        for n in range(len(waves)):
            new.append(sum((math.comb(n, k)*orbit.hop(waves[n-k], stage, k)
                            for k in range(n+1)), orbit.zero(stage+1)))
        waves = new
        if stage < 2:
            order, part = (r1, part1) if stage == 0 else (r2, part2)
            waves = [orbit.resolvent(w, stage+1, order, delta, omega, lam, part) for w in waves]
        if jets:
            gauss.append(orbit.gauss_norm(waves[0], stage+1))
    return [complex(np.vdot(final, w)) for w in waves], gauss


def propagation(psi, local_bits, up_number, delta, omega, lam, angle=0., jets=False):
    orbit = Orbit(up_number, angle)
    initial = endpoint(orbit, psi, local_bits, 0, omega, lam)
    final = endpoint(orbit, psi, local_bits, 3, omega, lam)
    terms = {}
    total = np.zeros(3 if jets else 1, complex)
    gauss_max = 0.
    for i, j, r1, r2 in itertools.product(range(3), repeat=4):
        if i+j+r1+r2 != 2:
            continue
        label = f'endpoint_left_{i}_right_{j}_R1_{r1}_R2_{r2}'
        values, norms = matrix_element(orbit, initial[i], final[j], r1, r2,
                                      delta, omega, lam, jets=jets)
        total += values
        terms[label] = float(values[0].real)
        gauss_max = max([gauss_max]+norms)
    bare, norms = matrix_element(orbit, initial[0], final[0], 0, 0, delta, omega, lam, jets=jets)
    split = {}
    for r in (0, 1):
        for part in ('energy', 'VV'):
            val, _ = matrix_element(orbit, initial[0], final[0], 2 if r == 0 else 0,
                                    2 if r == 1 else 0, delta, omega, lam,
                                    part1=part if r == 0 else None,
                                    part2=part if r == 1 else None)
            split[f'R{r+1}_{part}'] = float(val[0].real)
    # Actual backtracking return, including both endpoints and the same R.
    hopped = [[orbit.hop(endpoint_wave, 0, d) for d in range(3)]
              for endpoint_wave in initial]
    diagonal_bare = np.zeros(3, complex)
    diagonal_g2 = np.zeros(3, complex)
    for i, j, order in itertools.product(range(3), repeat=3):
        if i+j+order not in (0, 2):
            continue
        destination = diagonal_bare if i+j+order == 0 else diagonal_g2
        for derivative in range(3):
            destination[derivative] += sum(
                math.comb(derivative, k)*np.vdot(
                    hopped[j][k], orbit.resolvent(hopped[i][derivative-k], 1,
                                                order, delta, omega, lam))
                for k in range(derivative+1))
    return {'bare': np.array(bare), 'g2': total, 'terms': terms, 'split': split,
            'gauss_max': max([gauss_max]+norms),
            'diagonal_bare': diagonal_bare, 'diagonal_g2': diagonal_g2,
            'dimensions': [49*19*len(s) for s in orbit.bits]}


def same_state_observables(data, psi, g, omega, mass):
    lam = float(np.vdot(psi, data['g2'] @ psi).real)
    images = np.array([op @ psi for op in data['mixed']])
    energy = (omega-np.sqrt(omega**2+4*g*g*lam))/2
    alpha = 1/np.sqrt(1+g*g*lam/(energy-omega)**2)
    one = alpha*g/(energy-omega)*images
    states = [sum(1 << k for k in occupied) for occupied in itertools.combinations(range(12), 3)]
    axial, magnetic = [], []
    for bits in states:
        az, mz = 0., 0.
        for k in range(12):
            if bits >> k & 1:
                spin = 1 if (k//3) % 2 == 0 else -1
                az += spin*(1 if k//6 == 0 else -1)
                mz += spin*(2/3 if k//6 == 0 else -1/3)/(2*mass)
        axial.append(az); magnetic.append(mz)
    def expectation(diagonal):
        diagonal = np.array(diagonal)
        return float((alpha**2*np.vdot(psi, diagonal*psi)
                      +np.vdot(one, one*diagonal[None, :])).real)
    return {'lambda': lam, 'energy': float(energy), 'rho': float(np.vdot(one, one).real),
            'gA': expectation(axial), 'mu_quark_Pauli': expectation(magnetic)}


def threshold_data(g, omega, delta, nucleon_energy):
    spectra = []
    for n in (1, 2):
        states = [sum(1 << k for k in occupied) for occupied in itertools.combinations(range(12), n)]
        index = {b: j for j, b in enumerate(states)}
        mixed = [lift(np.kron(np.kron(np.eye(2), s), t), states, index) for t in COLOUR for s in SIGMA]
        g2 = sum((op @ op for op in mixed), sparse.csr_matrix((len(states),)*2, dtype=complex))
        values = np.linalg.eigvalsh(g2.toarray())
        lam = float(values[-1])
        energy = (omega-np.sqrt(omega**2+4*g*g*lam))/2
        spectra.append({'N': n, 'G2_eigenvalues': sorted(set(np.round(values, 10))),
                        'lambda_max': lam, 'unrestricted_local_min_energy': float(energy)})
    lower = delta+sum(x['unrestricted_local_min_energy'] for x in spectra)
    return {'local_spectra': spectra, 'electric_cost': delta,
            'physical_intermediate_threshold_lower_bound': float(lower),
            'compact_carrier_energy': nucleon_energy,
            'threshold_minus_carrier_lower_bound': float(lower-nucleon_energy),
            'meaning': 'Gauss restriction can only raise this unrestricted bound; positive margin proves ordering.'}


def pair(z):
    return [float(z.real), float(z.imag)]


def run():
    mass, spacing, k_e = 3., 1., .27
    delta = 4*k_e/3
    omega = 4*delta
    g = -1/(2*np.sqrt(3)*mass*spacing**2)
    t = 1/(2*mass*spacing**2)
    data, levels = matter_sector(3)
    proton = next(level['psi'] for level in levels if level['S'] == .5)
    delta_psi = next(level['psi'] for level in levels if level['S'] == 1.5)
    neutron = (data['isospin'][0]-1j*data['isospin'][1]) @ proton
    neutron /= np.linalg.norm(neutron)
    states = [sum(1 << k for k in occupied) for occupied in itertools.combinations(range(12), 3)]
    lam = float(np.vdot(proton, data['g2'] @ proton).real)
    p = same_state_observables(data, proton, g, omega, mass)
    n = same_state_observables(data, neutron, g, omega, mass)
    d = same_state_observables(data, delta_psi, g, omega, mass)
    propagation_zero = propagation(proton, states, 2, delta, omega, lam)
    lam_delta = float(np.vdot(delta_psi, data['g2'] @ delta_psi).real)
    propagation_delta = propagation(delta_psi, states, 2, delta, omega, lam_delta)
    c_delta = float((propagation_delta['g2'][0]/propagation_delta['bare'][0]).real)
    c = float((propagation_zero['g2'][0]/propagation_zero['bare'][0]).real)
    amplitude = t**3*(propagation_zero['bare'][0]+g*g*propagation_zero['g2'][0])
    gauged = []
    for name, psi, ups, q, source_lam in [
            ('proton', proton, 2, 1., lam), ('neutron', neutron, 1, 0., lam),
            ('Delta_plus', delta_psi, 2, 1., lam_delta)]:
        res = propagation(psi, states, ups, delta, omega, source_lam, angle=.37, jets=True)
        jets = t**3*(res['bare']+g*g*res['g2'])
        gauged.append({'source': name, 'Q': q, 'angle': .37, 'transfer_current_contact': [pair(z) for z in jets],
                       'current_covariance_residual': float(abs(jets[1]+1j*q*jets[0])),
                       'contact_covariance_residual': float(abs(jets[2]+q*q*jets[0])),
                       'intermediate_Gauss_max_norm_for_coefficient_terms': res['gauss_max']})
    rho2 = lam*g*g/omega**2
    ga2 = 5/3-5*rho2/9
    ds2 = 8*g*g/omega
    mu_p2 = 1/(2*mass)*(1-rho2/3)
    mu_n2 = -1/(3*mass)*(1-rho2/3)
    curve = []
    for curve_mass in (3., 6., 12.):
        curve_t = 1/(2*curve_mass*spacing**2)
        curve_g = -curve_t/np.sqrt(3)
        nucleon_obs = same_state_observables(data, proton, curve_g, omega, curve_mass)
        delta_obs = same_state_observables(data, delta_psi, curve_g, omega, curve_mass)
        bands = {}
        for name, source_lam, source in [('N', lam, propagation_zero),
                                         ('Delta', lam_delta, propagation_delta)]:
            hop = float((curve_t**3*(source['bare'][0]+curve_g**2*source['g2'][0])).real)
            sigma = float((curve_t**2*(source['diagonal_bare'][0]+curve_g**2*source['diagonal_g2'][0])).real)
            rest = -source_lam*curve_g**2/omega
            bands[name] = {'E_local_g2': rest, 'return_per_bond_g2t2': sigma,
                           'transfer_g2t3': hop, 'band_center': rest+2*sigma,
                           'bandwidth': 4*abs(hop), 'curvature': -2*hop*spacing**2,
                           'E_k0': rest+2*sigma+2*hop, 'E_kPi': rest+2*sigma-2*hop}
            exact_rest = nucleon_obs['energy'] if name == 'N' else delta_obs['energy']
            bands[name]['exact_25cell_rest_plus_truncated_kinetics'] = {
                'E_local_exact': exact_rest, 'band_center': exact_rest+2*sigma,
                'E_k0': exact_rest+2*sigma+2*hop, 'E_kPi': exact_rest+2*sigma-2*hop}
        center_gap = bands['Delta']['band_center']-bands['N']['band_center']
        cosine_gap = 2*(bands['Delta']['transfer_g2t3']-bands['N']['transfer_g2t3'])
        exact_rest_center_gap = (delta_obs['energy']-nucleon_obs['energy']
                                 +2*(bands['Delta']['return_per_bond_g2t2']
                                     -bands['N']['return_per_bond_g2t2']))
        curve.append({'m': curve_mass, 'a': spacing, 'kE': k_e, 'g': curve_g,
                      't_NR': curve_t, 't_over_Delta': curve_t/delta,
                      'recommended_weak_hopping_example': curve_mass == 6.,
                      'gA_g2': 5/3-5*lam*curve_g**2/(9*omega**2),
                      'mu_p_g2': (1-lam*curve_g**2/(3*omega**2))/(2*curve_mass),
                      'local_Delta_split_g2': (lam-lam_delta)*curve_g**2/omega,
                      'local_Delta_split_exact': delta_obs['energy']-nucleon_obs['energy'],
                      'bands': bands,
                      'curvature_Delta_minus_N': bands['Delta']['curvature']-bands['N']['curvature'],
                      'curvature_difference_over_bare_curvature': (
                          (bands['Delta']['curvature']-bands['N']['curvature'])
                          /(12*curve_t**3*spacing**2/delta**2)),
                      'Delta_N_gap_k': {'constant': center_gap, 'cos_ka_coefficient': cosine_gap,
                                        'minimum': center_gap-abs(cosine_gap),
                                        'maximum': center_gap+abs(cosine_gap)},
                      'Delta_N_gap_exact25rest_truncatedkinetics': {
                          'constant': exact_rest_center_gap, 'cos_ka_coefficient': cosine_gap,
                          'minimum': exact_rest_center_gap-abs(cosine_gap),
                          'maximum': exact_rest_center_gap+abs(cosine_gap),
                          'meaning': 'Exact finite-25-cell rest only; kinetic terms remain g^2t^2,g^2t^3. Not a complete canonical g^4 expansion or rigorous finite-t positivity bound.'},
                      'threshold_ordering': threshold_data(curve_g, omega, delta, nucleon_obs['energy'])})
    return {'title': 'Physical-resolvent weak-flux baryon transfer and same-source response',
            'plan': {'id': 'issue222-journal-construction-20261019',
                     'pinned_revision': '3b74f3515e8112c13e9a553e241ec6c2cec7c1b466d59d0566f7e88bbb30ae1b',
                     'work_item': 'journal-result'},
            'root': {'m': mass, 'a': spacing, 'kE': k_e, 'Delta': delta, 'Omega': omega, 'g': g,
                     't_NR': t, 't_over_Delta': t/delta,
                     'matching': 'Pauli/Foldy-Wouthuysen kinetic -D^2/(2m): nearest-neighbour t=1/(2ma^2); local magnetic g=-1/(2sqrt(3)ma^2).'},
            'propagation': {'order': 'g^2 t^3, independent local field cells; no scalar cloud',
                            'bare_over_t3': pair(propagation_zero['bare'][0]),
                            'g2_coefficient_over_t3': pair(propagation_zero['g2'][0]),
                            'relative_g2_coefficient_C': c,
                            'dimensionless_C_Omega2': c*omega**2,
                            'relative_correction': g*g*c,
                            'analytic_coefficient': 'C Omega^2=-4(r+2)(4r^2+9r+6)/(3(r+1)^2), r=Omega/Delta',
                            'operator_contractions': {
                                'same_resolvent_VV': float(propagation_zero['split']['R1_VV']/propagation_zero['bare'][0].real*delta*(delta+omega)),
                                'cross_resolvent_VV': float(propagation_zero['terms']['endpoint_left_0_right_0_R1_1_R2_1']/propagation_zero['bare'][0].real*(delta+omega)**2),
                                'adjacent_endpoint_V': float(propagation_zero['terms']['endpoint_left_1_right_0_R1_1_R2_0']/propagation_zero['bare'][0].real*omega*(delta+omega)),
                                'far_endpoint_V': float(propagation_zero['terms']['endpoint_left_1_right_0_R1_0_R2_1']/propagation_zero['bare'][0].real*omega*(delta+omega)**2/delta)},
                            'transfer_at_A0': pair(amplitude),
                            'one_dimensional_band': 'E_B(k)=-lambda_B g^2/Omega+2 Sigma_B_per_bond+2 M_B cos(ka)+O(t^4,g^4 t^2,g^4); energies omit common quark rest mass and discretization constant.',
                            'diagonal_return_per_bond': [pair(z) for z in t*t*(propagation_zero['diagonal_bare']+g*g*propagation_zero['diagonal_g2'])],
                            'diagonal_return_meaning': 'Self energy, its first and second A derivatives; 1D has two adjacent empty bonds.',
                            'band_curvature_at_zero': float(-2*amplitude.real*spacing**2),
                            'bandwidth': float(4*abs(amplitude)),
                            'insertion_terms_over_t3': propagation_zero['terms'],
                            'resolvent_second_order_split': propagation_zero['split'],
                            'orbit_dimensions': propagation_zero['dimensions']},
            'Delta_plus_propagation': {'source_S': 1.5, 'source_I': 1.5, 'source_Jz': .5,
                                       'source_Iz': .5,
                                       'source_Q': float(np.vdot(delta_psi, data['charge']*delta_psi).real),
                                       'lambda': lam_delta, 'relative_g2_coefficient_C': c_delta,
                                       'dimensionless_C_Omega2': c_delta*omega**2,
                                       'g2_coefficient_over_t3': pair(propagation_delta['g2'][0]),
                                       'insertion_terms_over_t3': propagation_delta['terms'],
                                       'resolvent_second_order_split': propagation_delta['split'],
                                       'not_half_N_coefficient': c_delta-c/2},
            'matched_mass_curve': curve,
            'same_source': {'proton': p, 'neutron': n, 'Delta': d,
                            'exact_local_Delta_split': d['energy']-p['energy'],
                            'order_g2': {'rho': rho2, 'gA': ga2, 'Delta_split': ds2,
                                         'mu_p': mu_p2, 'mu_n': mu_n2}},
            'relation': {'transfer_relative_shift_per_axial_deficit': c*omega**2/(lam*5/9),
                         'transfer_relative_shift_per_Delta_split_over_Omega': c*omega**2/8,
                         'primary_curvature_relation': '(kappa_D-kappa_N)/kappa_bare=(1096/25)(g/Omega)^2=(1233/250)(5/3-gA)=(137/25)delta_local/Omega through g^2t^3',
                         'momentum_gap_relation': 'gap_g2(k)=delta_local_g2*[1-4(t/Delta)^2-(411/25)(t/Delta)^3*cos(ka)] at Omega=4Delta; exact25rest alternative recorded separately',
                         'formula': '(M/M_bare-1)=C g^2=(9 C Omega^2/80)(5/3-gA)=(C Omega^2/8)(Delta_N/Omega)+O(g^4); mu_p=(3 gA/5)/(2m), mu_n=-(3 gA/5)/(3m) at g^2.',
                         'novelty_scope': 'New calculable transfer/axial/splitting correlation in this specified quantum-cell truncation, not a claim that axial recoupling or strong-coupling baryon hopping is new.'},
            'threshold_ordering': threshold_data(g, omega, delta, p['energy']),
            'gauged_resolvents': gauged,
            'limitations': ['Finite 25-state local field law, not a converged canonical loop tower.',
                            't/Delta=0.463: leading t^3 band coefficient is calculated; quantitative t-series convergence is not established.',
                            'Delta larger curvature is a nonrelativistic lattice-model/in-medium effect, not a Lorentz-invariant vacuum QCD mass-dispersion prediction.',
                            'Primary result is the curvature difference, which cancels spin-flavor-universal bare higher hopping corrections; no complete matched single-parameter heavy-mass series is claimed.',
                            'No fully eliminated sea or retarded sea has been spliced into the explicit tree flux.',
                            'Exact local observables are displayed for comparison but correlation is consistently truncated at g^2.',
                            'Tree link 6 and higher representations cannot reach final singlet in three forward hops: after two hops only 3bar can close with the third 3; local vertices preserve link irrep.',
                            'At g^2 two-field channels cannot close onto vacuum endpoints; two independent one-field endpoint sectors and all allowed insertions are retained.']}


if __name__ == '__main__':
    result = run()
    target = HERE / 'rich_journal_baryon_band_20261019_result.json'
    target.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'propagation': result['propagation'], 'relation': result['relation'],
                      'threshold_ordering': result['threshold_ordering'],
                      'gauged_resolvents': result['gauged_resolvents']}, indent=2))
