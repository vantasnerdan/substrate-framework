#!/usr/bin/env python3
"""Consume admitted N/material receipts; run only a new formal identifiability oracle.

No nuclear/material solver is imported or replayed. No physical rate or emitted
spectrum is inferred. See nuclear-medium-assessment-registration.json, written
before this implementation, for source identities, definitions and mutations.
Requires SymPy. Root runs this under the assigned CPU/thread restriction.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
REGISTRATION = HERE / 'nuclear-medium-assessment-registration.json'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def admit(binding, *, primary_dir=None):
    path = Path(binding['path'])
    if path.is_absolute():
        if primary_dir is not None:
            path = primary_dir / path.name
    else:
        path = ROOT / path
    require(path.stat().st_size == int(binding['bytes']), f'Byte count mismatch: {path}')
    require(digest(path) == binding['sha256'], f'SHA256 mismatch: {path}')
    return path


def json_bytes(path):
    return json.loads(path.read_bytes())


def require_receipts(receipts):
    required = {
        'nuclear-product-energy-ledger.json', 'measured-nuclear-channels-001.json',
        'nuclear-energy-transport-001.json', 'full-common-reference-004/summary.json',
        'full-common-reference-004/full-common-reference.json',
        'finiteq-material-001/result.json', 'full-field-aggregate-001/result.json',
        'ionic-gaussian-material-002/result.json', 'nuclear-oracle-command-002.json',
        'measured-nuclear-oracle-001.json', 'nuclear-transport-oracle-002.json',
    }
    require(set(receipts) == required, 'Missing or unregistered published receipt')


def formal_oracle():
    """Independent full-space evolution derivative versus reduced contraction.

    This finite Hamiltonian is a counterexample to identification from the
    medium alone, not a proposed material Hamiltonian or a golden-rule bath.
    """
    import sympy as s

    hbar = s.Symbol('hbar', positive=True)
    coupling = 2 + s.I
    bath = s.Matrix([[1, 2 + s.I], [3*s.I, 2]])
    rho_b = s.diag(s.Rational(1, 3), s.Rational(2, 3))
    nuclear_transition = s.Matrix([[0, 0], [1, 0]])
    entrance = s.diag(1, 0)
    exit_projector = s.kronecker_product(s.diag(0, 1), s.eye(2))
    rho = s.kronecker_product(entrance, rho_b)

    def full_coefficient(g):
        transition = s.kronecker_product(g*nuclear_transition, bath)
        hamiltonian = transition + transition.H
        commutator = hamiltonian*rho-rho*hamiltonian
        second = -(hamiltonian*commutator-commutator*hamiltonian)/(2*hbar**2)
        return s.simplify(s.trace(exit_projector*second))

    medium_norm = s.simplify(s.trace(bath*rho_b*bath.H))
    reduced = s.simplify(s.conjugate(coupling)*coupling*medium_norm/hbar**2)
    direct = full_coefficient(coupling)
    doubled = full_coefficient(2*coupling)
    require(s.simplify(direct-reduced) == 0, 'Full-space/reduced contraction mismatch')
    require(s.simplify(doubled-4*direct) == 0 and direct != 0,
            'Nuclear-coupling identifiability witness failed')

    # Column b is the initial medium state; row bp is the final state.
    energies = [s.Integer(0), s.Integer(5)]
    time = s.Symbol('t', real=True)
    propagator = s.diag(*(s.exp(-s.I*energy*time/hbar) for energy in energies))
    direct_correlation = s.trace(rho_b*propagator.H*bath.H*propagator*bath)
    spectral_correlation = s.Integer(0)
    delta = s.Symbol('Delta', real=True)
    transitions = []
    spectral = s.Integer(0)
    reversed_spectral = s.Integer(0)
    weight_sum = s.Integer(0)
    for b in range(2):
        for bp in range(2):
            weight = s.simplify(rho_b[b, b]*s.conjugate(bath[bp, b])*bath[bp, b])
            energy_transfer = energies[bp]-energies[b]
            transitions.append({'initial': b, 'final': bp,
                                'medium_energy_gain': str(energy_transfer),
                                'weight': str(weight)})
            spectral_correlation += weight*s.exp(-s.I*energy_transfer*time/hbar)
            weight_sum += weight
            spectral += weight*s.DiracDelta(delta-energy_transfer)
            reversed_spectral += weight*s.DiracDelta(delta+energy_transfer)
    require(s.simplify(direct_correlation-spectral_correlation) == 0,
            'Heisenberg correlation / signed spectral weights mismatch')
    require(s.simplify(weight_sum-medium_norm) == 0, 'Spectral weight normalization failed')
    require(s.simplify(spectral-reversed_spectral) != 0, 'Signed support mutation not exposed')

    # Delta(E-hbar*omega)=delta(omega-E/hbar)/hbar, hbar>0.
    n2, w = s.symbols('n2 w', positive=True)
    omega = s.Symbol('omega', real=True)
    energy = s.Symbol('energy', positive=True)
    delta_jacobian = s.integrate(s.DiracDelta(energy-hbar*omega),
                                (omega, -s.oo, s.oo))
    energy_line = s.simplify(2*s.pi*n2*w*delta_jacobian/hbar)
    frequency_line = n2*(2*s.pi*w)/hbar**2
    require(s.simplify(energy_line-frequency_line) == 0, 'Golden-rule Jacobian mismatch')
    mutations = {
        'drop_nuclear_coupling': s.simplify(direct-medium_norm/hbar**2) != 0,
        'transpose_instead_of_adjoint': s.simplify(
            direct-abs(coupling)**2*s.trace(bath*rho_b*bath.T)/hbar**2) != 0,
        'reverse_noncommuting_medium_order': s.simplify(
            direct-abs(coupling)**2*s.trace(rho_b*bath*bath.H)/hbar**2) != 0,
        'reverse_energy_transfer': s.simplify(spectral-reversed_spectral) != 0,
        'drop_hbar_jacobian': s.simplify(energy_line-2*s.pi*n2*w/hbar) != 0,
    }
    require(all(mutations.values()), 'An exposing mathematical mutation survived')
    return {
        'status': 'passed', 'kind': 'formal finite-dimensional witness, NOT material prediction',
        'fixture': {'B_dimensionless': str(bath), 'rho_medium': str(rho_b),
                    'g_energy_unit': str(coupling), 'H0': '0',
                    'nuclear_basis_order': ['initial', 'final'],
                    'medium_energies_arbitrary_common_unit': [str(e) for e in energies]},
        'direct_full_space_t2_coefficient': str(direct),
        'independent_reduced_t2_coefficient': str(reduced),
        'same_medium_doubled_g_t2_coefficient': str(doubled),
        'integrated_spectral_weight': str(weight_sum), 'spectral_lines': transitions,
        'mutations_rejected': {key: bool(value) for key, value in mutations.items()},
        'meaning': 'Same medium state/operator, different nuclear couplings, different transition probability coefficients. Neither an irreversible rate nor a value for Ni/H follows.',
        'golden_rule_limit': 'Energy/frequency normalization checked only; stationary continuum, weak coupling and physical nuclear matrix element applicability are NOT established by this finite example.',
    }


def admission_oracle(registration):
    changed = dict(registration['published_receipts']['nuclear-product-energy-ledger.json'])
    changed['sha256'] = '0'*64
    try:
        admit(changed)
    except ValueError:
        pass
    else:
        raise ValueError('Corrupt SHA256 mutation survived')
    omitted = dict(registration['published_receipts'])
    omitted.pop('nuclear-product-energy-ledger.json')
    try:
        require_receipts(omitted)
    except ValueError:
        pass
    else:
        raise ValueError('Missing receipt mutation survived')
    return {'sha256_mutation_rejected': True, 'required_receipt_omission_rejected': True}


def build_report(primary_dir=None):
    registration = json_bytes(REGISTRATION)
    require_receipts(registration['published_receipts'])
    bindings = []
    for group in ('primary_sources', 'original_nuclear_sources', 'implementation_dependencies'):
        for entry in registration[group]:
            resolved = admit(entry, primary_dir=primary_dir if group == 'primary_sources' else None)
            bindings.append({**entry, 'resolved_path': str(resolved), 'role': group})
    documents = {}
    for name, entry in registration['published_receipts'].items():
        documents[name] = json_bytes(admit(entry))
        bindings.append({**entry, 'role': 'published_receipt'})

    ledger = documents['nuclear-product-energy-ledger.json']
    measured = documents['measured-nuclear-channels-001.json']
    transport = documents['nuclear-energy-transport-001.json']
    reference = documents['full-common-reference-004/full-common-reference.json']
    summary = documents['full-common-reference-004/summary.json']
    finiteq = documents['finiteq-material-001/result.json']
    field = documents['full-field-aggregate-001/result.json']
    ionic = documents['ionic-gaussian-material-002/result.json']
    require(measured['existing_ledger_sha256'] == registration['published_receipts'][
        'nuclear-product-energy-ledger.json']['sha256'], 'Measured/energy provenance mismatch')
    require(documents['nuclear-oracle-command-002.json']['exit_code'] == 0,
            'Existing energy oracle did not exit successfully')
    require(documents['measured-nuclear-oracle-001.json']['status'] == 'pass',
            'Existing measured-channel oracle not passed')
    require(documents['nuclear-transport-oracle-002.json']['status'] == 'passed',
            'Existing transport oracle not passed')
    for report, keys in (
        (ledger, ('calculated_fusion_rate', 'calculated_heat_rate')),
        (measured, ('actual_material_rate', 'correlated_two_proton_amplitude',
                    'ni_receiver_intrinsic_matrix_elements', 'deposited_heat')),
        (transport, ('actual_material_heat_mev', 'actual_material_power_w')),
    ):
        require(all(report[key] is None for key in keys), 'Upstream scope changed: review required')
    require(summary['status'] == 'complete' and summary['full_material_total'] is True,
            'Prepared reference publication incomplete')
    require(reference['metadata']['complete_electronic_coverage'] is True,
            'Prepared electronic coverage incomplete')
    require(field['status'] == 'complete' and ionic['status'] == 'complete',
            'Required material component publication incomplete')

    # Authenticate actual material payloads, without replaying their construction.
    archives = [reference['archive'], finiteq['archive'], ionic['receipt']['archive']]
    for archive in archives:
        admit(archive)
        bindings.append({**archive, 'role': 'published_material_archive_raw_bytes'})

    oracle = formal_oracle()
    oracle['admission_mutations'] = admission_oracle(registration)
    products = [{key: row[key] for key in (
        'channel', 'mechanism', 'entrance', 'nuclear_products', 'radiation',
        'endpoint', 'q_atomic_ground_state', 'central_value_sign',
        'uncertainty_l1_upper_bound')} for row in ledger['reaction_rows']]
    isotopes = {label: {'A': entry['A'], 'Z': entry['Z'],
                        'nubase': entry['nubase'],
                        'isotope_abundance_scope': entry['isotope_abundance_scope']}
                for label, entry in ledger['isotope_inputs'].items()}
    material = {
        'prepared_reference': {
            'source_receipt': 'full-common-reference-004/full-common-reference.json',
            'source_count': summary['source_count'], 'slot_count': summary['slot_count'],
            'arrays': reference['arrays'],
            'complete_electronic_coverage': reference['metadata']['complete_electronic_coverage'],
            'held_density': reference['metadata']['held_density'],
            'held_gradient': reference['metadata']['held_gradient'],
            'limitations': reference['metadata']['limitations'],
            'nuclear_use': 'Prepared electronic/ionic reference and screening/force inputs only; neither nuclear pair wavefunction nor nuclear transition operator.'},
        'directed_finite_q_operator': {key: finiteq[key] for key in (
            'scope', 'source_dimensions', 'target_dimensions', 'q_fractional',
            'q_bohr_inv', 'weight', 'grid_measure_bohr3', 'ispin', 'limitations')},
        'field_moments': {
            'normalization': field['aggregate']['normalization'],
            'moments': field['aggregate']['moments'],
            'geometry_policy': field['aggregate']['geometry_policy'],
            'correction_scope': field['aggregate']['correction_scope'],
            'joint_scalar_limit': field['joint_scalar_limit'],
            'field_state': field['field_state'],
            'nuclear_use': 'Zero-field prepared moments, not a stationary nuclear-frequency spectral density or a radiative continuum.'},
        'ionic_preparation': {key: ionic[key] for key in (
            'preparation', 'covariance_published', 'dynamical_force_grid_modified',
            'propagated_trajectory')},
        'normalization_contract': {
            'electrons': 'Independent spin/mesh sectors; explicit slot weights, raw Gram/factor occupations, no extra degeneracy. Density bohr^-3, physical grid pairing dv in bohr^3.',
            'ionic': 'All 37 atom-order xyz coordinates (111); u bohr, positive density gradient Ry/bohr, unmassweighted Cq Ry/bohr^2, Kq maps density bohr^-3 to Ry; time tau=t Ry/hbar.',
            'field': 'Uniform periodic q=0 a=e A_SI/hbar in bohr^-1, Pi canonical; electric field E_Ha=-da/dtau/2. This vector-potential A is not DPV1988 electric-field A. No cavity frequency, nuclear receiver or vacuum inferred.',
            'statistics': 'A response/commutator is not automatically an unsymmetrized fluctuation spectrum. A moving prepared state need not be stationary or obey equilibrium FDT.'},
        'not_claimed': 'Whole many-region dynamics/comparison or a complete nuclear-frequency correlator; these are distinct from the component receipts consumed here.',
    }
    return {
        'schema': 'P254.A28.nuclear-medium-assessment.v1', 'status': 'source_operator_assessment_complete',
        'node': 'N', 'earned_scope': registration['scope'],
        'registration': {'path': str(REGISTRATION.relative_to(ROOT)), 'sha256': digest(REGISTRATION)},
        'consumer': {'path': str(Path(__file__).resolve().relative_to(ROOT)), 'sha256': digest(Path(__file__).resolve())},
        'verified_raw_byte_bindings': bindings,
        'mechanism_assessments': registration['mechanism_assessments'],
        'material_input_assessment': material,
        'reused_energy_product_routes': products,
        'reused_isotope_quantum_numbers_and_decay_properties': isotopes,
        'reused_decay_topology': ledger['decay_topology'],
        'ordinary_measured_channels': {
            'tables': measured['tables'], 'domains': measured['source_domains'],
            'printed_q_discrepancy': measured['source_printed_q_comparison'],
            'sequential_energy': measured['sequential_energy'],
            'decay_scope': measured['decay_scope'],
            'material_entrance_distribution': measured['actual_material_entrance_distribution'],
            'missing_inputs': measured['missing_inputs']},
        'conditional_product_transport': {
            'photon_tables': transport['photon_tables'], 'stopping_tables': transport['stopping_tables'],
            'source_gaps': transport['source_gaps'],
            'primary_table_limitations': transport['primary_table_limitations'],
            'missing_material_inputs': transport['missing_material_inputs'],
            'bridge': registration['conditional_energy_bridge'],
            'no_q_to_heat': True},
        'formal_conditional_bridge': registration['formal_bridge'], 'new_oracle': oracle,
        'existing_oracles_reused_not_replayed': [
            'nuclear-oracle-command-002.json', 'measured-nuclear-oracle-001.json',
            'nuclear-transport-oracle-002.json'],
        'actual_material_nuclear_rate_s_inv': None, 'actual_material_heat_w': None,
        'correlated_two_proton_amplitude': None, 'nuclear_receiver_matrix_element': None,
        'scope_decision': {
            'N_assessment': 'Complete exact source/operator/input/product assessment with missing prerequisites as its legitimate observable; not held open solely for absence of a material nuclear rate.',
            'after_D_and_C': 'Even an executed collective material response and matched controls do not supply pair preparation, a channel-resolved nuclear amplitude, nuclear-frequency matching or exclusive product spectra. Those require external physical inputs/construction, not rerunning D/C.',
            'full_LENR_parent': 'Not completed: predictive isotope-resolved nuclear rates, measured products and finite-sample thermal power still require the listed physical inputs.',
            'logical_limit': 'No missing operator is equated to a universal no-go or an actual zero physical rate.'},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--primary-source-dir', type=Path,
                        help='Relocate the six exact primary PDF/text files; hashes remain mandatory')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    # Fail closed before publishing a successful receipt on any mismatch/oracle failure.
    require(not args.output.exists(), 'Use a new receipt path; do not overwrite evidence')
    report = build_report(args.primary_source_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': report['status'], 'output': str(args.output),
                      'oracle': report['new_oracle']['status'],
                      'material_rate': None, 'material_heat': None}))


if __name__ == '__main__':
    main()
