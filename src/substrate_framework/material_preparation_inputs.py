"""Authenticate complete Bloch coverage and the fixed A27 preparation decision.

Requires attempts/0024 and0025 runtime modules on PYTHONPATH. No background,
H0 or accepted functional is regenerated; material decoding occurs only when
these entry points are called inside an admitted bounded producer.
"""
from dataclasses import asdict
import hashlib

from material_run_contract import bound_json, canonical
from .material_bloch_output import load_saved_bloch_coverage, preload_parsers
from .material_preparation import (FullGridSector, MaterialPreparationPolicy,
                                   SourceBinding, specify_material_preparation)
from .qe_nc_upf import read_nc_upf

CONTINUATION = 'cbb9f3698d8b432fb3825c640487095f9abef15004d36aeb1c311a38174cb92e'
IMPLEMENTATION = 'b59c25bd13d015ae17b2d81c8b98567f1358bacdc783b5c235f2f186aef268bc'
ORACLE = 'd9b80bb87772056394f1b2100fc871204f3d048d882271864927fad129c791c5'
FACTOR_TERMINAL = 'deb800eb26059fac797e497435b8ceddf869a9c4362fd26dc3d012804c207a2f'
RECIPE_NAME = 'neutral-free-M-saved-label-exact-block-empty-complement-v1'
KEYS = ('neutral-populations', 'spin1-gamma-preparation')


def validate_settings(inputs, settings):
    if set(inputs) != {'bloch_terminal', 'factor_terminal'}:
        raise ValueError('completed authenticated Bloch and original A25 factor terminals required')
    if set(settings) != {'continuation_registration', 'implementation_registration',
                         'preparation_oracle', 'preparation_recipe', 'selected_spin'}:
        raise ValueError('exact bound continuation/decision/oracle and selected recipe required')
    for key, digest in (('continuation_registration', CONTINUATION),
                        ('implementation_registration', IMPLEMENTATION), ('preparation_oracle', ORACLE)):
        if settings[key]['sha256'] != digest:
            raise ValueError('different settled preparation authority: ' + key)
        bound_json(settings[key])
    if settings['preparation_recipe'] != RECIPE_NAME or settings['selected_spin'] != 1:
        raise ValueError('only registered complete neutral / original spin1 Gamma consumer admitted')
    if inputs['factor_terminal']['sha256'] != FACTOR_TERMINAL or str(inputs['factor_terminal']['bytes']) != '1772580':
        raise ValueError('original admitted A25 factor terminal required')


def validate_grant(grant):
    validate_settings(grant['inputs'], grant['settings'])
    if grant.get('library_binding') is not None:
        raise ValueError('no PBE/backend library is invoked by preparation')
    receipts = grant['resume_receipts']
    if (not isinstance(receipts, dict) or not set(receipts) <= set(KEYS)
            or (not grant['resume'] and receipts)):
        raise ValueError('explicit immutable completed-boundary receipts required')


def preparation_policy(settings, coverage_identity):
    authority = 'A27 implementation-registration sha256=' + settings['implementation_registration']['sha256']
    return MaterialPreparationPolicy(authority,
        authority + '; actual atom-resolved NC UPF valence sum equals saved nelec',
        authority + '; Q=Psi S^-1/2; diagonal original saved Ry labels in Q',
        authority + '; exact supplied equality; scalar number-constrained boundary block',
        authority + '; XML-order unitary identity-first canonical extension, no little-group assertion',
        authority + '; independent explicit spin/full-k quasifree sectors',
        authority + '; exactly empty initial orthogonal complement, all future PW rows retained',
        'authenticated complete20-source/128-slot Bloch recipe sha256=' + coverage_identity)


def load_population_inputs(inputs, settings):
    """Load completed coverage, bound UPFs and saved metadata, then fill mean N."""
    validate_settings(inputs, settings)
    coverage = load_saved_bloch_coverage(inputs['bloch_terminal'])
    reference = coverage.grid.reference
    if (len(coverage.slot_records) != 128 or len(coverage.source_receipts) != 20
            or reference.nbnd_by_spin != (346, 346) or reference.nat != 37
            or reference.atom_names != ('Ni',) * 31 + ('H',) * 6):
        raise ValueError('requires actual complete Ni31H6 coverage, not fixture or partial Gamma')
    original = coverage.recipe['inputs']
    upfs = {name: read_nc_upf(record['path'], expected_bytes=record['bytes'],
                            expected_sha256=record['sha256']) for name, record in original['upfs'].items()}
    sectors = tuple(FullGridSector(r.ispin, r.mesh_index, r.source_ik,
        SourceBinding.from_mapping(r.source_binding), r.gauge_identity, r.basis_identity,
        r.n_pw, coverage.coverage_identity + ':' + r.gauge_identity) for r in coverage.slot_records)
    specification = specify_material_preparation(reference, upfs, sectors,
        policy=preparation_policy(settings, coverage.coverage_identity))
    if abs(specification.filling.weighted_number - reference.nelec) > 1e-10:
        raise ValueError('conditional full-grid weighted neutral number failed')
    return coverage, specification


def load_factor_terminal(inputs, coverage):
    terminal = bound_json(inputs['factor_terminal'])
    recipe = terminal['recipe']
    if (terminal['status'] != 'complete' or terminal['exit_code'] != 0
            or recipe['kind'] != 'physical-factors'
            or hashlib.sha256(canonical(recipe)).hexdigest() != terminal['recipe_sha256']
            or len(terminal['outputs']) != 2):
        raise ValueError('complete original two-spin factor recipe required')
    if recipe['inputs']['xml'] != coverage.recipe['inputs']['xml']:
        raise ValueError('Bloch and accepted factor original XML sources differ')
    if {o['key'] for o in terminal['outputs']} != {'spin1-functional', 'spin2-functional'}:
        raise ValueError('original factor spin outputs differ')
    return terminal


def population_metadata(inputs, settings, coverage, specification):
    return dict(schema='material-neutral-populations-v1', inputs=inputs, settings=settings,
        coverage_identity=coverage.coverage_identity, policy=asdict(specification.policy),
        neutral_target=specification.filling.target_mean_number,
        weighted_number=specification.filling.weighted_number,
        spin_numbers=[list(pair) for pair in specification.filling.spin_numbers],
        chemical_potential_interval=list(specification.filling.chemical_potential_interval),
        boundary_log_fugacity=specification.filling.boundary_log_fugacity,
        energy_units='Ry', number_units='electrons per primitive cell',
        degeneracy='exact supplied energy equality; no near-degeneracy clustering',
        state_scope='conditional mean-N independent full-k/spin quasifree; not exact-N or H0 equilibrium',
        complement='initially empty; full-row future actions retained',
        symmetry='declared canonical unitary extension; little-group invariance not established',
        reference_binding=asdict(specification.reference_binding),
        charge_bindings=[[name, asdict(record)] for name, record in specification.charge_bindings],
        valence_charges=[list(pair) for pair in specification.valence_charges])
