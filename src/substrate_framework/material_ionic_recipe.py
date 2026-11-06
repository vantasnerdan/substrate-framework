"""A27 actual ionic input/recipe authority; no density or electronic rebuild.

Requires the existing A24/A25 control modules on PYTHONPATH. The immutable
recipe binds original XML/NCUPFs/stdout and the settled A27 authorities. Only
ordinary periodic Ni31H6 is admitted, with an explicit numerical sum choice.
"""
from dataclasses import dataclass
import xml.etree.ElementTree as ET

import numpy as np
from material_run_contract import bound_json
from .qe_saved_reference import load_saved_reference
from .qe_nc_upf import read_nc_upf
from .qe_native_force_ledger import load_native_force_components
from .qe_ionic_ewald import EwaldSumContract

KIND = 'material-ionic-linear'
SCHEMA = 'material-ionic-linear-v1'
KEY = 'whole-ionic-baseline'
COORDINATES = 'Ni31H6 original XML atoms 1..37; atom-major Cartesian xyz; bohr/Ry'
SCOPE = '37-ion/111-coordinate finite periodic NC ion-ion scalar and gradient only; not total force or a prepared equilibrium'
EWALD = dict(alpha_bohr_inv2=1.0, real_cutoff_bohr=5.0, reciprocal_cutoff_bohr_inv=12.0)
RATIONALE = ('Fixed alpha=1 bohr^-2 is the initial candidate in QE force_ew; real radius '
             '5/sqrt(alpha) follows that force path. The reciprocal sphere is explicitly '
             '|G|<=12 bohr^-1 with all signed vectors. This is one declared finite scalar, '
             'not QE auto-selection or an error guarantee. No cutoff fitting or refinement is implicit.')
AUTHORITY_SHA = {
    'continuation_registration': 'cbb9f3698d8b432fb3825c640487095f9abef15004d36aeb1c311a38174cb92e',
    'root_decision': 'f0f36c5cb6e5d0c59ca9672d37be2a24f241a7017b03d8f0a4e063b484566b24',
    'implementation_registration': 'b59c25bd13d015ae17b2d81c8b98567f1358bacdc783b5c235f2f186aef268bc',
    'ionic_oracle': '40671f6cae4675faa8f1b956e51ba19f6d2699ef16fab2c8b7ca2ddab3d950ab',
    'native_ledger_oracle': '2e1fabcfa30c08ece4d988f35255cc839cf84fbe27da994b6fc4da5165b93f26',
}
INPUT_SHA = {
    'xml': '3831ae2b83231b02475931d50055415b1d9673cdb05edb4541d0f1115ec383da',
    'Ni': 'e6a7b028f02fccddb9158a150c7b16202af4741e8f5000e50df5ccc4b8122d3e',
    'H': 'a442cf8b2577354d11a7e126158ea8bb596e1323c018531990570425739e821a',
    'native_stdout': '6d9b4f4fd4c3b795b9bda39a8d4f3714f8476e928b0287607bac6eb488275b05',
}


def preload_parsers():
    """Preload XML/NPZ dispatch without reading material or evaluating science."""
    ET.fromstring('<source-inventory-preload/>')
    # material_run_contract -> material_pilot_io preloads actual NPZ dispatch.


def validate_recipe(recipe):
    """Authenticate additional A27 authorities under unchanged compound A25 IO."""
    inputs, settings = recipe['inputs'], recipe['settings']
    if recipe.get('kind', KIND) != KIND or recipe.get('library_binding') is not None:
        raise ValueError('ordinary ionic kind without a PBE backend required')
    if set(inputs) != {'xml', 'upfs', 'native_stdout'} or set(inputs['upfs']) != {'Ni', 'H'}:
        raise ValueError('exact original XML/two UPFs/Native stdout inputs required')
    records = dict(xml=inputs['xml'], native_stdout=inputs['native_stdout'], **inputs['upfs'])
    for name, record in records.items():
        if set(record) != {'path', 'bytes', 'sha256'} or record['sha256'] != INPUT_SHA[name]:
            raise ValueError('wrong original material source: '+name)
    if set(settings) != {'authorities', 'ewald_contract', 'native_evaluation_index', 'coordinate_identity'}:
        raise ValueError('explicit ionic recipe settings required')
    if settings['ewald_contract'] != EWALD or settings['coordinate_identity'] != COORDINATES:
        raise ValueError('only the declared finite sum/coordinate recipe is admitted')
    if type(settings['native_evaluation_index']) is not int or settings['native_evaluation_index'] != 0:
        raise ValueError('explicit original evaluation_index=0 required')
    authorities = settings['authorities']
    if set(authorities) != set(AUTHORITY_SHA):
        raise ValueError('all original A27 authorities required')
    documents = {}
    for name, digest in AUTHORITY_SHA.items():
        if authorities[name]['sha256'] != digest:
            raise ValueError('wrong original A27 authority: '+name)
        documents[name] = bound_json(authorities[name])
    if documents['native_ledger_oracle']['ledger']['source'] != inputs['native_stdout']:
        raise ValueError('comparison stdout differs from admitted original ledger')
    return documents


def validate_grant(grant):
    validate_recipe(grant)
    receipts = grant['resume_receipts']
    if (not isinstance(receipts, dict) or not set(receipts) <= {KEY}
            or (not grant['resume'] and receipts)):
        raise ValueError('fresh uses {}; resume explicitly binds the completed whole-baseline receipt')


@dataclass(frozen=True)
class IonicInputs:
    reference: object
    valence_charges: np.ndarray
    sources: tuple
    native_ledger: object


def load_ionic_inputs(inputs, settings):
    """Decode actual saved geometry and NCUPF VALENCE, plus original force ledger.

    No density, WFC, local-PBE construction or radial interpolation table is
    needed. read_nc_upf retains its original NC admissibility checks.
    """
    ref = load_saved_reference(inputs['xml'])
    if (ref.nat != 37 or ref.atom_names != ('Ni',)*31+('H',)*6
            or ref.atom_indices != tuple(range(1, 38))):
        raise ValueError('full original Ni31H6 atom order required')
    charges_by_species, sources = {}, [dict(ref.source)]
    for species in ref.species:
        record = inputs['upfs'][species.name]
        upf = read_nc_upf(record['path'], expected_bytes=record['bytes'], expected_sha256=record['sha256'])
        if upf.element != species.name:
            raise ValueError('UPF/XML element mismatch')
        charges_by_species[species.name] = upf.z_valence
        sources.append(dict(path=upf.path, bytes=str(len(upf.original_bytes)), sha256=upf.sha256))
    charges = np.asarray([charges_by_species[name] for name in ref.atom_names], dtype=float)
    charges.flags.writeable = False
    species_order = [species.name for species in ref.species]
    types = tuple(species_order.index(name)+1 for name in ref.atom_names)
    ledger = load_native_force_components(inputs['native_stdout'], atom_type_indices=types,
                                          evaluation_index=settings['native_evaluation_index'])
    return IonicInputs(ref, charges, tuple(sources), ledger)


def block_metadata(recipe, data):
    ledger = data.native_ledger
    return dict(schema=SCHEMA, kind=KIND, scope=SCOPE, coordinate_identity=COORDINATES,
                atom_indices=list(data.reference.atom_indices), atom_names=list(data.reference.atom_names),
                sources=list(data.sources), inputs=recipe['inputs'], settings=recipe['settings'],
                energy_term_order=['real', 'reciprocal', 'self', 'background'],
                energy_units='Ry', force_units='Ry/bohr', cell_vectors='columns',
                coordinate_order='atom-major Cartesian xyz', finite_sum_rationale=RATIONALE,
                native=dict(source=dict(ledger.source), atom_type_indices=list(ledger.atom_type_indices),
                            evaluation_index=ledger.evaluation_index, evaluation_count=ledger.evaluation_count,
                            summary_line=ledger.summary_line, summary_literal=ledger.summary_literal,
                            ion_ion=ledger.blocks['ion_ion'].as_dict()),
                comparison_scope='rounded original ion-ion block only; no processed total or cancellation',
                recovery_scope='whole completed baseline only; unfinished Ewald sums are lost',
                curvature_scope='use only linear F_II with existing joint C; no displaced Ewald curvature')


def sum_contract(settings):
    return EwaldSumContract(**settings['ewald_contract'])
