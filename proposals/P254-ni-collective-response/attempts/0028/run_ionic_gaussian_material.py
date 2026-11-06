#!/usr/bin/env python3
"""Publish conditional sudden-release ionic covariance from the complete A28 lattice.

Root executes after the independent oracle. No propagation, eigenmode repair,
spin/k reweighting or equilibrium claim. Missing positive trap data are reported
as missing preparation input, never replaced by a stabilized dynamical matrix.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
for directory in (ROOT/'src', HERE.parent/'0025', HERE.parent/'0024'):
    sys.path.insert(0, str(directory))

import numpy as np
from material_run_contract import binding, bound_json, canonical, recover_block, commit_block
from material_reference_output import load_arrays
from substrate_framework.qe_saved_reference import load_saved_reference
from substrate_framework.ionic_gaussian_preparation import prepare_local_trap_vacuum, MissingTrapPreparation


def source_paths():
    return [HERE/'ionic-gaussian-registration.json', Path(__file__),
            HERE/'run_ionic_gaussian_oracle.py',
            *[ROOT/'src/substrate_framework'/name for name in (
                'ionic_gaussian_preparation.py', 'physical_bilinear_preparation.py',
                'inertial_fermion_dynamics.py', 'material_modes.py',
                'material_force_lattice.py', 'qe_saved_reference.py')],
            HERE.parent/'0025/material_run_contract.py',
            HERE.parent/'0025/material_reference_output.py',
            HERE.parent/'0024/material_pilot_io.py']


def explicit_binding(path, digest):
    record = binding(path)
    if record['sha256'] != digest:
        raise ValueError('explicit input SHA256 differs: '+str(path))
    return record


def authenticated_reference(registration):
    authority = registration['reference']
    receipt = bound_json(authority)
    archive = Path(receipt['archive']['path'])
    recovered = recover_block(archive.parent, archive.stem,
                              receipt['recipe_sha256'], receipt['metadata'])
    if recovered != receipt:
        raise ValueError('complete common reference receipt/archive differs')
    metadata = receipt['metadata']
    if (metadata['complete_electronic_coverage'] is not True or len(metadata['slots']) != 128
            or len(metadata['source_records']) != 20):
        raise ValueError('complete twenty-source/128-slot common reference required')
    manifest = bound_json(metadata['recipe'])
    if metadata['source_records'] != manifest['sources']:
        raise ValueError('source records disagree with exact full-reference recipe')
    preparation = bound_json(manifest['preparation_terminal'])
    bloch_binding = preparation['recipe']['inputs']['bloch_terminal']
    bloch = bound_json(bloch_binding)
    xml = bloch['recipe']['inputs']['xml']
    original = load_saved_reference(xml)
    arrays = load_arrays(receipt)
    species = {item.name: item.mass_amu for item in original.species}
    masses = np.asarray([species[name] for name in original.atom_names])
    for key, expected in [('atomic_masses_amu', masses), ('positions_bohr', original.positions_bohr),
                          ('cell_bohr', original.cell_bohr)]:
        if not np.array_equal(arrays[key], expected):
            raise ValueError('original XML mass/coordinate identity differs: '+key)
    shape = tuple(registration['grid_shape'])+(original.nat, 3, original.nat, 3)
    if any(arrays[name].shape != shape for name in
           ('common_force_lattice_ry_bohr2', 'common_q_grid_ry_bohr2')):
        raise ValueError('complete original seven-axis 4^3 force lattice required')
    if tuple(original.monkhorst_pack['nk'+str(a)] for a in (1,2,3)) != tuple(registration['grid_shape']):
        raise ValueError('original mesh differs from registered 64-cell convention')
    identity = dict(xml=xml, atom_indices_one_based=list(original.atom_indices),
                    atom_names=list(original.atom_names), coordinate_identity=registration['coordinate_identity'],
                    preparation_terminal=manifest['preparation_terminal'], bloch_terminal=bloch_binding)
    return arrays, identity, metadata['recipe']


def publish(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--oracle', type=Path, required=True)
    parser.add_argument('--oracle-sha256', required=True)
    parser.add_argument('--initial-ions', type=Path, required=True)
    parser.add_argument('--initial-ions-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    sources = [binding(path) for path in source_paths()]
    registration_binding = binding(HERE/'ionic-gaussian-registration.json')
    registration = bound_json(registration_binding)
    oracle_binding = explicit_binding(args.oracle, args.oracle_sha256)
    oracle = bound_json(oracle_binding)
    if oracle.get('status') != 'passed' or oracle.get('schema') != 'P254-ionic-gaussian-oracle-v1':
        raise ValueError('passed ionic Gaussian independent oracle required')
    if oracle['source_bindings'] != sources:
        raise ValueError('oracle does not bind current ionic preparation consumers/registration')
    ions_binding = explicit_binding(args.initial_ions, args.initial_ions_sha256)
    ions = bound_json(ions_binding)
    arrays, original, recipe_binding = authenticated_reference(registration)
    n = 3*len(arrays['atomic_masses_amu'])
    if ions['coordinate_identity'] != registration['coordinate_identity']:
        raise ValueError('initial means have wrong original coordinate identity')
    for key in ('u_bohr', 'p_canonical'):
        values = np.asarray(ions[key])
        if values.shape != (n,) or not np.array_equal(values, np.zeros(n)):
            raise ValueError('declared sudden-release preparation requires supplied exact zero '+key)
    metadata = dict(schema='P254-ionic-gaussian-material-v1', registration=registration_binding,
        hypothesis=registration, reference=registration['reference'], full_reference_recipe=recipe_binding,
        original_coordinate_and_mass_identity=original, initial_means=ions_binding,
        oracle=oracle_binding, source_bindings=sources, time_ry_inverse=0.,
        electronic_ionic_cross=dict(ordered='identically zero', reverse_ordered='identically zero',
            symmetrized='identically zero', commutator='identically zero',
            authority='explicit original-electronic-state tensor-product ionic-trap-vacuum hypothesis'),
        covariance_storage='exact factorized local block plus Kronecker delta in cell or q/-q indices; not a truncation',
        dynamical_force_grid_modified=False, extra_spin_or_k_measure=False)
    args.output.mkdir(parents=True, exist_ok=False)
    try:
        state = prepare_local_trap_vacuum(arrays['common_force_lattice_ry_bohr2'],
            arrays['atomic_masses_amu'], state_authority=registration_binding['sha256'],
            reality_relative_tolerance=registration['reality_relative_tolerance'])
    except MissingTrapPreparation as error:
        publish(args.output/'result.json', dict(status='missing-preparation-input',
            failure=str(error), metadata=metadata, covariance_published=False))
        return 1
    payload = state.archive_arrays()
    payload.update(positions_bohr=arrays['positions_bohr'], cell_bohr=arrays['cell_bohr'],
                   atomic_masses_amu=arrays['atomic_masses_amu'])
    metadata['diagnostics'] = dict(
        stiffness_imaginary_max=float(np.max(np.abs(state.stiffness_raw.imag))),
        stiffness_real_min=float(np.min(state.stiffness_raw.real)),
        exact_ground_state_uncertainty='Vuu_aa*Vpp_aa=1/4; ordered up=+i/2',
        ncells=state.ncells, nion_coordinates_per_cell=state.ndof,
        initial_force_nonzero=bool(np.any(arrays['force_at_origin_ry_bohr'] != 0)),
        initial_force_norm=float(np.linalg.norm(arrays['force_at_origin_ry_bohr'])),
        q_conjugacy='q_minus_index labels adjoints, not independent real canonical coordinates',
        mean_CCR='i/64; mean_ordered must NOT be passed as canonical IonicPreparation')
    identity = hashlib.sha256(canonical(metadata)).hexdigest()
    receipt = commit_block(args.output, 'ionic-gaussian', payload, metadata, identity, 64 << 20)
    publish(args.output/'result.json', dict(status='complete', covariance_published=True,
        receipt=receipt, preparation='explicit conditional sudden-release trap vacuum',
        dynamical_force_grid_modified=False, propagated_trajectory=False))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
