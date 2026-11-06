#!/usr/bin/env python3
"""Finite saved-primary intake, complete mesh ledger, ONE selected preparation.

No Native executable/library, SCF, PH or dynamics is invoked. Source header/PW
admission is exhaustive, but only the requested source/target Gram and D are
instantiated. Run allocation and execution belong to Root.
Use --require-orbital-conjugation for dispersed spin1 mesh(3,0,0). The same
explicit manifest/binding CLI also accepts selected unitary slots of controls.
"""
import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np

from substrate_framework.material_control_inputs import (
    load_control_inputs, iter_control_source_admissions,
    control_population_specification, control_population_arrays,
)
from substrate_framework.material_reference_state import prepare_source
from substrate_framework.qe_saved_native import load_wavefunctions
from substrate_framework.qe_bloch_frame import admit_source
from substrate_framework.qe_real_bloch_grid import transport_real_state


def binding(path, expected_sha256=None):
    path = Path(path).resolve()
    with path.open('rb') as stream:
        digest = __import__('hashlib').file_digest(stream, 'sha256').hexdigest()
    if expected_sha256 is not None and digest != expected_sha256:
        raise ValueError('explicit input SHA256 mismatch: ' + str(path))
    return dict(path=str(path), bytes=path.stat().st_size, sha256=digest)


def save_json(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
    return binding(path)


def save_arrays(path, arrays):
    with path.open('xb') as stream:
        np.savez(stream, **arrays)
    return binding(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('case-manifest', 'native-receipt', 'admission'):
        parser.add_argument('--' + name, type=Path, required=True)
        parser.add_argument('--' + name + '-sha256', required=True)
    parser.add_argument('--spin', type=int, choices=(1, 2), required=True)
    parser.add_argument('--mesh-index', nargs=3, type=int, required=True)
    parser.add_argument('--require-orbital-conjugation', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    authorities = [binding(getattr(args, name), getattr(args, name + '_sha256'))
                   for name in ('case_manifest', 'native_receipt', 'admission')]
    control = load_control_inputs(*authorities)
    grid, ref = control.grid, control.reference
    selected = next((s for s in grid.slots if s.ispin == args.spin
                     and s.mesh_index == tuple(args.mesh_index)), None)
    if selected is None:
        raise ValueError('explicit selected slot is absent from the full mesh')
    if args.require_orbital_conjugation and not selected.orbital_conjugation:
        raise ValueError('requested slot does not exercise orbital conjugation')
    registration = binding(Path(__file__).with_name('control-preparation-registration.json'))
    decision = 'A28 control-preparation-registration sha256=' + registration['sha256']
    descriptors, admissions = [], []
    for source, audit, records in iter_control_source_admissions(control):
        descriptors.extend(records)
        admissions.append(dict(binding=dict(source.source), header=dict(ik=source.ik,
            ispin=source.ispin, rows=source.igwx, columns=source.nbnd, npol=source.npol,
            gamma_only=source.gamma_only, scalef=source.scalef), admission=dict(audit)))
        del source
    descriptors.sort(key=lambda d: (tuple(d['mesh_index']), d['spin']))
    coverage = save_json(output / 'mesh-source-admission.json', dict(
        schema='real-collinear-control-coverage-v1', provenance=control.provenance(),
        source_admissions=admissions, slots=descriptors,
        scope='all original source bytes/headers/rows/columns admitted; all mesh bases validated; not all-slot Gram/state execution'))
    specification = control_population_specification(control, descriptors,
        coverage_authority='complete source/slot admission sha256=' + coverage['sha256'],
        decision_authority=decision)
    ledger = control_population_arrays(control, specification)
    population_receipt = save_arrays(output / 'population-ledger.npz', ledger)
    index = next(i for i, slot in enumerate(grid.slots) if slot is selected)
    start, stop = specification.band_slices[index]
    population = specification.filling.populations[start:stop]
    # Same state on every route of this original spin/source orbit.
    for i, slot in enumerate(grid.slots):
        if (slot.ispin, slot.source_ik) == (selected.ispin, selected.source_ik):
            a, b = specification.band_slices[i]
            if tuple(population) != specification.filling.populations[a:b]:
                raise ValueError('source-orbit filling is not identical')
    source = load_wavefunctions(selected.source_binding)
    audit = admit_source(ref, selected, source)
    source_state = prepare_source(source.coefficients, population,
        source_identity=source.source['sha256'], basis_identity=audit['source_rows_sha256'], authority=decision)
    frame, target_state = transport_real_state(grid, selected, source, source_state, source_admission=audit)
    # Independent raw-frame recomputation of BOTH target S and symmetric-frame D.
    direct = prepare_source(frame.coefficients, population,
        source_identity=source.source['sha256'], basis_identity=frame.basis.identity, authority=decision)
    gram_error = float(np.max(np.abs(direct.gram - target_state.gram)))
    d_error = float(np.max(np.abs(direct.density_coefficients - target_state.density_coefficients)))
    tolerance = 512 * np.finfo(float).eps * max(frame.coefficients.shape)
    if gram_error > tolerance or d_error > tolerance:
        raise ValueError('direct target Gram/D differs from source conjugation')
    source_number = np.einsum('ij,ji->', source_state.density_coefficients, source_state.gram)
    target_number = np.einsum('ij,ji->', target_state.density_coefficients, target_state.gram)
    expected_number = source_number.conjugate() if selected.orbital_conjugation else source_number
    if abs(target_number - expected_number) > tolerance:
        raise ValueError('source/target number normalization differs')
    selected_receipt = save_arrays(output / 'selected-real-frame-state.npz', dict(
        coefficients=frame.coefficients, miller=frame.basis.miller, source_miller=source.miller,
        row_map=frame.basis.row_map, reciprocal_bohr_inv=frame.basis.reciprocal_bohr_inv,
        k_cart_bohr_inv=frame.basis.k_cart_bohr_inv, source_k_cart_bohr_inv=source.k_cart_bohr_inv,
        cell_bohr=ref.cell_bohr, positions_bohr=ref.positions_bohr,
        atom_names=np.asarray(ref.atom_names), gram=target_state.gram,
        density_coefficients=target_state.density_coefficients,
        density_eigenvalues=np.asarray(target_state.admissibility['density_eigenvalues']),
        source_gram=source_state.gram, source_density_coefficients=source_state.density_coefficients,
        populations=np.asarray(population), eigenvalues_ry=frame.eigenvalues_ry,
        raw_mv_occupations=frame.raw_occupations,
        number_pair=np.asarray([source_number, target_number]),
        direct_gram_max_error=np.asarray(gram_error), direct_d_max_error=np.asarray(d_error)))
    filling = specification.filling
    implementation = [binding(Path(__file__))]
    import substrate_framework.qe_real_bloch_grid as grid_module
    import substrate_framework.material_control_inputs as input_module
    implementation += [binding(Path(m.__file__)) for m in (grid_module, input_module)]
    result = dict(schema='material-control-preparation-v1', status='completed',
        provenance=control.provenance(), registration=registration, implementation=implementation,
        coverage=coverage, population_ledger=population_receipt, selected_arrays=selected_receipt,
        selected=dict(ispin=selected.ispin, mesh_index=selected.mesh_index, source_ik=selected.source_ik,
            orbital_conjugation=selected.orbital_conjugation, basis_identity=frame.basis.identity,
            frame_identity=frame.gauge_identity, rows=source.igwx, columns=source.nbnd,
            weight=selected.weight, gram_max_error=gram_error, d_max_error=d_error,
            direct_comparison_allowance=tolerance),
        mesh_slots=len(grid.slots), original_sources=len(admissions),
        policy=asdict(specification.policy), neutral_target=filling.target_mean_number,
        weighted_number=filling.weighted_number, spin_numbers=filling.spin_numbers,
        chemical_potential_interval_ry=filling.chemical_potential_interval,
        claim='Complete source/header/PW/mesh/population ledger and one full-row/full-column real-policy preparation; not all-slot physical observables',
        limits=['stationary-primary only; no equivalence to current Gamma continuation asserted',
            'raw MV remains numerical smearing metadata, not probabilities or physical temperature',
            'conditional saved-label preparation is not reconstructed-H0 equilibrium',
            'no Native compute, dynamics, force/current material result, convergence or ranking'])
    receipt = save_json(output / 'result.json', result)
    print(json.dumps(dict(result=receipt, case=control.case, mesh_slots=len(grid.slots),
                          original_sources=len(admissions)), sort_keys=True))


if __name__ == '__main__':
    main()
