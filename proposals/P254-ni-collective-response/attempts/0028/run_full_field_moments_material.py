#!/usr/bin/env python3
"""Original-source full-column physical coefficients and complete slot aggregate.

No new populations/Gram matrices. Run oracle first. Source checkpoints are
unweighted and immutable. Aggregate requires all 20 original checkpoints and
explicitly retains the declared saved-geometry approximation for 128 slots.
Estimate reads admitted metadata/background, but applies no PW operator.
"""
import argparse
from collections.abc import Mapping
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE.parent/'0025'), str(HERE.parent/'0024')]

import numpy as np
from material_run_contract import binding, bound_json, recover_block, runtime_context
from material_reference_output import load_arrays
from substrate_framework.material_bloch_output import preload_parsers
from substrate_framework.material_preparation_output import load_saved_material_preparation
from substrate_framework.material_bloch_admission import restore_accepted, wave_bindings
from substrate_framework.material_reference_state import source_population_indices, restore_source
from substrate_framework.qe_saved_reference import load_saved_reference
from substrate_framework.qe_saved_native import load_wavefunctions
from substrate_framework.qe_bloch_grid import build_bloch_grid
from substrate_framework.qe_bloch_frame import transport_bloch_frame
from substrate_framework.qe_bloch_dyads import BlochDyadCoupling
from substrate_framework.qe_material_hamiltonian import MaterialHamiltonian
from substrate_framework import material_physical_field_moments as moments_module
from substrate_framework import qe_material_field_derivatives as derivative_module
from substrate_framework import material_reference_transport as transport_module


def pack(value):
    if isinstance(value, Mapping):
        return {k: pack(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [pack(v) for v in value]
    if isinstance(value, np.ndarray) or isinstance(value, (complex, np.complexfloating)):
        a = np.asarray(value)
        return dict(real=a.real.tolist(), imag=a.imag.tolist())
    if isinstance(value, np.generic):
        return value.item()
    return value


def unpack(value):
    if isinstance(value, dict):
        if set(value) == {'real', 'imag'}:
            return np.asarray(value['real'])+1j*np.asarray(value['imag'])
        return {k: unpack(v) for k, v in value.items()}
    if isinstance(value, list):
        return [unpack(v) for v in value]
    return value


def receipt_arrays(record):
    receipt = bound_json(record)
    archive = Path(receipt['archive']['path'])
    actual = recover_block(archive.parent, archive.stem, receipt['recipe_sha256'], receipt['metadata'])
    if actual != receipt:
        raise ValueError('authenticated completed source receipt differs')
    return load_arrays(actual), actual


def implementation_bindings():
    return [binding(p) for p in (Path(__file__), Path(moments_module.__file__),
                                Path(derivative_module.__file__), Path(transport_module.__file__),
                                HERE/'run_full_field_moments_oracle.py',
                                HERE/'run_field_derivative_oracle.py')]


def load_inputs(aggregate_path):
    aggregate_binding = binding(aggregate_path)
    if aggregate_binding['sha256'] != '9b7be6a210219069f8f46aae38c7ae1ba0b5c150cc35d953db8c6632084678f9':
        raise ValueError('settled full common reference bytes differ')
    aggregate = bound_json(aggregate_binding)
    manifest_binding = aggregate['metadata']['recipe']
    manifest = bound_json(manifest_binding)
    saved = load_saved_material_preparation(manifest['preparation_terminal'])
    bloch = bound_json(saved.recipe['inputs']['bloch_terminal'])['recipe']
    original = bloch['inputs']
    reference = load_saved_reference(original['xml'])
    grid = build_bloch_grid(reference, wave_bindings(original), symmetry_policy=bloch['settings']['symmetry_policy'])
    keys = [(r['spin'], r['source_ik']) for r in manifest['sources']]
    if (len(keys) != 20 or len(set(keys)) != 20 or set(keys) != set(wave_bindings(original))
            or len(grid.slots) != 128):
        raise ValueError('complete 20-source/128-slot authority required')
    return aggregate_binding, manifest_binding, manifest, saved, original, grid


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('estimate', 'source', 'aggregate'), required=True)
    parser.add_argument('--full-reference', type=Path, default=HERE/'full-common-reference-004/full-common-reference.json')
    parser.add_argument('--output', type=Path, required=True, help='fresh directory')
    parser.add_argument('--output-max-bytes', type=int, default=4194304)
    parser.add_argument('--oracle', type=Path, required=True)
    parser.add_argument('--spin', type=int)
    parser.add_argument('--source-ik', type=int)
    parser.add_argument('--source-checkpoint', type=Path, action='append', default=[])
    parser.add_argument('--geometry-policy', choices=('declared-geometric',))
    parser.add_argument('--row-block', type=int, default=256)
    args = parser.parse_args()
    if args.row_block < 1 or args.output_max_bytes < 65536:
        raise ValueError('positive row block and at least 65536 output bytes required')
    started = time.monotonic()
    preload_parsers()
    registration_binding = binding(HERE/'full-field-moments-registration.json')
    registration = bound_json(registration_binding)
    for parent in registration['sources']:
        inherited = bound_json(parent)
        for source in inherited['source_bindings']:
            actual = binding(source['path'])
            if actual['sha256'] != source['sha256'] or int(actual['bytes']) != source['size_bytes']:
                raise ValueError('original registered source bytes changed: '+source['path'])
    implementation = implementation_bindings()
    oracle_binding = binding(args.oracle)
    oracle = bound_json(oracle_binding)
    if (oracle.get('status') != 'passed' or oracle.get('registration') != registration_binding
            or oracle.get('implementation_sources') != implementation):
        raise ValueError('passed oracle bound to these exact implementation bytes required')
    aggregate_binding, manifest_binding, manifest, saved, original, grid = load_inputs(args.full_reference)
    report = dict(schema='P254-full-field-moments-material-v1', mode=args.mode,
                  registration=registration_binding, oracle=oracle_binding,
                  implementation_sources=implementation, full_reference=aggregate_binding,
                  manifest=manifest_binding, runtime_context=runtime_context(),
                  field_state=None, h0_transport_supported=False,
                  joint_scalar_limit=registration['joint_scalar_limit'])
    if args.mode == 'aggregate':
        if args.geometry_policy is None or len(args.source_checkpoint) != 20:
            raise ValueError('all 20 source checkpoint paths and explicit geometry policy required')
        sources, checkpoints = {}, []
        for path in args.source_checkpoint:
            record = binding(path)
            checkpoint = unpack(bound_json(record))
            if (checkpoint.get('mode') != 'source' or checkpoint.get('status') != 'complete'
                    or checkpoint['registration'] != registration_binding
                    or checkpoint['implementation_sources'] != implementation
                    or checkpoint['manifest'] != manifest_binding):
                raise ValueError('source checkpoint construction/implementation/state authority differs')
            key = (checkpoint['spin'], checkpoint['source_ik'])
            if key in sources or checkpoint['source_binding'] != wave_bindings(original)[key]:
                raise ValueError('duplicate or mismatched original source checkpoint')
            sources[key] = checkpoint['moments']
            checkpoints.append(record)
        report.update(source_checkpoints=checkpoints,
                      aggregate=moments_module.aggregate_slots(grid, sources, geometry_policy=args.geometry_policy))
    else:
        background, background_receipt = restore_accepted(original)
        report['background_receipt'] = background_receipt
        if args.mode == 'estimate':
            columns = 29*sum(sum(name == species.name for name in grid.reference.atom_names)
                             *background.inputs.tables[species.name].upf.nh
                             for species in grid.reference.species)
            counts = []
            for source in manifest['sources']:
                point, energies, _ = grid.reference.select(source['source_ik'], source['spin'])
                p, n = int(point.npw), int(energies.size)
                counts.append(dict(spin=source['spin'], source_ik=source['source_ik'], pw_rows=p,
                                   raw_columns=n, frame_bytes=16*p*n, wide_projector_columns=columns,
                                   overlap_complex_multiply_adds=p*n*columns,
                                   kinetic_complex_multiply_adds=p*n*n))
            report['estimate'] = dict(sources=counts,
                original_frame_bytes=sum(c['frame_bytes'] for c in counts),
                full_slot_frame_bytes=13137492672,
                overlap_complex_multiply_adds=sum(c['overlap_complex_multiply_adds'] for c in counts),
                kinetic_complex_multiply_adds=sum(c['kinetic_complex_multiply_adds'] for c in counts),
                overlap_workspace_upper_bytes=16*max(c['raw_columns'] for c in counts)*columns,
                row_workspace='O(row_block*(29*max_nh+raw_columns)); existing species beta bank retained',
                radial_work='once per species/PW row: 3 exact retained radial Simpson contractions per l',
                output_real_scalars_per_source=2*(29+15*grid.reference.nat),
                output_cap_bytes_per_invocation=args.output_max_bytes,
                elapsed_time_prediction=None, throughput='not measured; no runtime extrapolation invented',
                full_pw_operator_applications=0, local_fft_count=0)
        else:
            key = (args.spin, args.source_ik)
            selected = [r for r in manifest['sources'] if (r['spin'], r['source_ik']) == key]
            if len(selected) != 1:
                raise ValueError('--spin and --source-ik must select one complete original source')
            record = selected[0]
            population, pairs = source_population_indices(grid, saved.population_ledger, *key)
            rep = next(s for _, s in pairs if s.operation.index == 1 and not np.any(s.reciprocal_shift))
            source = load_wavefunctions(rep.source_binding)
            frame = transport_bloch_frame(grid.reference, rep, source)
            if record['reuse_selected_gamma']:
                if key != (1, 1) or record['state'] is not None:
                    raise ValueError('selected Gamma source identity differs')
                arrays = {name: saved.selected_gamma[name] for name in
                          ('gram', 'populations', 'density_coefficients', 'density_eigenvalues')}
                state_authority = manifest['preparation_terminal']
            else:
                arrays, receipt = receipt_arrays(record['state'])
                metadata = receipt['metadata']
                if (metadata['source']['sha256'] != source.source['sha256']
                        or metadata['original_basis_identity'] != frame.basis.identity
                        or metadata['original_frame_identity'] != frame.gauge_identity):
                    raise ValueError('state and complete original frame gauge differ')
                state_authority = record['state']
            if not np.array_equal(arrays['populations'], population):
                raise ValueError('saved state and full slot ledger differ')
            state = restore_source(arrays, source_identity=source.source['sha256'])
            coupling = BlochDyadCoupling(background.inputs, frame.basis, block_size=8)
            hamiltonian = MaterialHamiltonian(background, coupling.factors, key[0], block_size=8)
            result = moments_module.full_column_moments(hamiltonian, frame.coefficients,
                                                        state.density_coefficients, row_block=args.row_block)
            report.update(spin=key[0], source_ik=key[1], source_binding=rep.source_binding,
                          state_authority=state_authority, basis_identity=frame.basis.identity,
                          frame_identity=frame.gauge_identity, moments=result,
                          original_source_weight=grid.reference.select(key[1], key[0])[0].weight,
                          original_source_weight_applied=False, full_material_total=False)
    report.update(status='complete', elapsed_seconds=time.monotonic()-started)
    raw = (json.dumps(pack(report), indent=2, allow_nan=False)+'\n').encode()
    if len(raw) > args.output_max_bytes:
        raise ValueError('output exceeds explicit byte capacity')
    args.output.mkdir(parents=True, exist_ok=False)
    with (args.output/'result.json').open('xb') as stream:
        stream.write(raw)
    print(json.dumps(dict(status='complete', output=binding(args.output/'result.json'))))


if __name__ == '__main__':
    main()
