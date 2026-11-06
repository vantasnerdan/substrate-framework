#!/usr/bin/env python3
"""Additional mixed source contractions and complete conditional q0 scalar join.

Run the independent oracle first. No DFT, propagated trajectories or noise.
Source mode performs genuinely NEW full-column overlaps, not receipt reuse.
Aggregate mode consumes all20 mixed sources and complete optical receipts.
"""
import argparse
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'),str(HERE.parent/'0025'),str(HERE.parent/'0024')]

import numpy as np
import run_full_field_moments_material as full
from material_run_contract import binding, bound_json, recover_block, runtime_context
from material_reference_output import load_arrays
from substrate_framework import material_mixed_field_ionic as mixed_module
from substrate_framework import material_joint_field_scalar as scalar_module
from substrate_framework.material_modes import AMU_RY
from substrate_framework.physical_field_environment import UniformMaxwellField


def implementation_bindings():
    names = ('material_joint_field_scalar.py','material_mixed_field_ionic.py',
             'prepared_factor_dynamics.py','prepared_reference.py','physical_field_environment.py',
             'qe_material_field_derivatives.py','material_physical_field_moments.py',
             'material_reference_transport.py','material_modes.py')
    paths = [Path(__file__),HERE/'run_joint_field_scalar_oracle.py',
             HERE/'run_full_field_moments_material.py',
             *[ROOT/'src/substrate_framework'/name for name in names],
             HERE.parent/'0025/material_run_contract.py',
             HERE.parent/'0025/material_reference_output.py']
    return [binding(path) for path in paths]


def authenticate_sources():
    record = binding(HERE/'joint-field-scalar-registration.json')
    registration = bound_json(record)
    for parent in registration['sources']:
        source_registration = bound_json(parent)
        for source in source_registration.get('source_bindings',source_registration.get('original_bytes',[])):
            actual = binding(source['path'])
            if (actual['sha256'] != source['sha256']
                    or int(actual['bytes']) != int(source.get('bytes',source.get('size_bytes')))):
                raise ValueError('original registered source bytes changed: '+source['path'])
        for source in source_registration.get('accepted_material',{}).get('inputs',{}).get('upfs',{}).values():
            actual = binding(source['path'])
            if actual['sha256'] != source['sha256'] or int(actual['bytes']) != int(source['bytes']):
                raise ValueError('original UPF bytes changed: '+source['path'])
    return record, registration


def recover_receipt(receipt):
    archive = Path(receipt['archive']['path'])
    if recover_block(archive.parent,archive.stem,receipt['recipe_sha256'],receipt['metadata']) != receipt:
        raise ValueError('completed array receipt differs on disk')
    return load_arrays(receipt)


def load_source(key, manifest, saved, original, grid, background):
    selected = [r for r in manifest['sources'] if (r['spin'],r['source_ik']) == key]
    if len(selected) != 1:
        raise ValueError('select one exact original spin/source pair')
    record = selected[0]
    population, pairs = full.source_population_indices(grid,saved.population_ledger,*key)
    rep = next(s for _,s in pairs if s.operation.index == 1 and not np.any(s.reciprocal_shift))
    source = full.load_wavefunctions(rep.source_binding)
    frame = full.transport_bloch_frame(grid.reference,rep,source)
    if record['reuse_selected_gamma']:
        if key != (1,1) or record['state'] is not None:
            raise ValueError('selected Gamma identity differs')
        arrays = {name:saved.selected_gamma[name] for name in
                  ('gram','populations','density_coefficients','density_eigenvalues')}
        state_authority = manifest['preparation_terminal']
    else:
        arrays, receipt = full.receipt_arrays(record['state'])
        metadata = receipt['metadata']
        if (metadata['source']['sha256'] != source.source['sha256']
                or metadata['original_basis_identity'] != frame.basis.identity
                or metadata['original_frame_identity'] != frame.gauge_identity):
            raise ValueError('original complete frame and fixed D authority differ')
        state_authority = record['state']
    if not np.array_equal(arrays['populations'],population):
        raise ValueError('source populations disagree with complete slot ledger')
    state = full.restore_source(arrays,source_identity=source.source['sha256'])
    coupling = full.BlochDyadCoupling(background.inputs,frame.basis,block_size=8)
    hamiltonian = full.MaterialHamiltonian(background,coupling.factors,key[0],block_size=8)
    return hamiltonian,frame,state,rep,state_authority


def optical_moments(path, manifest_binding, reference_binding, original, grid):
    record = binding(path)
    aggregate = full.unpack(bound_json(record))
    expected_implementation = full.implementation_bindings()
    expected_registration = binding(HERE/'full-field-moments-registration.json')
    if (aggregate.get('status') != 'complete' or aggregate.get('mode') != 'aggregate'
            or aggregate['implementation_sources'] != expected_implementation
            or aggregate['registration'] != expected_registration
            or aggregate['manifest'] != manifest_binding or aggregate['full_reference'] != reference_binding
            or aggregate['aggregate']['geometry_policy'] != 'declared-geometric'
            or len(aggregate['source_checkpoints']) != 20):
        raise ValueError('complete matching original optical aggregate required')
    oracle = bound_json(aggregate['oracle'])
    if (oracle.get('status') != 'passed' or oracle['implementation_sources'] != expected_implementation
            or oracle['registration'] != expected_registration):
        raise ValueError('optical aggregate independent oracle authority differs')
    sources = {}
    for source_record in aggregate['source_checkpoints']:
        source = full.unpack(bound_json(source_record))
        key = (source['spin'],source['source_ik'])
        if (key in sources or source.get('status') != 'complete' or source.get('mode') != 'source'
                or source['implementation_sources'] != expected_implementation
                or source['registration'] != expected_registration or source['manifest'] != manifest_binding
                or source['source_binding'] != full.wave_bindings(original)[key]):
            raise ValueError('optical source checkpoint coverage/authority differs')
        sources[key] = source['moments']
    recomputed = full.moments_module.aggregate_slots(grid,sources,geometry_policy='declared-geometric')
    for name,value in recomputed['moments'].items():
        if not np.array_equal(value,aggregate['aggregate']['moments'][name]):
            raise ValueError('optical scalar aggregate disagrees with complete source receipts: '+name)
    return recomputed['moments'],record


def aggregate_join(args, report, manifest_binding, reference_binding, original, grid, implementation):
    if len(args.mixed_source) != 20 or args.geometry_policy != 'declared-geometric':
        raise ValueError('exactly20 mixed receipts and explicit geometry policy required')
    if args.optical_aggregate is None or args.ionic_preparation is None:
        raise ValueError('complete optical aggregate and actual ionic preparation required')
    sources, records = {}, []
    for path in args.mixed_source:
        record = binding(path)
        item = full.unpack(bound_json(record))
        key = (item['spin'],item['source_ik'])
        if (key in sources or item.get('mode') != 'source' or item.get('status') != 'complete'
                or item['registration'] != report['registration']
                or item['implementation_sources'] != implementation
                or item['manifest'] != manifest_binding or item['full_reference'] != reference_binding
                or item['source_binding'] != full.wave_bindings(original)[key]):
            raise ValueError('mixed source checkpoint coverage/authority differs')
        sources[key] = item['mixed']
        records.append(record)
    mixed = mixed_module.aggregate_mixed(grid,sources,geometry_policy=args.geometry_policy)
    optical,optical_record = optical_moments(args.optical_aggregate,manifest_binding,reference_binding,original,grid)
    reference_receipt = bound_json(reference_binding)
    arrays = recover_receipt(reference_receipt)
    metadata = reference_receipt['metadata']
    if (metadata['complete_electronic_coverage'] is not True or len(metadata['slots']) != 128
            or len(metadata['source_records']) != 20):
        raise ValueError('complete common reference required')
    background,_ = full.restore_accepted(original)
    names = tuple(grid.reference.atom_names)
    charges = np.asarray([background.inputs.tables[name].upf.z_valence for name in names])
    masses_by_name = {s.name:s.mass_amu for s in grid.reference.species}
    masses_amu = np.asarray([masses_by_name[name] for name in names])
    if (not np.array_equal(masses_amu,arrays['atomic_masses_amu'])
            or not np.array_equal(grid.reference.cell_bohr,arrays['cell_bohr'])
            or not np.array_equal(grid.reference.positions_bohr,arrays['positions_bohr'])):
        raise ValueError('original atom masses/cell/positions differ from common scalar')
    masses = np.repeat(masses_amu*AMU_RY,3)
    mesh = tuple(grid.reference.monkhorst_pack['nk'+str(i)] for i in (1,2,3))
    if mesh != (4,4,4):
        raise ValueError('this actual consumer requires original BvK4^3 mesh')
    ncells = 64
    field = UniformMaxwellField(float(abs(np.linalg.det(arrays['cell_bohr'])))*ncells,
                                np.eye(3),reference_binding['sha256']+'; BvK4^3 not cavity')
    coefficients = scalar_module.JointFieldCoefficients(
        scalar_module.onsite_matrix(optical['delta_nonlocal_onsite_ry_bohr2']),
        optical['contact'],mixed['mixed_ionic_field'].reshape(-1,3),report['registration']['sha256'])
    ionic_record = binding(args.ionic_preparation)
    ionic = bound_json(ionic_record)
    if (ionic.get('status') != 'complete' or ionic.get('covariance_published') is not True
            or ionic['receipt']['metadata']['reference'] != reference_binding):
        raise ValueError('actual complete ionic preparation for this reference required')
    ionic_arrays = recover_receipt(ionic['receipt'])
    means = bound_json(ionic['receipt']['metadata']['initial_means'])
    ndof = 3*len(charges)
    u,p = (np.asarray(means[name]) for name in ('u_bohr','p_canonical'))
    if u.shape != (ndof,) or p.shape != (ndof,) or np.any(u) or np.any(p):
        raise ValueError('registered sudden-release hypothesis requires actual zero ionic means')
    if not np.array_equal(ionic_arrays['masses_ry'],masses):
        raise ValueError('actual ionic covariance and original Ry masses differ')
    correction = scalar_module.coordinate_correction(coefficients,u,p,np.zeros(3),masses,charges,
        optical['delta_nonlocal_positive_gradient_ry_bohr'].reshape(-1),optical['gradient'])
    q0 = arrays['common_q_grid_ry_bohr2'][0,0,0].reshape(ndof,ndof)
    if not np.allclose(q0,q0.T,rtol=1e-10,atol=1e-10):
        raise ValueError('original common q0 curvature is not symmetric; no repair')
    report.update(mixed_source_checkpoints=records,optical_aggregate=optical_record,
        ionic_preparation=ionic_record,mixed_aggregate=mixed,
        original_atoms=dict(names=names,indices_one_based=list(grid.reference.atom_indices),
                            valence_charges_e=charges,masses_amu=masses_amu,masses_ry=masses),
        maxwell=field.normalization(),periodic_cells=ncells,
        coefficients=dict(contact=coefficients.contact,mixed=coefficients.mixed,
            delta_curvature=coefficients.delta_curvature,
            common_q0_corrected=q0+coefficients.delta_curvature,
            spatial_correction='same onsite DeltaC added at R=0, equivalently at every existing matter q; no ASR repair'),
        initial_join=dict(force_per_cell_ry_bohr=arrays['force_at_origin_ry_bohr']+correction['force_correction'],
            original_force_per_cell_ry_bohr=arrays['force_at_origin_ry_bohr'],
            electronic_current_total_hartree=ncells*correction['electronic_current_hartree'],
            ionic_current_total_hartree=ncells*correction['ionic_current_hartree'],
            field_momentum_rhs=-ncells*correction['field_gradient'],
            field_coordinate_rhs=np.zeros(3),ionic_velocity=correction['velocity'],
            unanchored_energy_correction_total_ry=ncells*optical['delta_nonlocal_energy_ry'],
            energy_scope='mean scalar values below; product covariance correction reported separately',
            ionic_product_curvature_correction_total_ry=ncells*np.trace(
                coefficients.delta_curvature @ ionic_arrays['symmetrized_uu'])/2,
            absolute_total_energy_ry=None,
            absolute_energy_missing='actual full-slot H0/K evaluation; original-source representative energy is not transported H0',
            scalar_join='E_total=64*(existing full-slot PreparedReference scalar+unanchored correction)+Maxwell; runtime JointPreparedFieldDynamics evaluates absolute scalar and every RHS together',
            current_density_total_hartree_per_bohr3=ncells*correction['electronic_current_hartree']/field.volume_bohr3),
        ionic_preparation_changed=False,product_cross_covariance='initially zero by explicit tensor-product hypothesis, not coupled equilibrium')
    try:
        report['field_preparation'] = scalar_module.conditional_field_gaussian(field,coefficients.contact,
            masses,charges,ncells,state_authority=report['registration']['sha256'])
        report['field_preparation_status'] = 'conditional-preparation-defined'
        report['initial_join']['product_energy_correction_total_ry'] = (
            report['initial_join']['unanchored_energy_correction_total_ry']
            +report['initial_join']['ionic_product_curvature_correction_total_ry']
            +report['field_preparation']['quadratic_zero_point_energy_ry'])
    except scalar_module.MissingFieldPreparation as error:
        report['field_preparation'] = None
        report['field_preparation_status'] = 'missing-preparation-input'
        report['field_preparation_missing'] = str(error)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('estimate','source','aggregate'),required=True)
    parser.add_argument('--full-reference',type=Path,default=HERE/'full-common-reference-004/full-common-reference.json')
    parser.add_argument('--oracle',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--output-max-bytes',type=int,default=16777216)
    parser.add_argument('--spin',type=int)
    parser.add_argument('--source-ik',type=int)
    parser.add_argument('--row-block',type=int,default=256)
    parser.add_argument('--mixed-source',type=Path,action='append',default=[])
    parser.add_argument('--optical-aggregate',type=Path)
    parser.add_argument('--ionic-preparation',type=Path)
    parser.add_argument('--geometry-policy',choices=('declared-geometric',))
    args = parser.parse_args()
    if args.row_block < 1 or args.output_max_bytes < 65536:
        raise ValueError('positive row block and explicit sufficient output capacity required')
    started = time.monotonic()
    registration_binding,registration = authenticate_sources()
    implementation = implementation_bindings()
    oracle_binding = binding(args.oracle)
    oracle = bound_json(oracle_binding)
    if (oracle.get('status') != 'passed' or oracle.get('registration') != registration_binding
            or oracle.get('implementation_sources') != implementation):
        raise ValueError('passed independent oracle bound to current implementation required')
    full.preload_parsers()
    reference_binding,manifest_binding,manifest,saved,original,grid = full.load_inputs(args.full_reference)
    report = dict(schema='P254-joint-field-material-v1',mode=args.mode,
        registration=registration_binding,implementation_sources=implementation,oracle=oracle_binding,
        full_reference=reference_binding,manifest=manifest_binding,runtime_context=runtime_context(),
        scope=registration['scope'],scalar=registration['scalar'],h0_transport_supported=False,
        finite_q_photons=False,propagated_noise=False,propagated_trajectory=False)
    if args.mode == 'aggregate':
        aggregate_join(args,report,manifest_binding,reference_binding,original,grid,implementation)
    else:
        background,background_receipt = full.restore_accepted(original)
        report['background_receipt'] = background_receipt
        if args.mode == 'estimate':
            columns = 16*sum(sum(name == species.name for name in grid.reference.atom_names)
                             *background.inputs.tables[species.name].upf.nh
                             for species in grid.reference.species)
            counts = []
            for source in manifest['sources']:
                point,energies,_ = grid.reference.select(source['source_ik'],source['spin'])
                rows,n = int(point.npw),int(energies.size)
                counts.append(dict(spin=source['spin'],source_ik=source['source_ik'],pw_rows=rows,
                    raw_columns=n,frame_bytes=16*rows*n,wide_projector_columns=columns,
                    complex_overlap_multiply_adds=rows*n*columns,overlap_workspace_upper_bytes=16*n*columns))
            report['estimate'] = dict(sources=counts,
                total_complex_overlap_multiply_adds=sum(v['complex_overlap_multiply_adds'] for v in counts),
                row_workspace='O(row_block*(16*max_nh+raw_columns)); original background/species beta bank additional',
                kinetic_density_pass=False,local_fft_count=0,pw_square_matrices=0,
                output_complex_scalars_per_source=9*grid.reference.nat,
                wall_time_prediction=None,existing_moment_receipts_supply_overlaps=False)
        else:
            key = (args.spin,args.source_ik)
            h,frame,state,rep,state_authority = load_source(key,manifest,saved,original,grid,background)
            result = mixed_module.full_column_mixed(h,frame.coefficients,state.density_coefficients,row_block=args.row_block)
            report.update(spin=key[0],source_ik=key[1],source_binding=rep.source_binding,
                state_authority=state_authority,basis_identity=frame.basis.identity,frame_identity=frame.gauge_identity,
                mixed=result,original_source_weight_applied=False,full_material_total=False)
    report.update(status='complete',elapsed_seconds=time.monotonic()-started)
    raw = (json.dumps(full.pack(report),indent=2,allow_nan=False)+'\n').encode()
    if len(raw) > args.output_max_bytes:
        raise ValueError('output exceeds explicit byte capacity')
    args.output.mkdir(parents=True,exist_ok=False)
    with (args.output/'result.json').open('xb') as stream:
        stream.write(raw)
    print(json.dumps(dict(status='complete',output=binding(args.output/'result.json'))))


if __name__ == '__main__':
    main()
