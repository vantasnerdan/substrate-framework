#!/usr/bin/env python3
"""New physical joint mean: initialize, full initial reduction, bounded refinement.

Example (after the independent oracle passes; allocate CPU6 explicitly):
  python run_joint_streamed_material.py --mode initialize --output /var/tmp/joint-init \
    --oracle /var/tmp/p254-joint-streamed-oracle-001.json --cpu 6 \
    --max-actions 100000 --max-seconds 7200 --max-output-bytes 100000000
  python run_joint_streamed_material.py --mode reduce --output /var/tmp/joint-reduce \
    --oracle /var/tmp/p254-joint-streamed-oracle-001.json --cpu 6 \
    --max-actions 100000 --max-seconds 7200 --max-output-bytes 100000000
  python run_joint_streamed_material.py --mode refine --output /var/tmp/joint-refine \
    --oracle /var/tmp/p254-joint-streamed-oracle-001.json --cpu 6 \
    --max-actions 10000000 --max-seconds 86400 --max-output-bytes 45000000000 \
    --interval .001 --steps 1 --max-steps 2 --max-iterations 30 --tolerance 1e-11

These are bounded execution requests, not predictions of completion/accuracy.
Serial slot streaming uses one explicitly allocated CPU and BLAS1. Whole local-q
response and propagated quantum covariance remain separate consumers.
"""
import argparse
from contextlib import contextmanager
from hashlib import sha256
import json
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0025'), str(HERE.parent/'0024')]
for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[variable] = '1'
import numpy as np
from run_streamed_material_evolution import restore_inputs, slot_identity
from run_clamped_array_material import authenticated_arrays
from run_full_field_moments_material import unpack
from material_run_contract import bound_json
from substrate_framework.material_bloch_admission import restore_accepted
from substrate_framework.prepared_reference import PreparedReference
from substrate_framework.material_physical_operators import physical_hamiltonian, physical_coupling, PhysicalFieldVertices
from substrate_framework.material_joint_field_scalar import JointFieldCoefficients
from substrate_framework.physical_field_environment import UniformMaxwellField
from substrate_framework.streamed_material_evolution import real_physical
from substrate_framework.streamed_material_state import file_binding, authenticate
from substrate_framework.streamed_joint_state import DiskJointMean, JointBudget
from substrate_framework.streamed_joint_mean import StreamedJointMean


def implementation_bindings():
    paths = [Path(__file__), HERE/'run_joint_streamed_oracle.py', HERE/'joint-streamed-registration.json']
    paths += [ROOT/'src/substrate_framework'/name for name in (
        'streamed_joint_state.py', 'streamed_joint_reduction.py', 'streamed_joint_mean.py',
        'material_joint_field_scalar.py', 'material_physical_operators.py',
        'streamed_material_state.py', 'streamed_material_evolution.py')]
    return [file_binding(p) for p in paths]


def load_settled():
    registration = json.loads((HERE/'joint-streamed-registration.json').read_text())
    for name, digest in registration['sources'].items():
        if file_binding(ROOT/name)['sha256'] != digest:
            raise ValueError('settled source changed: '+name)
    intake = json.loads((HERE/'full-state-initialization-intake-001.json').read_text())
    original_binding = intake['metadata'][1]['snapshot']
    original_model = bound_json(intake['metadata'][0]['snapshot'])
    spec = original_model['specification']
    restored = restore_inputs(spec['aggregate'], spec['manifest'])
    joint_binding = file_binding(HERE/'joint-field-aggregate-001/result.json')
    joint = unpack(bound_json(joint_binding))
    ionic_binding = joint['ionic_preparation']
    ionic = bound_json(ionic_binding)
    ions = bound_json(spec['initial_ions'])
    if (joint['status'] != 'complete' or joint['periodic_cells'] != 64
            or joint['full_reference'] != spec['aggregate']
            or ionic['status'] != 'complete' or not ionic['covariance_published']
            or ionic['receipt']['metadata']['reference'] != spec['aggregate']
            or joint['field_preparation_status'] != 'conditional-preparation-defined'):
        raise ValueError('settled joint/reference/ionic preparations differ')
    preparation = joint['field_preparation']
    for key in ('coordinate_mean', 'momentum_mean'):
        preparation[key] = real_physical(np.asarray(preparation[key]), key,
                                         preparation['reality_relative_tolerance'])
    ionic_arrays = authenticated_arrays(ionic['receipt'])
    if not np.array_equal(ionic_arrays['masses_ry'], joint['original_atoms']['masses_ry']):
        raise ValueError('ionic original masses differ')
    field_data = joint['maxwell']
    field = UniformMaxwellField(field_data['volume_bohr3'], np.asarray(field_data['polarizations']), field_data['cell_authority'])
    if field.mass != field_data['mass']:
        raise ValueError('field normalization changed')
    c = joint['coefficients']
    coeff = JointFieldCoefficients(c['delta_curvature'], c['contact'], c['mixed'], joint_binding['sha256'])
    new_spec = dict(aggregate=spec['aggregate'], manifest=spec['manifest'], preparation=spec['preparation'],
        initial_ions=spec['initial_ions'], joint=joint_binding, ionic=ionic_binding,
        original_electronic_boundary=original_binding,
        scalar='physical-H0/g plus field-only-L; common-C plus fixed-DeltaC-once; canonical ions/global Maxwell')
    identity = sha256(json.dumps(new_spec, sort_keys=True).encode()).hexdigest()
    return restored, original_binding, ions, joint, coeff, field, new_spec, identity


def build_model(restored, ions, joint, coefficients, field, identity, args, budget):
    arrays, receipt, manifest, saved, original, grid, ordered, pairs = restored
    background, _ = restore_accepted(original)
    if not np.array_equal(arrays['saved_density_bohr3'], background.rho_updown):
        raise ValueError('saved common reference differs from actual H0/K density')
    ndof = 3*grid.reference.nat
    real = lambda value, name: real_physical(value, name, args.reality_tolerance)
    support = np.asarray(ions.get('support_force_ry_bohr', np.zeros(ndof)))
    r = PreparedReference(real(arrays['saved_density_bohr3'], 'nbar'),
        real(arrays['ionic_force_ry_bohr'], 'ionic force').reshape(ndof),
        np.zeros_like(arrays['saved_density_bohr3'], dtype=float), np.zeros(ndof),
        real(arrays['common_q_grid_ry_bohr2'][0,0,0], 'common Gamma C').reshape(ndof, ndof),
        np.asarray(joint['original_atoms']['masses_ry']), support, background.apply_kernel, 64,
        receipt['recipe_sha256'], manifest['ionic_terminal']['sha256'], saved.coverage_identity,
        ions.get('support_authority', 'explicit zero support'), ions['coordinate_identity'],
        manifest['preparation_terminal']['sha256'])
    slots = {slot_identity(s): s for s in ordered}
    @contextmanager
    def factory(record):
        slot = slots[record.identity]
        original_frame = grid.load_frame(slot)
        if (original_frame.basis.identity != record.basis_identity
                or original_frame.gauge_identity != record.gauge_identity
                or original_frame.coefficients.shape != record.shape
                or slot.source_binding['sha256'] != record.source_identity):
            raise ValueError('original full raw basis/gauge/rows/columns changed')
        basis = original_frame.basis
        del original_frame
        h0 = physical_hamiltonian(background, basis, slot.ispin, args.block_size)
        coupling = physical_coupling(background.inputs, basis, h0.factors, args.block_size)
        yield coupling, h0, PhysicalFieldVertices(h0)
    measure = abs(np.linalg.det(grid.reference.cell_bohr))/int(np.prod(grid.reference.fft_smooth))
    return StreamedJointMean(r, factory, coefficients, field, joint['original_atoms']['valence_charges_e'],
        model_identity=identity, expected_slots=tuple(slot_identity(s) for s in ordered),
        grid_measure=measure, budget=budget, reality_tolerance=args.reality_tolerance)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--mode', choices=('initialize', 'reduce', 'evolve', 'refine'), required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--oracle', type=Path, required=True)
    p.add_argument('--checkpoint', type=Path)
    p.add_argument('--checkpoint-sha256')
    p.add_argument('--cpu', type=int, choices=range(6,11), required=True)
    p.add_argument('--max-actions', type=int, required=True)
    p.add_argument('--max-seconds', type=float, required=True)
    p.add_argument('--max-output-bytes', type=int, required=True)
    p.add_argument('--interval', type=float, default=.001)
    p.add_argument('--steps', type=int, default=1)
    p.add_argument('--max-steps', type=int, default=2)
    p.add_argument('--max-iterations', type=int, default=30)
    p.add_argument('--tolerance', type=float, default=1e-11)
    p.add_argument('--reality-tolerance', type=float, default=1e-10)
    p.add_argument('--block-size', type=int, default=8)
    args = p.parse_args()
    os.sched_setaffinity(0, {args.cpu})
    budget = JointBudget(max_actions=args.max_actions, max_seconds=args.max_seconds, max_bytes=args.max_output_bytes)
    budget.reserve(65536)  # A terminal partial/failure receipt remains writable after a cap.
    args.output.mkdir(parents=True, exist_ok=False)
    result = dict(schema='P254-joint-streamed-material-v1', status='partial', mode=args.mode,
        scope='homogeneous coupled physical mean; whole local-q dynamics and quantum response separate',
        requested_interval=args.interval if args.mode in ('evolve','refine') else 0., cpu=args.cpu, blas_threads=1)
    try:
        bindings = implementation_bindings()
        oracle = bound_json(file_binding(args.oracle))
        if oracle['status'] != 'passed' or oracle['schema'] != 'P254-joint-streamed-oracle-v1':
            raise ValueError('executed joint disk/dynamics oracle required')
        expected = {b['path']: b['sha256'] for b in oracle['source_bindings']}
        for b in bindings:
            # The finite oracle exercises these numerical modules, not this
            # saved-material CLI. Bind the current CLI in model.json and prove
            # its data path through the actual material run instead.
            if Path(b['path']) == Path(__file__).resolve():
                continue
            if expected.get(b['path']) != b['sha256']:
                raise ValueError('oracle does not bind current source: '+b['path'])
        result['numerical_oracle'] = file_binding(args.oracle)
        result['verification_scope'] = 'unchanged numerical modules plus this actual saved-material execution'
        restored, original_binding, ions, joint, coefficients, field, spec, identity = load_settled()
        ordered = restored[6]
        if args.checkpoint:
            b = file_binding(args.checkpoint)
            if b['sha256'] != args.checkpoint_sha256:
                raise ValueError('checkpoint explicit digest differs')
            state = DiskJointMean(args.checkpoint, expected_identity=identity)
        else:
            state = DiskJointMean.wrap_original(original_binding, args.output/'initial', identity=identity,
                coordinate=joint['field_preparation']['coordinate_mean'], momentum=joint['field_preparation']['momentum_mean'],
                expected_slots=tuple(slot_identity(s) for s in ordered),
                evidence=dict(new_scalar=spec, field_preparation=spec['joint']), budget=budget)
            if not np.array_equal(state.u, ions['u_bohr']) or not np.array_equal(state.p, ions['p_canonical']):
                raise ValueError('original ionic means differ from settled preparation')
        budget.json(args.output/'model.json', dict(identity=identity, specification=spec, implementation=bindings))
        result['initial_boundary'] = file_binding(state.path)
        if args.mode != 'initialize':
            model = build_model(restored, ions, joint, coefficients, field, identity, args, budget)
            initial = model.reduce(state)
            record = initial.save(args.output/'initial-reduction', budget, metadata=dict(kind='full-initial-reduction', time_ry_inverse=state.time))
            budget.json(args.output/'initial-reduction.json', record)
            result['initial_energy_total_ry'] = initial.energy
            result['initial_number_per_cell'] = initial.number
            result['initial_force_aggregate_difference_max'] = float(np.max(np.abs(initial.force-np.asarray(joint['initial_join']['force_per_cell_ry_bohr']).reshape(-1))))
            result['initial_field_rhs_aggregate_difference_max'] = float(np.max(np.abs(initial.field_force-joint['initial_join']['field_momentum_rhs'])))
            if args.mode in ('evolve', 'refine'):
                branches = (1,2) if args.mode == 'refine' else (args.steps,)
                results = []
                for steps in branches:
                    final, history = model.evolve(state, args.output/('steps-%d' % steps), physical_interval=args.interval,
                        steps=steps, max_steps=args.max_steps, max_iterations=args.max_iterations,
                        tolerance=args.tolerance, refinement_label='steps-'+str(steps))
                    results.append(final)
                result['final_boundaries'] = [file_binding(s.path) for s in results]
                if len(results) == 2:
                    a, b = results
                    frame_error = max(float(np.linalg.norm(a.arrays(s)[0]-b.arrays(t)[0]))
                                      for s,t in zip(a.slots,b.slots))
                    result['refinement'] = dict(max_frame_frobenius=frame_error,
                        coordinate_max=max(float(np.max(np.abs(getattr(a,k)-getattr(b,k)))) for k in ('u','p','coordinate','momentum')),
                        accuracy_certified=False, interval_not_an_accuracy_license=True)
        result['status'] = 'complete'
    except Exception as error:
        result['error'] = type(error).__name__+': '+str(error)
    result['completed_boundaries'] = [str(p.resolve()) for p in args.output.glob('**/boundary.json')]
    result['budget'] = budget.snapshot()
    encoded = (json.dumps(result, indent=2, allow_nan=False)+'\n').encode()
    if len(encoded) > 65536:
        raise RuntimeError('reserved terminal report capacity exceeded')
    with (args.output/'result.json').open('xb') as stream:
        stream.write(encoded)
    if result['status'] != 'complete':
        raise SystemExit(2)


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()
