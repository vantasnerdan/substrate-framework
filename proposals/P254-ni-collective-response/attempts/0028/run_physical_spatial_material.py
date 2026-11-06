#!/usr/bin/env python3
"""Physical spatial/field cutover witnesses, not a filled response or evolution.

Actual full-PW saved/complement operator probes; original full frames and D stay
source-bound. Reuse completed scalar coefficients without rerunning20 sources.
Large witness arrays must be written outside Git. Root owns all execution.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE.parent/'0025'), str(HERE.parent/'0024')]

import numpy as np
import run_joint_field_material as joint
import run_full_field_moments_material as full
from material_run_contract import bound_json, runtime_context
from substrate_framework.material_physical_operators import (
    physical_hamiltonian, physical_coupling, physical_finiteq,
    PhysicalFieldVertices, PhysicalFrozenAction, spatial_coefficient_handoff)
from substrate_framework.material_mixed_field_ionic import UPFPreparedVertices
from substrate_framework.material_joint_field_scalar import JointFieldCoefficients, onsite_matrix
from substrate_framework.qe_bloch_edge import BlochEdge
from substrate_framework.qe_bloch_dyads import BlochDyadCoupling


def binding(path):
    p = Path(path).resolve()
    h, size = hashlib.sha256(), 0
    with p.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
            size += len(block)
    return dict(path=str(p), bytes=str(size), sha256=h.hexdigest())


def implementation_bindings():
    names = ('material_physical_operators.py', 'qe_material_hamiltonian.py',
             'qe_material_field_derivatives.py', 'material_mixed_field_ionic.py',
             'qe_bloch_dyads.py', 'qe_bloch_ionic_action.py', 'qe_bloch_finiteq_ionic.py',
             'material_frozen_action.py', 'material_joint_field_scalar.py')
    return [binding(p) for p in (Path(__file__), HERE/'run_physical_spatial_oracle.py',
            *[ROOT/'src/substrate_framework'/n for n in names])]


def authenticate_registration():
    record = binding(HERE/'physical-spatial-registration.json')
    registration = bound_json(record)
    for source in registration['source_bindings']:
        if binding(source['path']) != source:
            raise ValueError('registered original bytes changed: '+source['path'])
    return record, registration


def relative(actual, expected):
    a, b = np.asarray(actual), np.asarray(expected)
    return float(np.linalg.norm(a-b)/max(1., float(np.linalg.norm(b))))


def compare(records, name, actual, expected, tolerance=3e-10):
    residual = relative(actual, expected)
    records[name] = residual
    if not np.isfinite(residual) or residual > tolerance:
        raise AssertionError(name+': relative residual '+str(residual))


def enforce_runtime():
    affinity = set(os.sched_getaffinity(0))
    if not affinity or not affinity <= set(range(6, 11)):
        raise ValueError('Root must authorize execution on CPUs6-10 only')
    for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        if os.environ.get(name) != '1':
            raise ValueError(name+'=1 required')


def fresh_output(path):
    p = Path(path).resolve()
    if p == ROOT or ROOT in p.parents:
        raise ValueError('large operator witness arrays must remain outside Git/repository')
    p.mkdir(parents=True, exist_ok=False)
    return p


def probe_columns(frame, state, columns):
    """Original saved columns plus unnormalized raw-span complement witnesses.

    Projection uses all raw columns and the ORIGINAL Gram; it does not repair
    the frame or alter its state. D[:2m,:2m] is a declared operator coefficient
    witness on these probes, never a population assignment to the complement.
    """
    psi = frame.coefficients
    rng = np.random.default_rng(2542800+frame.basis.slot.source_ik)
    trial = rng.normal(size=(len(psi), columns))+1j*rng.normal(size=(len(psi), columns))
    trial /= np.sqrt(len(psi))
    overlap = psi.conj().T @ trial
    complement = trial-psi @ np.linalg.solve(state.gram, overlap)
    probes = np.column_stack((psi[:, :columns], complement))
    d = state.density_coefficients[:2*columns, :2*columns]
    return probes, d, dict(complement_raw_span_overlap_norm=float(np.linalg.norm(psi.conj().T @ complement)),
        complement_norms=np.linalg.norm(complement, axis=0).tolist(), saved_columns=list(range(columns)),
        probe_coefficient='unchanged principal2m block of original D; operator witness ONLY, not filled state')


def same_basis_witness(legacy, frame, state, columns):
    inputs, basis = legacy.background.inputs, frame.basis
    h = physical_hamiltonian(legacy.background, basis, legacy.ispin, block_size=2)
    coupling = physical_coupling(inputs, basis, h.factors, block_size=2)
    old = BlochDyadCoupling(inputs, basis, factors=legacy.factors, block_size=2)
    psi, d, probe_report = probe_columns(frame, state, columns)
    vertices = PhysicalFieldVertices(h, projector_block_size=64)
    reference = UPFPreparedVertices(legacy)
    u = np.zeros((inputs.reference.nat, 3))
    # Include displaced atoms away from origin; no change to original geometry.
    u[0] = (.13, -.07, .19)
    u[-1] = (-.11, .17, .03)
    a = np.array([.017, -.013, .011])
    zero = reference.evaluate(psi, d, np.zeros_like(u), np.zeros(3))
    shifted = reference.evaluate(psi, d, u, a)
    ionic = reference.evaluate(psi, d, u, np.zeros(3))
    hphys, hold = h.apply(psi), legacy.apply(psi)
    gphys = coupling.ionic_action(u.ravel(), psi)
    gold = old.ionic_action(u.ravel(), psi)
    lphys = vertices.prepare(a).apply(psi)
    density, gradient = coupling.moments(psi, psi, d)
    old_density, old_gradient = old.moments(psi, psi, d)
    moments = vertices.moments(psi, d)
    frozen = PhysicalFrozenAction(h, coupling, direction=u.ravel(), shift=a,
                                  fft_block_size=2, projector_block_size=64)
    records = {}
    compare(records, 'physical_H0_equals_A25_plus_DeltaV0', hphys, hold+zero['action'])
    compare(records, 'physical_gL_equals_A25_plus_independent_vertices',
            gphys+lphys, gold+shifted['action']-zero['action'])
    compare(records, 'physical_g_equals_A25_plus_Delta_g_action',
            gphys, gold+ionic['action']-zero['action'])
    compare(records, 'physical_L_equals_independent_UPF_action',
            lphys, shifted['action']-ionic['action'])
    compare(records, 'positive_gradient_equals_A25_plus_Delta_g', gradient,
            old_gradient+zero['delta_gradient'])
    compare(records, 'density_unchanged', density, old_density)
    compare(records, 'physical_field_gradient', moments['field_gradient'], zero['field_gradient'])
    compare(records, 'physical_current_sign', moments['electronic_current_hartree'], -zero['field_gradient']/2)
    compare(records, 'fused_physical_first_order', frozen.apply(psi), hold+gold+shifted['action'])
    field_only = vertices.evaluate(psi, d, u, a)
    compare(records, 'no_second_zero_field_correction', field_only['action'], lphys)
    compare(records, 'exact_zero_delta_energy', field_only['delta_energy'], 0.)
    compare(records, 'exact_zero_delta_gradient', field_only['delta_gradient'], np.zeros_like(gradient))
    compare(records, 'gradient_moment', u.ravel() @ gradient, np.trace(d @ (psi.conj().T @ gphys)))
    compare(records, 'field_moment', a @ moments['field_gradient'], np.trace(d @ (psi.conj().T @ lphys)))
    double = relative(hphys+zero['action'], hold+zero['action'])
    if double <= 1e-12:
        raise AssertionError('actual zero-field double correction not resolved by these probes')
    arrays = dict(probes=psi, probe_D=d, original_full_D=state.density_coefficients,
        original_full_gram=state.gram, source_miller=basis.miller, k_bohr_inv=basis.k_cart_bohr_inv,
        physical_H0_action=hphys, A25_H0_action=hold, correction_DeltaV0=zero['action'],
        physical_ionic_action=gphys, A25_ionic_action=gold, physical_field_action=lphys,
        A25_plus_UPFPreparedVertices_action=hold+gold+shifted['action'],
        independent_UPF_field_action=shifted['action']-ionic['action'],
        independent_Delta_g_action=ionic['action']-zero['action'],
        positive_gradient=gradient, A25_positive_gradient=old_gradient,
        correction_Delta_g=zero['delta_gradient'], field_gradient=moments['field_gradient'],
        electronic_current_hartree=moments['electronic_current_hartree'], displacement=u, field_shift=a)
    report = dict(basis_identity=basis.identity, frame_identity=frame.gauge_identity,
        original_source_binding=dict(basis.slot.source_binding), original_full_frame_shape=list(frame.coefficients.shape),
        source_identity=state.source_identity, probe=probe_report, residuals=records,
        double_zero_field_correction_exposure=double)
    return h, psi, d, arrays, report


def rectangular_witness(inputs, source, target, source_h, target_h, right, left, d):
    edge = BlochEdge.create(inputs.reference, source.basis, target.basis, (0., 0., .25))
    if tuple(edge.reciprocal_fold) != (0, 0, -1):
        raise ValueError('actual ik2->ik3 edge fold changed')
    coupling = physical_finiteq(inputs, edge, source_h.factors, target_h.factors, block_size=2)
    direction = np.zeros((inputs.reference.nat, 3), complex)
    direction[0] = (.17+.11j, -.09j, .13)
    direction[-1] = (-.21j, .19, -.07+.03j)
    prepared = coupling.prepare(direction)
    forward, adjoint = prepared.apply(right), prepared.adjoint_apply(left)
    density, gradient = coupling.moments(left, right, d, weight=source.basis.slot.weight)
    overlap = left.conj().T @ forward
    records = {}
    compare(records, 'rectangular_hilbert_adjoint', overlap, adjoint.conj().T @ right)
    compare(records, 'complex_original_D_moment', np.vdot(direction.ravel(), gradient),
            source.basis.slot.weight*np.conj(np.vdot(d, overlap)))
    # The q0 limit must agree with the very same factors on the diagonal leg.
    diagonal_edge = BlochEdge.create(inputs.reference, source.basis, source.basis, (0., 0., 0.))
    diagonal = physical_finiteq(inputs, diagonal_edge, source_h.factors, source_h.factors, block_size=2)
    gamma = physical_coupling(inputs, source.basis, source_h.factors, block_size=2)
    compare(records, 'rectangular_q0_equals_same_basis', diagonal.prepare(direction).apply(right),
            gamma.ionic_action(direction.ravel(), right))
    return dict(source_probe=right, target_probe=left, original_D_probe_coefficient=d,
        source_miller=source.basis.miller, target_miller=target.basis.miller,
        source_k=source.basis.k_cart_bohr_inv, target_k=target.basis.k_cart_bohr_inv,
        q_fractional=edge.q_fractional, q_bohr_inv=edge.q_cart_bohr_inv,
        reciprocal_fold=edge.reciprocal_fold, direction=direction,
        physical_forward_action=forward, physical_adjoint_action=adjoint,
        density_bohr3=density, positive_gradient=gradient), dict(residuals=records,
        source_basis_identity=source.basis.identity, target_basis_identity=target.basis.identity,
        source_frame_identity=source.gauge_identity, target_frame_identity=target.gauge_identity,
        source_binding=dict(source.basis.slot.source_binding), target_binding=dict(target.basis.slot.source_binding),
        source_full_frame_shape=list(source.coefficients.shape), target_full_frame_shape=list(target.coefficients.shape),
        scope='all actual endpoint PW rows; selected saved/complement operator probes, NOT populated response')


def coefficient_join(joint_record, reference_record, manifest_record):
    """Initial scalar join against completed full-source receipts; no PW actions."""
    data = full.unpack(bound_json(joint_record))
    if (data.get('status') != 'complete' or data.get('mode') != 'aggregate'
            or data['full_reference'] != reference_record or data['manifest'] != manifest_record
            or len(data['mixed_source_checkpoints']) != 20 or data['periodic_cells'] != 64):
        raise ValueError('settled complete joint-field aggregate required')
    optical = full.unpack(bound_json(data['optical_aggregate']))
    if (optical.get('status') != 'complete' or optical.get('mode') != 'aggregate'
            or optical['full_reference'] != reference_record or optical['manifest'] != manifest_record
            or len(optical['source_checkpoints']) != 20):
        raise ValueError('matching completed optical aggregate required')
    # Authentication only: do not recompute any source or geometric aggregate.
    for record in data['mixed_source_checkpoints']+optical['source_checkpoints']:
        bound_json(record)
    common = joint.recover_receipt(bound_json(reference_record))
    c = data['coefficients']
    coefficients = JointFieldCoefficients(c['delta_curvature'], c['contact'], c['mixed'], joint_record['sha256'])
    ndof = coefficients.mixed.shape[0]
    original_q = common['common_q_grid_ry_bohr2'].reshape(4, 4, 4, ndof, ndof)
    corrected, residual = spatial_coefficient_handoff(original_q, coefficients)
    optical_values = optical['aggregate']['moments']
    initial = data['initial_join']
    records = {}
    compare(records, 'onsite_correction_from_original_optical', coefficients.delta_curvature,
            onsite_matrix(optical_values['delta_nonlocal_onsite_ry_bohr2']))
    compare(records, 'contact_from_original_optical', coefficients.contact, optical_values['contact'])
    compare(records, 'mixed_from_completed_mixed_aggregate', coefficients.mixed,
            data['mixed_aggregate']['mixed_ionic_field'].reshape(ndof, 3))
    compare(records, 'physical_q0_matches_joint', corrected[0, 0, 0], c['common_q0_corrected'])
    compare(records, 'initial_positive_gradient_force_join',
            common['force_at_origin_ry_bohr']-optical_values['delta_nonlocal_positive_gradient_ry_bohr'].reshape(-1),
            initial['force_per_cell_ry_bohr'])
    compare(records, 'initial_current_join', -32*optical_values['gradient'], initial['electronic_current_total_hartree'])
    compare(records, 'initial_Maxwell_join', -64*optical_values['gradient'], initial['field_momentum_rhs'])
    compare(records, 'initial_ionic_current', np.zeros(3), initial['ionic_current_total_hartree'])
    compare(records, 'initial_unanchored_energy_join', 64*optical_values['delta_nonlocal_energy_ry'],
            initial['unanchored_energy_correction_total_ry'])
    # Independent polynomial identity for the two mutually exclusive C routes.
    u = np.linspace(-.017, .023, ndof)
    a = np.array([.011, -.007, .013])
    old_route = (u @ original_q[0, 0, 0] @ u+u @ coefficients.delta_curvature @ u)/2
    new_route = u @ corrected[0, 0, 0] @ u/2
    compare(records, 'once_only_C_handoff', new_route, old_route)
    fixed = a @ coefficients.contact @ a/2+u @ coefficients.mixed @ a
    compare(records, 'contact_mixed_retained_once',
            a @ residual.contact @ a/2+u @ residual.mixed @ a, fixed)
    compare(records, 'residual_DeltaC_exact_zero', residual.delta_curvature, np.zeros((ndof, ndof)))
    arrays = dict(common_q_original=original_q, physical_q_corrected=corrected,
        delta_curvature=coefficients.delta_curvature, contact=coefficients.contact,
        mixed=coefficients.mixed, residual_delta_curvature=residual.delta_curvature,
        physical_initial_force=initial['force_per_cell_ry_bohr'],
        physical_initial_electronic_current=initial['electronic_current_total_hartree'],
        physical_initial_field_momentum_rhs=initial['field_momentum_rhs'])
    report = dict(joint_aggregate=joint_record, optical_aggregate=data['optical_aggregate'],
        full_reference=reference_record, manifest=manifest_record, residuals=records,
        geometry_policy='inherited declared-geometric; no H0 transport authorized',
        corrected_C_imaginary_max=float(np.max(np.abs(corrected.imag))),
        prepared_reference_route='Keep original strictly-real PreparedReference.curvature and full fixed joint DeltaC; replace only H0/g and use PhysicalFieldVertices. No DeltaV0/Delta_g added.',
        spatial_route='Use returned raw complex Cphysical(q) and zero residual DeltaC. Do not feed raw complex C into the strictly-real PreparedReference API or silently discard its imaginary part.',
        invariant_inputs='nbar,K,original complex D,raw frames,FII,held/support terms,ions and geometry unchanged',
        field_preparation_status=data['field_preparation_status'],
        limits='initial join consistency only; no new full-slot H0 energy, trajectory or propagated noise')
    return arrays, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--oracle', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--columns', type=int, default=2, choices=(1, 2, 3, 4))
    parser.add_argument('--array-max-bytes', type=int, default=512*1024**2)
    parser.add_argument('--json-max-bytes', type=int, default=16*1024**2)
    args = parser.parse_args()
    enforce_runtime()
    registration_record, registration = authenticate_registration()
    implementation = implementation_bindings()
    oracle_record = binding(args.oracle)
    oracle = bound_json(oracle_record)
    if (oracle.get('status') != 'passed' or oracle['registration'] != registration_record
            or oracle['implementation_sources'] != implementation):
        raise ValueError('passed independent oracle bound to current implementation required')
    full.preload_parsers()
    ref_record, manifest_record, manifest, saved, original, grid = full.load_inputs(
        HERE/'full-common-reference-004/full-common-reference.json')
    background, background_receipt = full.restore_accepted(original)
    output = fresh_output(args.output)
    archives, used = {}, 0
    def save(name, arrays):
        nonlocal used
        size = sum(np.asarray(value).nbytes for value in arrays.values())
        # Reserve ZIP headers as well as array payload before writing.
        if used+size+4096*len(arrays) > args.array_max_bytes:
            raise ValueError('external witness output byte cap exceeded before allocation/write')
        path = output/(name+'.npz')
        with path.open('xb') as stream:
            np.savez(stream, **arrays)
        record = binding(path)
        used += int(record['bytes'])
        archives[name] = record
    report = dict(schema='P254-physical-spatial-material-v1', registration=registration_record,
        implementation_sources=implementation, oracle=oracle_record, runtime_context=runtime_context(),
        background_receipt=background_receipt, full_reference=ref_record, manifest=manifest_record,
        full_material_response=False, propagated_trajectory=False, h0_transport_supported=False,
        memory_scope='two full raw frames, endpoint physical factors and one legacy comparison bank; no PW-square matrices; bounded column FFT and64-column field GEMM')
    settled = next(r for r in registration['source_bindings'] if r['path'].endswith('joint-field-aggregate-001/result.json'))
    arrays, report['scalar_handoff'] = coefficient_join(settled, ref_record, manifest_record)
    save('scalar-handoff', arrays)
    del arrays
    legacy, frame, state, _, state_authority = joint.load_source((1, 1), manifest, saved, original, grid, background)
    h, probes, d, arrays, report['Gamma'] = same_basis_witness(legacy, frame, state, args.columns)
    report['Gamma']['state_authority'] = state_authority
    save('Gamma-witness', arrays)
    del arrays, legacy, frame, state, h, probes, d
    endpoints = []
    for ik in (2, 3):
        legacy, frame, state, _, state_authority = joint.load_source((1, ik), manifest, saved, original, grid, background)
        h, probes, d, arrays, item = same_basis_witness(legacy, frame, state, args.columns)
        item['state_authority'] = state_authority
        report['source-'+str(ik)] = item
        save('source-'+str(ik)+'-witness', arrays)
        endpoints.append((frame, h, probes, d))
        del arrays, legacy, state
    source, sh, right, d = endpoints[0]
    target, th, left, _ = endpoints[1]
    arrays, report['finite_q'] = rectangular_witness(background.inputs, source, target, sh, th, right, left, d)
    save('finite-q-witness', arrays)
    report.update(status='complete-operator-probes-and-scalar-handoff', archives=archives,
                  external_array_bytes=used, no_scientific_trajectory_claim=True)
    text = json.dumps(full.pack(report), indent=2, allow_nan=False)+'\n'
    if len(text.encode()) > args.json_max_bytes:
        raise ValueError('JSON output byte cap exceeded')
    (output/'result.json').write_text(text)
    print(json.dumps(dict(status=report['status'], result=binding(output/'result.json'))))


if __name__ == '__main__':
    main()
